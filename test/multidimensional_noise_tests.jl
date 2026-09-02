using RungeKutta: Tableau
using StochasticIntegrators
using Random
using Test

include("utils.jl")

# The test problem: a Kubo oscillator driven by **two** independent Wiener processes.
#
# Both diffusion columns are multiples of the drift, B^{·r} = ν_r v, so the solution is the
# deterministic rotation evaluated at the random time θ(t) = t + ν₁ W¹(t) + ν₂ W²(t), and ‖q‖ is
# conserved exactly along every path — for any increments at all, not only for a genuine Wiener
# path. That last point is what makes the invariant usable for the weak methods, whose increments
# are discrete random variables rather than a sampled Brownian path.
#
# Everything else in the suite is driven by one-dimensional noise, which leaves the branches that
# distinguish one noise dimension from another never executed: `qdiff3` in `WERK` and `WIRK`, and
# `qdiff2` in `WERK`, which additionally needs a non-zero ΔZ to be reached at all.

"Drift of the two-noise Kubo oscillator."
function kubo2_v(v, t, q, params)
    v[1] = q[2]
    v[2] = -q[1]
end

"Diffusion of the two-noise Kubo oscillator, both columns multiples of the drift."
function kubo2_B(B, t, q, params)
    B[1, 1] = +params.ν1 * q[2]
    B[2, 1] = -params.ν1 * q[1]
    B[1, 2] = +params.ν2 * q[2]
    B[2, 2] = -params.ν2 * q[1]
end

"The same diffusion restricted to a single noise dimension."
function kubo1_B(B, t, q, params)
    B[1, 1] = +params.ν1 * q[2]
    B[2, 1] = -params.ν1 * q[1]
end

const NDIMS = 2
const TSPAN = (0.0, 0.5)
const STEP = 0.01
const PARAMS2 = (ν1 = 0.1, ν2 = 0.15)
const Q0 = [0.5, 0.0]

"The two-dimensional-noise problem, driven by `process`."
function kubo2problem(process)
    SDEProblem(kubo2_v, kubo2_B, process, TSPAN, STEP, Q0; parameters = PARAMS2)
end

"A `GridProcess` of `m` prescribed noise columns for this problem's time grid, with ``ΔZ = 0``."
function gridnoise(m, seed)
    local nw = grid_length(TSPAN, STEP)
    GridProcess(randn(Xoshiro(seed), m, nw) .* sqrt(STEP), zeros(m, nw))
end

"""
The same, but with a **non-zero** second increment.

`gridnoise` prescribes ``ΔZ = 0``, which zeroes `WERK`'s `ydiff2` identically whatever the tableau
says, so the ``\\hat H^{(l)}`` stages are unreachable through it. The two-point law used here is
the one a weak method draws for ``\\tilde I``: ``±\\sqrt{Δt}`` with probability one half each.
"""
function gridnoise_nonzero_Z(m, seed)
    local nw = grid_length(TSPAN, STEP)
    local rng = Xoshiro(seed)
    GridProcess(randn(rng, m, nw) .* sqrt(STEP),
        rand(rng, (-1.0, 1.0), m, nw) .* sqrt(STEP))
end

@testset "$(rpad("Multidimensional noise", 80))" begin
    prob = kubo2problem(WienerProcess(NDIMS))
    @test noisedims(prob) == NDIMS

    # every method that accepts an SDEProblem must cope with m > 1
    for method in (BurrageE1(), StochasticGLRK(1), RoesslerRS1(), RoesslerRS2(), SRKw1())
        sol = integrate(prob, method)
        @test all(isfinite, sol.q[end])
        @test rel_energy_err(sol) < 1E-4
    end

    # the symplectic weak method conserves the quadratic invariant to round-off
    @test rel_energy_err(integrate(prob, SRKw1())) < 1E-13
end

@doc raw"""
Zero the `qdiff3` block of a weak tableau, leaving everything else alone.

`qdiff3` is the block applied to the noise dimensions *other than* the one whose stage family is
being formed, so it is exactly the off-diagonal coupling. Zeroing it gives a method that agrees
with the original whenever there is only one noise dimension and differs as soon as there is more
than one — which is what the tests below assert.
"""
function without_offdiagonal(tab::TableauWIRK)
    local z = Tableau(:no_offdiagonal, 0, zero(Matrix(tab.qdiff3.a)),
        Vector(tab.qdiff3.b), Vector(tab.qdiff3.c))
    TableauWIRK(:no_offdiagonal, tab.qdrift0, tab.qdrift1, tab.qdiff0, tab.qdiff1, z)
end

function without_offdiagonal(tab::TableauWERK)
    local z = Tableau(:no_offdiagonal, 0, zero(Matrix(tab.qdiff3.a)),
        Vector(tab.qdiff3.b), Vector(tab.qdiff3.c))
    TableauWERK(:no_offdiagonal, tab.qdrift0, tab.qdrift1, tab.qdrift2,
        tab.qdiff0, tab.qdiff1, tab.qdiff2, z)
end

@doc raw"""
Zero the `qdiff2` block of a `WERK` tableau, leaving everything else alone.

`qdiff2` is the only block that multiplies ``\Delta Z``, and it appears solely in the
``\hat H^{(l)}`` stages, whose terms carry ``\hat I_{r,l}`` and vanish for ``l = r``. Zeroing it
therefore gives a method that agrees with the original unless *both* noise dimensions and a
non-zero ``\Delta Z`` are present — which is what the tests below assert. `WIRK` has no
counterpart: it stores ``\Delta Z`` in its cache and never reads it.
"""
function without_deltaZ(tab::TableauWERK)
    local z = Tableau(:no_deltaZ, 0, zero(Matrix(tab.qdiff2.a)),
        Vector(tab.qdiff2.b), Vector(tab.qdiff2.c))
    TableauWERK(:no_deltaZ, tab.qdrift0, tab.qdrift1, tab.qdrift2,
        tab.qdiff0, tab.qdiff1, z, tab.qdiff3)
end

# These are the tests that give the coverage its meaning. Running a weak method on a
# two-dimensional problem exercises the off-diagonal branch, but on its own it does not show that
# the branch *does* anything — a scheme that ignored `qdiff3` entirely would still produce finite,
# energy-conserving output. So each method is run twice on one prescribed path, once with its own
# tableau and once with `qdiff3` zeroed:
#
#   * with two noise dimensions the two must **disagree** — the off-diagonal block is reached and
#     its value matters;
#   * with one noise dimension they must agree **exactly** — there is no other dimension to couple
#     to, so the branch is never taken.
#
# Together those pin down both that the branch runs and when.

@testset "$(rpad("Off-diagonal noise coupling is reached and matters", 80))" begin
    local gp2 = gridnoise(NDIMS, 4242)
    local prob2 = kubo2problem(gp2)

    for (name, tab, wrap) in (("WIRK/SRKw1", TableauSRKw1(), WIRK),
        ("WERK/RoesslerRS1", TableauRoesslerRS1(), WERK),
        ("WERK/RoesslerRS2", TableauRoesslerRS2(), WERK))
        @testset "$name" begin
            @test any(!iszero, tab.qdiff3.a)      # otherwise the comparison is vacuous

            local full = integrate(prob2, wrap(tab))
            local zeroed = integrate(prob2, wrap(without_offdiagonal(tab)))

            @test full.q[end] != zeroed.q[end]
            @test maximum(abs.(full.q[end] .- zeroed.q[end])) > 1E-10
        end
    end
end

@testset "$(rpad("Off-diagonal coupling is inert for one-dimensional noise", 80))" begin
    # The same comparison with m = 1: with no second noise dimension the off-diagonal block is
    # never reached, so zeroing it must change nothing at all.
    local nw = grid_length(TSPAN, STEP)
    local gp1 = GridProcess(randn(Xoshiro(99), 1, nw) .* sqrt(STEP), zeros(1, nw))

    local prob1 = SDEProblem(kubo2_v, kubo1_B, gp1, TSPAN, STEP, Q0; parameters = PARAMS2)
    @test noisedims(prob1) == 1

    for (tab, wrap) in ((TableauSRKw1(), WIRK), (TableauRoesslerRS1(), WERK))
        @test integrate(prob1, wrap(tab)).q[end] ==
              integrate(prob1, wrap(without_offdiagonal(tab))).q[end]
    end
end

# `qdiff3` is not the only block that needs more than one noise dimension. `WERK`'s Ĥ^(l) stages
# carry the terms in Î_{r,l}, which are gated twice over: they vanish for l = r, and they are
# multiplied by ΔZ. The tests above prescribe ΔZ = 0, so that block contributes nothing to them —
# a `qdiff2` of any value at all would pass.
#
# The Kubo oscillator cannot close that gap either, and the reason is worth recording. Both of its
# diffusion columns are multiples of the drift, so the two diffusion vector fields commute, and for
# commutative noise the Î_{r,l} contributions cancel: zeroing `qdiff2` moves the solution by 5.6e-17
# — round-off — however the block is weighted. The commutativity is exactly what makes ‖q‖ an
# invariant for any increments, so the property the tests above rely on for their assertion is the
# one that blinds them to this block. A non-commutative diffusion is required.

"""
Diffusion whose two columns do not commute: ``[B^{·1}, B^{·2}] = ν₁ν₂ (-q₁, q₂) ≠ 0``.

Nothing is conserved along this problem, which is fine — the assertions below compare two runs
against each other rather than against an invariant.
"""
function noncommuting_B(B, t, q, params)
    B[1, 1] = params.ν1 * q[2]
    B[2, 1] = 0.0
    B[1, 2] = 0.0
    B[2, 2] = params.ν2 * q[1]
end

"The first column of `noncommuting_B` alone."
function noncommuting_B1(B, t, q, params)
    B[1, 1] = params.ν1 * q[2]
    B[2, 1] = 0.0
end

@testset "$(rpad("The ΔZ coupling of WERK is reached and matters", 80))" begin
    local probz = SDEProblem(kubo2_v, noncommuting_B, gridnoise_nonzero_Z(NDIMS, 8484),
        TSPAN, STEP, Q0; parameters = PARAMS2)

    for (name, tab) in (("WERK/RoesslerRS1", TableauRoesslerRS1()),
        ("WERK/RoesslerRS2", TableauRoesslerRS2()))
        @testset "$name" begin
            @test any(!iszero, tab.qdiff2.a)      # otherwise the comparison is vacuous

            local full = integrate(probz, WERK(tab))
            local zeroed = integrate(probz, WERK(without_deltaZ(tab)))

            @test full.q[end] != zeroed.q[end]
            @test maximum(abs.(full.q[end] .- zeroed.q[end])) > 1E-10
        end
    end

    # one noise dimension: Î_{r,l} has no second dimension to couple to, so the block is inert
    local prob1 = SDEProblem(kubo2_v, noncommuting_B1, gridnoise_nonzero_Z(1, 555),
        TSPAN, STEP, Q0; parameters = PARAMS2)

    for tab in (TableauRoesslerRS1(), TableauRoesslerRS2())
        @test integrate(prob1, WERK(tab)).q[end] ==
              integrate(prob1, WERK(without_deltaZ(tab))).q[end]
    end

    # two noise dimensions but ΔZ ≡ 0: also inert, which is why the tests above could not see it
    local prob0 = SDEProblem(kubo2_v, noncommuting_B, gridnoise(NDIMS, 8484),
        TSPAN, STEP, Q0; parameters = PARAMS2)

    for tab in (TableauRoesslerRS1(), TableauRoesslerRS2())
        @test integrate(prob0, WERK(tab)).q[end] ==
              integrate(prob0, WERK(without_deltaZ(tab))).q[end]
    end
end
