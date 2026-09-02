using RungeKutta: Tableau
using StochasticIntegrators
using Random
using Test

include("utils.jl")

@doc raw"""
Kubo oscillator driven by **two** independent Wiener processes.

Both diffusion columns are multiples of the drift, ``B^{\cdot r} = \nu_r \, v``, so the solution
is the deterministic rotation evaluated at the random time
``\theta(t) = t + \nu_1 W^1(t) + \nu_2 W^2(t)`` and ``\lVert q \rVert`` is conserved **exactly**
along every path — for any increments at all, not only for a genuine Wiener path. That last point
is what makes the invariant usable for the weak methods, whose increments are discrete random
variables rather than a sampled Brownian path.

Everything else in the suite is driven by one-dimensional noise, which leaves the branches that
distinguish one noise dimension from another — `qdiff3` in [`WERK`](@ref) and [`WIRK`](@ref) —
never executed.
"""
function kubo2_v(v, t, q, params)
    v[1] = q[2]
    v[2] = -q[1]
end

function kubo2_B(B, t, q, params)
    B[1, 1] = +params.ν1 * q[2]
    B[2, 1] = -params.ν1 * q[1]
    B[1, 2] = +params.ν2 * q[2]
    B[2, 2] = -params.ν2 * q[1]
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

"A `GridProcess` of `m` prescribed noise columns for this problem's time grid."
function gridnoise(m, seed)
    local nw = grid_length(TSPAN, STEP)
    GridProcess(randn(Xoshiro(seed), m, nw) .* sqrt(STEP), zeros(m, nw))
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

    function kubo1_B(B, t, q, params)
        B[1, 1] = +params.ν1 * q[2]
        B[2, 1] = -params.ν1 * q[1]
    end

    local prob1 = SDEProblem(kubo2_v, kubo1_B, gp1, TSPAN, STEP, Q0; parameters = PARAMS2)
    @test noisedims(prob1) == 1

    for (tab, wrap) in ((TableauSRKw1(), WIRK), (TableauRoesslerRS1(), WERK))
        @test integrate(prob1, wrap(tab)).q[end] ==
              integrate(prob1, wrap(without_offdiagonal(tab))).q[end]
    end
end
