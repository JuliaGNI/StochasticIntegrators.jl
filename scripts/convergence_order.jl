#
# Measured mean-square convergence order of the strong methods on the Kubo oscillator.
#
# Companion to `scripts/tableau_conditions.jl`. That script checks the coefficients; this one
# checks the integrator that consumes them, by measuring the order the theory predicts.
#
# The reference is
#
#   M. Kraus and T. M. Tyranowski,
#   "Variational integrators for stochastic dissipative Hamiltonian systems",
#   IMA Journal of Numerical Analysis (2021),
#
# Theorem (Mean-square convergence theorem) and its Corollary. Two caveats about what that theory
# says, both of which this script depends on:
#
#   * Order 1.0 requires the noise to be **commutative**, Eq. (Commutation conditions); otherwise
#     the theorem gives only 0.5. The Kubo oscillator is driven by a one-dimensional Wiener
#     process, for which commutativity is automatic — that is the Corollary — so 1.0 is the
#     expected order here and would not be for a general multi-dimensional problem.
#
#   * These schemes use only the increments ΔW, which caps them at mean-square order 1.0 however
#     many stages they have. Order 1.5 would additionally need the iterated integrals ΔZ, which
#     among the methods here only the explicit Burrage schemes carry.
#
# The paper measures no convergence orders itself — its orders are theoretical, and the only
# error figures it reports for the Kubo oscillator are Monte-Carlo sampling errors. So this is a
# new measurement against its theorems rather than a reproduction of a published table.
#
# Method. The Kubo oscillator has a closed-form solution: its diffusion is proportional to its
# drift, so the solution is the deterministic one evaluated at the random time θ(t) = t + νW(t).
# Driving the integrator with a `GridProcess` of prescribed increments therefore gives a
# *pathwise* exact reference — no fine-grid reference solution, and no Monte-Carlo error in the
# reference itself. The mean-square error
#
#     e(Δt) = sqrt( E |q(T) - q_N|² )
#
# is estimated over a sample of paths, and its order is the slope of log e against log Δt.
#
# Usage:  julia --project=. scripts/convergence_order.jl

using GeometricProblems.KuboOscillator
using Printf
using Random
using Statistics
using StochasticIntegrators
using Test

import GeometricProblems.KuboOscillator as Kubo

const T_FINAL = 1.0
const NPATHS = 400
const REFINEMENTS = [2^k for k in 3:7]      # 8 … 128 steps, i.e. Δt from 1/8 to 1/128
const PARAMS = (ν = 0.5,)                   # a large noise intensity, to expose the noise order
const Q0 = [1.0]
const P0 = [0.5]

@doc raw"""
Build a `GridProcess` on `nsteps` steps of size `Δt` by coarsening one fixed fine-grid Brownian
path, so that every refinement sees *the same* realisation of the driving process.

This is what makes the comparison a strong-convergence measurement rather than a comparison of
unrelated random runs: the error at each `Δt` is the pathwise error along one common path.

Both families of increments are formed, not just `ΔW`:
```math
\Delta W_n = \int_{t_n}^{t_{n+1}} dW , \qquad
\Delta Z_n = \int_{t_n}^{t_{n+1}} \bigl( W(s) - W(t_n) \bigr) \, ds ,
```
the second by the trapezoidal rule over the fine subintervals. Leaving `ΔZ` at zero would silently
cost a scheme that uses it — `BurrageE1` and `BurrageG5` carry a second diffusion tableau against
exactly these terms — the very order those terms buy, and it would look like a convergence
failure rather than a badly posed measurement.
"""
function refine_path(fine::AbstractVector, nfine::Int, nsteps::Int, Δt::Float64)
    stride = nfine ÷ nsteps
    δ = Δt / stride

    ΔW = zeros(1, nsteps)
    ΔZ = zeros(1, nsteps)

    for n in 1:nsteps
        local w = 0.0          # W(s) - W(tₙ) at the right end of the current fine subinterval
        local z = 0.0
        for j in 1:stride
            local dw = fine[(n - 1) * stride + j]
            z += (w + (w + dw)) / 2 * δ      # trapezoid over this subinterval
            w += dw
        end
        ΔW[1, n] = w
        ΔZ[1, n] = z
    end

    (GridProcess(ΔW, ΔZ), sum(@view fine[1:(nsteps * stride)]))
end

"""
The Kubo oscillator as an SDE, with `q = (q₁, q₂)` playing the role of `(q, p)`, driven by the
prescribed increments `gp`.
"""
function sde_problem(gp, Δt)
    SDEProblem(Kubo.kubo_oscillator_sde_v, Kubo.kubo_oscillator_sde_B, gp,
        (0.0, T_FINAL), Δt, [Q0[begin], P0[begin]]; parameters = PARAMS)
end

"The same oscillator as a partitioned SDE."
function psde_problem(gp, Δt)
    PSDEProblem(Kubo.kubo_oscillator_psde_v, Kubo.kubo_oscillator_psde_f,
        Kubo.kubo_oscillator_psde_B, Kubo.kubo_oscillator_psde_G,
        gp, (0.0, T_FINAL), Δt, Q0, P0; parameters = PARAMS)
end

sde_endpoint(sol) = (sol.q[end][1], sol.q[end][2])
psde_endpoint(sol) = (sol.q[end][begin], sol.p[end][begin])

"Mean-square error of `method` on the Kubo oscillator at `nsteps` steps."
function mean_square_error(method_builder, problem_builder, endpoint, nsteps::Int, rng)
    local Δt = T_FINAL / nsteps
    local nfine = maximum(REFINEMENTS)
    local Δt_fine = T_FINAL / nfine

    local errs = zeros(NPATHS)

    for k in 1:NPATHS
        fine = randn(rng, nfine) .* sqrt(Δt_fine)
        gp, W_T = refine_path(fine, nfine, nsteps, Δt)

        sol = integrate(problem_builder(gp, Δt), method_builder())
        q, p = endpoint(sol)

        qex, pex = exact_solution(T_FINAL, W_T, Q0[begin], P0[begin], PARAMS)
        errs[k] = (q - qex)^2 + (p - pex)^2
    end

    sqrt(mean(errs))
end

"Least-squares slope of log(error) against log(Δt)."
function fitted_order(steps, errors)
    x = log.(T_FINAL ./ steps)
    y = log.(errors)
    x̄ = mean(x)
    ȳ = mean(y)
    sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end

function run(name, builder, problem_builder, endpoint)
    errors = Float64[]
    for n in REFINEMENTS
        push!(errors,
            mean_square_error(builder, problem_builder, endpoint, n, Xoshiro(20260902)))
    end

    @printf("\n%s\n", name)
    @printf("  %8s  %12s  %10s\n", "steps", "rms error", "local order")
    for (k, n) in enumerate(REFINEMENTS)
        if k == 1
            @printf("  %8d  %12.4e  %10s\n", n, errors[k], "-")
        else
            local p = log(errors[k - 1] / errors[k]) /
                      log(REFINEMENTS[k] / REFINEMENTS[k - 1])
            @printf("  %8d  %12.4e  %10.3f\n", n, errors[k], p)
        end
    end
    local order = fitted_order(REFINEMENTS, errors)
    @printf("  fitted order: %.3f\n", order)

    (order, errors)
end

const METHODS = [
    ("StochasticGLRK(1)", () -> StochasticGLRK(1), sde_problem, sde_endpoint),
    ("StochasticDIRK(0.5)", () -> StochasticDIRK(0.5), sde_problem, sde_endpoint),
    ("BurrageE1", () -> BurrageE1(), sde_problem, sde_endpoint),
    ("StochasticStoermerVerlet", () -> StochasticStoermerVerlet(), psde_problem,
        psde_endpoint)
]

function main()
    results = Dict{String, Tuple{Float64, Vector{Float64}}}()

    for (name, builder, pbuilder, endpoint) in METHODS
        results[name] = run(name, builder, pbuilder, endpoint)
    end

    println("\n", "="^78, "\n")

    @testset "Mean-square convergence order on the Kubo oscillator" begin
        for (name, _, _, _) in METHODS
            order, errors = results[name]
            @testset "$name" begin
                # Order 1.0 is what the theorem gives for one-dimensional (hence commutative)
                # noise. The band is wide enough to absorb the Monte-Carlo error of a 400-path
                # sample and narrow enough to separate 1.0 from 0.5, which is the distinction
                # that matters.
                @test 0.85 < order < 1.25

                # and the error must actually decrease under refinement
                @test issorted(errors; rev = true)
            end
        end
    end

    #
    # An open question, recorded rather than asserted.
    #
    # `BurrageE1` carries a second diffusion tableau against ΔZ, and its docstring — following
    # Burrage & Burrage (2000) — claims strong order 1.5 for one-dimensional Brownian motion.
    # Measured here with iterated integrals formed from the driving path it comes out at ≈ 1.05,
    # and its error is within a few percent of `StochasticGLRK(1)` at every refinement, which
    # says both are limited by the same O(Δt) term rather than that one is a half-order better.
    #
    # This is not a regression introduced by the rewrite: the arithmetic of the ΔZ terms is
    # carried over unchanged from the previous implementation, and the tableau coefficients are
    # untouched. It is either a property of this test problem that the published claim does not
    # cover, or a defect in the ΔZ terms that predates this work. It is left as an open issue in
    # CHANGELOG.md rather than papered over with a tolerance that would accept either answer.
    #
    println("\nNote: BurrageE1 measures order ", @sprintf("%.3f", results["BurrageE1"][1]),
        " against a documented claim of 1.5 — see CHANGELOG.md, Open Issues.")
end

main()
