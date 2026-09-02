#
# Per-step allocations of the six method families on the Kubo problems.
#
# This is the check behind the allocation table in `CHANGELOG.md`, *Performance*. The claim it
# establishes is that `integrate_step!` allocates **nothing** for any of the six families — the
# tableau fields carry one type parameter per Butcher tableau, so a stage loop infers
# `tab.qdrift.a[i,j]` concretely and does not materialise the `SMatrix` on the heap.
#
# Two measurements, because neither is sufficient on its own:
#
#   A. `@allocated` around `integrate_step!`, which is the quantity the claim is about. It is
#      taken after a warm-up, since the first call in a process measures compilation.
#
#   B. `d(alloc)/d(nsteps)` across a full `integrate` run — a 100- and a 1100-step problem, the
#      difference divided by the 1000 extra steps. This catches a per-step cost that sits in the
#      framework's loop rather than in the method, which A cannot see. It does not go to zero,
#      and is not meant to: it measures the solution storage, is identical for every method of
#      the same problem type, and is outside this package.
#
# B minus A is the framework's share. That it comes out the same for every method — and the same
# on the code before the change, where A was in the tens of kilobytes — is what says the residual
# is not method-dependent.
#
# Run with `--check-bounds=auto`. `Pkg.test()` defaults to `--check-bounds=yes`, which changes
# allocation behaviour; a figure measured under it does not establish the claim.
#
# Usage:  julia --project=scripts --check-bounds=auto scripts/step_allocations.jl
#

using GeometricProblems.KuboOscillator
using Printf
using StochasticIntegrators
using Test

import GeometricIntegratorsBase as GIB
import GeometricProblems.KuboOscillator as Kubo

# `problem` is the constructor rather than a built problem, so B can rebuild it at another length.
const CASES = [
    (name = "SERK", problem = sdeproblem, method = BurrageE1()),
    (name = "SIRK", problem = sdeproblem, method = StochasticGLRK(1)),
    (name = "WERK", problem = sdeproblem, method = RoesslerRS1()),
    (name = "WIRK", problem = sdeproblem, method = SRKw1()),
    (name = "SIPRK", problem = psdeproblem, method = StochasticStoermerVerlet()),
    (name = "SISPRK", problem = spsdeproblem, method = ModifiedStochasticStoermerVerlet())
]

"Median of `n` `@allocated` measurements of one `integrate_step!`, after `warmup` steps."
function step_allocations(problem, method; warmup = 10, n = 30)
    int = GeometricIntegrator(problem, method)
    solstep = GIB.solutionstep(int, GIB.Solution(problem)[0])
    current, history = GIB.current(solstep), GIB.history(solstep)
    params = GIB.parameters(solstep)

    step!() = GIB.integrate_step!(current, history, params, int)

    for _ in 1:warmup
        step!()
    end

    sort([@allocated(step!()) for _ in 1:n])[(n + 1) ÷ 2]
end

"Marginal allocation per step of a full `integrate`, from runs of `n₁` and `n₂` steps."
function marginal_allocations(constructor, method; n₁ = 100, n₂ = 1100, samples = 5)
    p₁ = constructor(; timespan = (0.0, Kubo.Δt * n₁))
    p₂ = constructor(; timespan = (0.0, Kubo.Δt * n₂))

    integrate(p₁, method)
    integrate(p₂, method)

    a₁ = minimum([@allocated(integrate(p₁, method)) for _ in 1:samples])
    a₂ = minimum([@allocated(integrate(p₂, method)) for _ in 1:samples])

    (a₂ - a₁) / (n₂ - n₁)
end

@printf("julia %s, check_bounds=%d\n\n", VERSION, Base.JLOptions().check_bounds)
@printf("%-8s %14s %16s %16s\n", "method", "per step (B)", "per run step (B)",
    "framework (B)")

results = map(CASES) do case
    per_step = step_allocations(case.problem(), case.method)
    per_run = marginal_allocations(case.problem, case.method)
    @printf("%-8s %14d %16.2f %16.2f\n", case.name, per_step, per_run, per_run - per_step)
    (case.name, per_step, per_run - per_step)
end

@testset "step allocations" begin
    @testset "$name" for (name, per_step, _) in results
        @test per_step == 0
    end

    # The framework's share is a property of the problem type, not of the method: the four SDE
    # methods must agree with each other, and the two that carry a momentum with each other.
    sde = [f for (n, _, f) in results if n in ("SERK", "SIRK", "WERK", "WIRK")]
    momentum = [f for (n, _, f) in results if n in ("SIPRK", "SISPRK")]

    @test all(≈(first(sde)), sde)
    @test all(≈(first(momentum)), momentum)
end
