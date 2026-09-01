@doc raw"""
    StochasticIntegratorCache{DT} <: IntegratorCache{DT}

Cache of a stochastic integrator. Beyond the usual internal-stage storage it holds the increments
`ΔW` and `ΔZ` of the current time step, drawn once per step by [`sample_noise!`](@ref).
"""
abstract type StochasticIntegratorCache{DT} <: IntegratorCache{DT} end

"Cache of an integrator for a `SDEProblem`."
abstract type SDEIntegratorCache{DT} <: StochasticIntegratorCache{DT} end

"Cache of an integrator for a `PSDEProblem` or `SPSDEProblem`."
abstract type PSDEIntegratorCache{DT} <: StochasticIntegratorCache{DT} end

"The increments `(ΔW, ΔZ)` of the current time step."
increments(cache::StochasticIntegratorCache) = (cache.ΔW, cache.ΔZ)

"`S` internal stage vectors of length `D`."
create_internal_stage_vector(DT, D, S) = [zeros(DT, D) for _ in 1:S]

"`S` internal stage matrices of size `D × M`."
create_internal_stage_matrix(DT, D, M, S) = [zeros(DT, D, M) for _ in 1:S]

@doc raw"""
    sample_noise!(int, sol)

Draw the increments for the step that `int` is about to take and store them in its cache.

Called at the top of every `integrate_step!`. The increments go into the cache belonging to the
problem's own data type — `cache(int)`, never `cache(int, ST)` — because they are ordinary real
numbers rather than the dual numbers a Newton solve propagates, and because every arithmetic path
within one step must see the *same* draw. A `residual!` evaluated in dual arithmetic therefore
reads its increments from `cache(int)` too.

The step index handed to the process is recovered from the current time. On entry `sol.t` already
holds ``t_{n+1}``, so with a constant time step ``n = (t_{n+1} - t_0) / \Delta t`` is the
one-based index of the step being taken.
"""
function sample_noise!(int::GeometricIntegrator{<:StochasticMethod}, sol)
    local Δt = timestep(int)
    local c = cache(int)
    local n = round(Int, (sol.t - initialtime(problem(int))) / Δt)

    sample_noise!(c.ΔW, c.ΔZ, noise(problem(int)),
        Val(convergence(method(int))), Δt, n, rng(method(int)))

    truncate_increments!(c.ΔW, truncation_bound(truncation(method(int)), Δt))

    c
end

@doc raw"""
    noise_increments(int, ST)

The increments of the current step, for use inside code specialised on the solver's element type
`ST`. Always the real-valued ones from `cache(int)`; see [`sample_noise!`](@ref).
"""
@inline noise_increments(int::GeometricIntegrator{<:StochasticMethod}) = increments(cache(int))
