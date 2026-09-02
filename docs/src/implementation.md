```@meta
CurrentModule = StochasticIntegrators
```

# Implementation

How the package fits into `GeometricIntegratorsBase`, and what to write to add a method of your
own.

## The step

`GeometricIntegratorsBase` drives the integration. Per step it does roughly

```julia
reset!(solstep, t)                              # shift history, advance t
initial_guess!(current, history, params, int)   # framework hook, unused here
integrate_step!(current, history, params, int)  # the method's own work
copy_internal_variables!(solstep, cache(int))
enforce_periodicity!(solstep)
compute_vectorfields!(current(solstep), problem(int))
```

so a method contributes an `integrate_step!` and the framework supplies everything else — the
loop, the solution storage, the NaN and non-convergence handling, ensembles.

Two conventions matter when reading the step implementations:

* **On entry `sol.t` already holds ``t_{n+1}`` while `sol.q` still holds ``q_n``.** So `sol.q` is
  the previous solution, and the previous time is `t̄ = sol.t - timestep(int)`.
* **`sol.q` is a `StateWithError`**, so the increment goes in with `GeometricBase.add!`, which is
  the operation that performs the compensated summation. `sol.q .+= Δq` does **not**: broadcasting
  falls through to the generic `AbstractArray` path, writes through `setindex!` and never touches
  the error field. Each mathematical increment is therefore accumulated in full into a scratch
  vector before it is added, rather than added term by term; splitting one increment over several
  `add!` calls would compensate each piece separately for no benefit. A step still makes more than
  one `add!` — the weights `b` and the correction weights `b̂` that `RungeKutta.Tableau` carries are
  applied separately, which is a higher-precision splitting rather than one increment broken up.

A stochastic step adds one thing at the front: drawing the increments for the step.

```@docs
sample_noise!(::GeometricIntegrator{<:StochasticMethod}, ::Any)
```

## Where the increments live

The increments are stored in the integrator's cache, and always in the cache belonging to the
problem's own data type — `cache(int)`, never `cache(int, ST)`.

This matters because a Newton solve propagates dual numbers, so `CacheDict` holds a separate cache
per element type, and a `residual!` evaluated in dual arithmetic sees a *different* cache from the
one `integrate_step!` sampled into. The increments are ordinary real numbers that multiply dual
quantities, and every arithmetic path within one step must see the same draw, so they are read
from the real-valued cache throughout:

```julia
function residual!(b, sol, params, int::GeometricIntegrator{<:SIRK})
    c  = cache(int, ST)      # ST-typed working arrays
    ΔW = cache(int).ΔW       # the increments, shared and real
    ...
end
```

The random number generator lives on the method, since the cache is per-element-type and a
per-cache generator would split the stream.

## History and extrapolation

A `SolutionStep` keeps states at ``t_{-1}, t_{-2}, \ldots`` that no initial condition provides,
and the framework fills them by running the problem *backwards* from ``t_0``. For a deterministic
problem that is well defined; for a stochastic one it is not, because the past of a sample path is
not a function of its present.

Stochastic methods therefore declare

```julia
GeometricIntegratorsBase.default_extrapolation(::StochasticMethod) = NoExtrapolation()
```

which copies the initial state into the history slots instead. This is harmless because every
method here is a one-step method and none reads the history.

## Caches

```@docs
StochasticIntegratorCache
SDEIntegratorCache
PSDEIntegratorCache
```

Each family's cache holds `ΔW` and `ΔZ` for the current step, scratch for the final update, and
its internal stage storage. Helpers:

```@docs
create_internal_stage_vector
create_internal_stage_matrix
```

## Adding a method

To add a family for, say, an `SDEProblem`:

1. **Define the tableau and the method.** The method subtypes `SDEMethod`, `PSDEMethod` or
   `SPSDEMethod`, and carries its tableau and a generator:

   ```julia
   struct MyMethod{TT, RNG} <: SDEMethod
       tableau::TableauMine{TT}
       rng::RNG
   end

   convergence(::MyMethod) = :strong
   ```

2. **Define a cache** subtyping `SDEIntegratorCache` with at least `ΔW` and `ΔZ` fields, and the
   `Cache`/`CacheType` constructors:

   ```julia
   function GeometricIntegratorsBase.Cache{ST}(problem::AbstractProblemSDE,
           method::MyMethod; kwargs...) where {ST}
       MyCache{ST, length(vec(initial_conditions(problem).q)),
               noisedims(problem), nstages(method)}(; kwargs...)
   end
   ```

   Note `noisedims(problem)`: the noise dimension comes from the problem, via its noise process.

3. **For an implicit method**, add `nlsolution(cache) = cache.x`, `default_solver`, `solversize`,
   and the `components!` / `residual!` pair.

4. **Define `integrate_step!`**, beginning with `sample_noise!(int, sol)`:

   ```julia
   function GeometricIntegratorsBase.integrate_step!(sol, history, params,
           int::GeometricIntegrator{<:MyMethod, <:AbstractProblemSDE})
       sample_noise!(int, sol)
       ...
       update_solution!(sol.q, ...)
   end
   ```

If a method needs an initial guess that depends on the increments — as [`SIRK`](@ref) and
[`SIPRK`](@ref) do — it must **not** use the framework's `initial_guess!` hook, which runs before
`integrate_step!` and therefore before the increments exist. Call a private predictor from inside
`integrate_step!` instead; `stage_initial_guess!` is the name used here.

## Final updates

```@docs
update_solution!
update_solution_weak!
```

These carry the numerical core of every scheme and are shared between families — [`WIRK`](@ref)
uses the same final update as [`SIRK`](@ref), for instance. The weak update carries its own name
because it is a different formula rather than an overload: it takes two families of diffusion
stages and a `√Δt` term that no strong update has. Each is called twice by its caller,
once with the tableau weights `b` and once with the correction weights `b̂` that
`RungeKutta.Tableau` carries; that pair is a higher-precision splitting of the same weights and
keeping the two additions separate is deliberate.

## Upstream requirements

The package needs three things of its dependencies that a purely deterministic integrator would
not, all of them added for this purpose:

* `GeometricBase` declares the interface stubs `noise` and `noisedims`.
* `GeometricEquations` provides the concrete `WienerProcess` and `GridProcess`, and the accessors
  on the stochastic equations and problems.
* `GeometricIntegratorsBase` dispatches `default_extrapolation` on the method, so that a family
  for which running the problem backwards is meaningless can decline it.
