@doc raw"""
    StochasticMethod <: GeometricMethod

Root of the stochastic method hierarchy.

`GeometricIntegratorsBase` splits `GeometricMethod` into `DeterministicMethod` and nothing else;
the stochastic branch lives here. A stochastic method differs from a deterministic one in that
each step consumes a realisation of the driving noise, drawn according to
[`convergence`](@ref) from the `AbstractStochasticProcess` carried by
the problem.
"""
abstract type StochasticMethod <: GeometricMethod end

"Stochastic method for a `SDEProblem`."
abstract type SDEMethod <: StochasticMethod end

"Stochastic method for a `PSDEProblem`."
abstract type PSDEMethod <: StochasticMethod end

"Stochastic method for a `SPSDEProblem`."
abstract type SPSDEMethod <: StochasticMethod end

"""
    issdemethod(method)

Whether `method` integrates an `SDEProblem`. Accepts an instance or a type, and answers `false`
for any other `GeometricMethod`, so it can be asked of a method of unknown kind.
"""
issdemethod(::Union{SDEMethod, Type{<:SDEMethod}}) = true
issdemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false

"""
    ispsdemethod(method)

Whether `method` integrates a `PSDEProblem`. See [`issdemethod`](@ref).
"""
ispsdemethod(::Union{PSDEMethod, Type{<:PSDEMethod}}) = true
ispsdemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false

"""
    isspsdemethod(method)

Whether `method` integrates an `SPSDEProblem`, the split form in which the forcing is given as two
separate terms. See [`issdemethod`](@ref).
"""
isspsdemethod(::Union{SPSDEMethod, Type{<:SPSDEMethod}}) = true
isspsdemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false

@doc raw"""
    convergence(method)

Whether a method needs increments that approximate the driving Wiener process in the **strong**
(mean-square) or the **weak** sense — `:strong` or `:weak`.

This is a property of the scheme, not a choice the caller makes. A strong scheme is derived
against the true increments ``\Delta W^r = \int_{t_n}^{t_{n+1}} dW^r``, and needs them drawn from
the correct Gaussian law, because it is asked to track an individual sample path. A weak scheme is
derived only against the *moments* of the increments, and deliberately replaces them with a cheap
discrete random variable that matches those moments — which is valid for computing expectations
and wrong for tracking a path.

Feeding a scheme the other kind of increment produces plausible-looking numbers of the wrong
accuracy, silently, which is why this is derived from the method rather than left as an option.
See [`sample_noise!`](@ref) for the two laws.

```jldoctest
julia> convergence(BurrageE1())
:strong

julia> convergence(StochasticStoermerVerlet())
:strong

julia> convergence(SRKw1())
:weak
```
"""
function convergence end

"Number of internal stages of a stochastic Runge-Kutta method."
@inline nstages(method::StochasticMethod) = tableau(method).s
@inline eachstage(method::StochasticMethod) = Base.OneTo(nstages(method))

GeometricBase.tableau(method::StochasticMethod) = method.tableau
GeometricBase.name(method::StochasticMethod) = tableau(method).name

"Random number generator a stochastic method draws its increments from."
rng(method::StochasticMethod) = method.rng

@doc raw"""
    default_extrapolation(::StochasticMethod)

`NoExtrapolation()` — a stochastic problem has no history to extrapolate to.

`SolutionStep` keeps states at ``t_{-1}, t_{-2}, \ldots``, and for a deterministic problem the
framework fills them by integrating backwards from ``t_0``. That is meaningless here: the past of
a sample path is not determined by its present, and running an SDE backwards is not the inverse of
running it forwards. `NoExtrapolation` copies the initial state into those slots instead, which is
harmless because none of the schemes in this package read the history — every one of them is a
one-step method.
"""
GeometricIntegratorsBase.default_extrapolation(::StochasticMethod) = NoExtrapolation()

@doc raw"""
    truncation(method)

The integer ``K`` parameterising the truncation of the Wiener increments before an implicit solve,
or `0` for no truncation. The bound ``A`` itself also needs the time step, and is computed by
[`truncation_bound`](@ref) — which is internal, so reach it as
`StochasticIntegrators.truncation_bound`.

For an implicit scheme the nonlinear system need not be solvable for arbitrarily large increments,
and a Gaussian increment is unbounded. Milstein & Tretyakov's remedy is to clip ``\Delta W`` at
``A = \sqrt{2 K \Delta t \, \lvert \log \Delta t \rvert}``, which leaves the mean-square order
intact because the discarded tail is of higher order than the scheme. `K = 0` disables it, and is
the default for every method.

```jldoctest
julia> truncation(StochasticGLRK(1))
0

julia> truncation(StochasticGLRK(1; K = 2))
2
```
"""
truncation(::StochasticMethod) = 0

"Milstein-Tretyakov truncation bound from the parameter `K` and the time step."
truncation_bound(K::Integer, Δt) = K == 0 ? zero(Δt) : sqrt(2 * K * Δt * abs(log(Δt)))

function Base.show(io::IO, method::StochasticMethod)
    print(io, nameof(typeof(method)), " method with tableau ", name(method))
end
