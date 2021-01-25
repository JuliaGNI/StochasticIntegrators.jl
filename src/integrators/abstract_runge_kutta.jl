
abstract type StochasticIntegratorRK{dType, tType, D, M, S} <: SDEIntegrator{dType, tType} end
abstract type StochasticIntegratorPRK{dType, tType, D, M, S} <: PSDEIntegrator{dType, tType} end

AbstractStochasticIntegratorRK{DT,TT,D,M,S} = Union{StochasticIntegratorRK{DT,TT,D,M,S}, StochasticIntegratorPRK{DT,TT,D,M,S}}

@inline Integrators.parameters(integrator::AbstractStochasticIntegratorRK) = integrator.params
@inline Integrators.equation(integrator::AbstractStochasticIntegratorRK, i::Symbol) = integrator.params.equs[i]
@inline Integrators.equations(integrator::AbstractStochasticIntegratorRK) = integrator.params.equs
@inline Integrators.timestep(integrator::AbstractStochasticIntegratorRK) = integrator.params.Δt
@inline Integrators.tableau(integrator::AbstractStochasticIntegratorRK)  = integrator.params.tab

@inline Integrators.nstages(::AbstractStochasticIntegratorRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = S
@inline noisedims(::AbstractStochasticIntegratorRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = M
@inline Base.ndims(::AbstractStochasticIntegratorRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = D
@inline Base.eltype(::AbstractStochasticIntegratorRK{DT}) where {DT} = DT

@inline eachstage(integrator::AbstractStochasticIntegratorRK) = 1:nstages(integrator)
@inline eachnoise(integrator::AbstractStochasticIntegratorRK) = 1:noisedims(integrator)
@inline eachdim(integrator::AbstractStochasticIntegratorRK) = 1:ndims(integrator)
