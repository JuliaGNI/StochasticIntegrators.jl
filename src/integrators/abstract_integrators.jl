
abstract type StochasticIntegrator{dType, tType} <: GeometricIntegrator{dType, tType} end

abstract type SDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type PSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type SPSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end

noisedims(integrator::StochasticIntegrator) = error("noisedims() not implemented for ", typeof(integrator))

"Create SolutionStep for SDE."
function Solutions.SolutionStep(solution::SolutionSDE{AT,TT}, integrator::Integrator) where {DT, TT, AT <: AbstractArray{DT}}
    SolutionStepSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm), get_internal_variables(integrator))
end

"Create SolutionStep for partitioned SDE."
function Solutions.SolutionStep(solution::SolutionPSDE{AT,TT}, integrator::Integrator) where {DT, TT, AT <: AbstractArray{DT}}
    SolutionStepPSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm), get_internal_variables(integrator))
end

