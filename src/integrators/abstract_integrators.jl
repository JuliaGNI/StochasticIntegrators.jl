
abstract type StochasticIntegrator{dType, tType} <: Integrator{dType, tType} end

abstract type SDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type PSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type SPSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end

noisedims(integrator::StochasticIntegrator) = error("noisedims() not implemented for ", typeof(integrator))

"Create AtomicSolution for SDE."
function Solutions.AtomicSolution(solution::SolutionSDE{AT,TT}, integrator::Integrator) where {DT, TT, AT <: AbstractArray{DT}}
    AtomicSolutionSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm), get_internal_variables(integrator))
end

"Create AtomicSolution for partitioned SDE."
function Solutions.AtomicSolution(solution::SolutionPSDE{AT,TT}, integrator::Integrator) where {DT, TT, AT <: AbstractArray{DT}}
    AtomicSolutionPSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm), get_internal_variables(integrator))
end

