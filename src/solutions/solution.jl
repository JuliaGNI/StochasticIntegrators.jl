
abstract type StochasticSolution{dType,tType,wType,NQ} <: AbstractSolution{dType,tType} end

conv(sol::StochasticSolution) = error("conv() not implemented for ", typeof(sol))

GeometricBase.nsamples(sol::StochasticSolution) = sol.ns


"Create solution for SDE."
function GeometricSolutions.GeometricSolutions(equation::SDE, Δt, ntime::Int; kwargs...)
    SolutionSDE(equation, Δt, ntime; kwargs...)
end

function GeometricSolutions.GeometricSolutions(equation::SDE, Δt, dW, dZ, ntime::Int; kwargs...)
    SolutionSDE(equation, Δt, dW, dZ, ntime; kwargs...)
end

"Create solution for PSDE."
function GeometricSolutions.GeometricSolutions(equation::Union{PSDE,SPSDE}, Δt, ntime::Int; kwargs...)
    SolutionPSDE(equation, Δt, ntime; kwargs...)
end

function GeometricSolutions.GeometricSolutions(equation::Union{PSDE,SPSDE}, Δt, dW, dZ, ntime::Int; kwargs...)
    SolutionPSDE(equation, Δt, dW, dZ, ntime; kwargs...)
end

# "Create ensemble solution for SDE."
# function GeometricSolutions.EnsembleSolution(equation::SDE, Δt, ntime::Int; kwargs...)
#     PSolutionSDE(equation, Δt, ntime; kwargs...)
# end

# "Create parallel solution for PSDE."
# function GeometricSolutions.EnsembleSolution(equation::Union{PSDE,SPSDE}, Δt, ntime::Int; kwargs...)
#     PSolutionPSDE(equation, Δt, ntime; kwargs...)
# end
