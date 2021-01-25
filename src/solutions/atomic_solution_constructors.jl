
"Create AtomicSolution for SDE."
function Solutions.AtomicSolution(equation::AbstractEquationSDE{DT,TT}) where {DT,TT}
    AtomicSolutionSDE(equation.t₀, equation.q₀[begin], zeros(DT,equation.m), zeros(DT,equation.m))
end

"Create AtomicSolution for SDE."
function Solutions.AtomicSolution(solution::SolutionSDE{AT,TT}) where {DT, TT, AT <: AbstractArray{DT}}
    AtomicSolutionSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm))
end

"Create AtomicSolution for partitioned SDE."
function Solutions.AtomicSolution(equation::AbstractEquationPSDE{DT,TT}) where {DT,TT}
    AtomicSolutionPSDE(equation.t₀, equation.q₀[begin], equation.p₀[begin], zeros(DT,equation.m), zeros(DT,equation.m))
end

"Create AtomicSolution for partitioned SDE."
function Solutions.AtomicSolution(solution::SolutionPSDE{AT,TT}) where {DT, TT, AT <: AbstractArray{DT}}
    AtomicSolutionPSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm))
end

