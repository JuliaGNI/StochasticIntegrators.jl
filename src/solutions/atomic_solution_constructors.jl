
"Create SolutionStep for SDE."
function Solutions.SolutionStep(equation::AbstractEquationSDE{DT,TT}) where {DT,TT}
    SolutionStepSDE(equation.t₀, equation.q₀[begin], zeros(DT,equation.m), zeros(DT,equation.m))
end

"Create SolutionStep for SDE."
function Solutions.SolutionStep(solution::SolutionSDE{AT,TT}) where {DT, TT, AT <: AbstractArray{DT}}
    SolutionStepSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm))
end

"Create SolutionStep for partitioned SDE."
function Solutions.SolutionStep(equation::AbstractEquationPSDE{DT,TT}) where {DT,TT}
    SolutionStepPSDE(equation.t₀, equation.q₀[begin], equation.p₀[begin], zeros(DT,equation.m), zeros(DT,equation.m))
end

"Create SolutionStep for partitioned SDE."
function Solutions.SolutionStep(solution::SolutionPSDE{AT,TT}) where {DT, TT, AT <: AbstractArray{DT}}
    SolutionStepPSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm))
end

