
#*****************************************************************************#
# Initialization functions for stochastic integrators                         #
#*****************************************************************************#

"Create integrator for stochastic explicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauSERK, Δt; kwargs...)
    IntegratorSERK(equation, tableau, Δt; kwargs...)
end

"Create integrator for weak explicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauWERK, Δt; kwargs...)
    IntegratorWERK(equation, tableau, Δt; kwargs...)
end

"Create integrator for stochastic fully implicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauSIRK, Δt; kwargs...)
    IntegratorSIRK(equation, tableau, Δt; kwargs...)
end

"Create integrator for stochastic fully implicit partitioned Runge-Kutta tableau."
function Integrator(equation::PSDE, tableau::TableauSIPRK, Δt; kwargs...)
    IntegratorSIPRK(equation, tableau, Δt; kwargs...)
end

"Create integrator for stochastic fully implicit split partitioned Runge-Kutta tableau."
function Integrator(equation::SPSDE, tableau::TableauSISPRK, Δt; kwargs...)
    IntegratorSISPRK(equation, tableau, Δt; kwargs...)
end

"Create integrator for weak fully implicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauWIRK, Δt; kwargs...)
    IntegratorWIRK(equation, tableau, Δt; kwargs...)
end


#*****************************************************************************#
# Integration functions for stochastic integrators                         #
#*****************************************************************************#

function integrate!(int::StochasticIntegrator{DT,TT}, sol::Solution{AT,TT}, asol::AtomicSolution{DT,TT}, m::Int, n::Int) where {DT, TT, AT <: AbstractArray{DT}}
    # copy the increments of the Brownian Process
    get_increments!(sol, asol, n, m)

    integrate_common!(int, sol, asol, m, n)
end
