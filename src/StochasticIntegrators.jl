module StochasticIntegrators

using LinearAlgebra: dot, mul!
using Random: Random
using Reexport

@reexport using GeometricBase
@reexport using GeometricEquations
@reexport using GeometricIntegratorsBase
@reexport using GeometricSolutions

# `GeometricBase` and `GeometricSolutions` each export a `solution`, and the two are different
# generic functions rather than one extended by the other. Both are re-exported above, so the
# name would resolve to neither and `solution` would be unusable after `using
# StochasticIntegrators`. It is bound to the `GeometricSolutions` one, which acts on the
# `GeometricSolution` that `integrate` returns; the other, on a `State` or a `SolutionStep`, is
# reachable as `GeometricBase.solution`.
using GeometricSolutions: solution

import GeometricBase: equations, name, noise, noisedims, tableau, timestep

import GeometricEquations: AbstractProblemSDE, AbstractProblemPSDE, AbstractProblemSPSDE,
                           AbstractStochasticProcess, WienerProcess, GridProcess,
                           initial_conditions

import GeometricSolutions: initialtime

import GeometricIntegratorsBase: GeometricIntegrator, GeometricMethod,
                                 Cache, IntegratorCache,
                                 cache, nlsolution, solver, solverstate,
                                 problem, method,
                                 integrate_step!, components!, residual!, update!,
                                 default_solver, check_solver_status,
                                 default_extrapolation, NoExtrapolation,
                                 isexplicit, isimplicit

import RungeKutta: AbstractTableau, Tableau, TableauGauss,
                   TableauLobattoIIIA, TableauLobattoIIIB, TableauLobattoIIID,
                   istrilstrict

import SimpleSolvers: Newton, solve_with_status!

export StochasticMethod, SDEMethod, PSDEMethod, SPSDEMethod
export convergence, truncation, nstages, noisedims
export issdemethod, ispsdemethod, isspsdemethod

include("methods.jl")

export sample_noise!, truncate_increments!

include("processes.jl")

export StochasticIntegratorCache, SDEIntegratorCache, PSDEIntegratorCache

include("integrators/cache.jl")
include("integrators/updates.jl")

export SERK, TableauSERK
export SIRK, TableauSIRK
export WERK, TableauWERK
export WIRK, TableauWIRK
export SIPRK, TableauSIPRK
export SISPRK, TableauSISPRK
export hasdiffusion2

include("integrators/serk.jl")
include("integrators/sirk.jl")
include("integrators/werk.jl")
include("integrators/wirk.jl")
include("integrators/siprk.jl")
include("integrators/sisprk.jl")

export TableauStochasticEuler, TableauStochasticHeun, TableauPlaten,
       TableauBurrageR2, TableauBurrageCL, TableauBurrageE1, TableauBurrageG5
export StochasticEuler, StochasticHeun, Platen,
       BurrageR2, BurrageCL, BurrageE1, BurrageG5

include("tableaus/tableaus_serk.jl")

export TableauStochasticGLRK, TableauStochasticDIRK
export StochasticGLRK, StochasticDIRK

include("tableaus/tableaus_sirk.jl")

export TableauStochasticSymplecticEuler, TableauStochasticStoermerVerlet
export StochasticSymplecticEuler, StochasticStoermerVerlet

include("tableaus/tableaus_siprk.jl")

export TableauStochasticLobattoIIIABD2, TableauModifiedStochasticStoermerVerlet
export StochasticLobattoIIIABD2, ModifiedStochasticStoermerVerlet

include("tableaus/tableaus_sisprk.jl")

export TableauRoesslerRS1, TableauRoesslerRS2
export RoesslerRS1, RoesslerRS2

include("tableaus/tableaus_werk.jl")

export TableauSRKw1, TableauSRKw2
export SRKw1, SRKw2

include("tableaus/tableaus_wirk.jl")

end
