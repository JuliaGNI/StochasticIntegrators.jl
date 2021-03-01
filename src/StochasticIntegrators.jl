module StochasticIntegrators

    using HDF5
    using LinearAlgebra: dot, mul!
    using OffsetArrays
    using Reexport
    using SimpleSolvers


    @reexport using GeometricIntegrators.Common
    @reexport using GeometricIntegrators.Config
    @reexport using GeometricIntegrators.Utils


    import GeometricIntegrators.Equations: AbstractEquationSDE, AbstractEquationPSDE,
                                           SDE, PSDE, SPSDE, get_function_tuple


    import GeometricIntegrators.Integrators
    import GeometricIntegrators.Integrators: Integrator, Parameters
    
    import GeometricIntegrators.Integrators: integrate!, integrate_common!, get_internal_variables

    import GeometricIntegrators.Integrators: IntegratorCache, CacheDict, CacheType
    
    import GeometricIntegrators.Integrators: AbstractTableauERK, AbstractTableauIRK, CoefficientsRK

    import GeometricIntegrators.Integrators: create_internal_stage_vector, create_internal_stage_matrix,
                                             create_internal_stage_vector_with_zero, create_nonlinear_solver


    import GeometricIntegrators.Solutions

    import GeometricIntegrators.Solutions: AtomicSolution, Solution, ParallelSolution
    
    import GeometricIntegrators.Solutions: TimeSeries, DataSeries, SDataSeries, PDataSeries
    
    import GeometricIntegrators.Solutions: hdf5, timesteps, nsave, counter, offset, lastentry,
                                           DEFAULT_NSAVE, DEFAULT_NWRITE

    import GeometricIntegrators.Solutions: get_initial_conditions, get_initial_conditions!, set_initial_conditions!,
                                           get_solution, get_solution!, set_solution!,
                                           get_data!, set_data!, update!, compute_timeseries!,
                                           createHDF5, create_hdf5, create_hdf5!, save_attributes,
                                           init_timeteps_in_hdf5, init_solution_in_hdf5, copy_solution_to_hdf5


    import GeometricIntegrators.Tableaus: CoefficientsGLRK, CoefficientsLobattoIIIA, CoefficientsLobattoIIIB, CoefficientsLobattoIIID


    const DEFAULT_SCONV = :strong
    
    export SemiMartingale, WienerProcess, generate_wienerprocess!, conv

    include("solutions/wienerprocess.jl")

    export AtomicSolution, AtomicSolutionSDE, AtomicSolutionPSDE
    export update!, cut_periodic_solution!,
           get_increments, get_increments!, set_increments!

    include("solutions/atomic_solution_sde.jl")
    include("solutions/atomic_solution_psde.jl")

    export Solution, StochasticSolution
    export SolutionSDE, SSolutionSDE, PSolutionSDE, SolutionPSDE, SSolutionPSDE, PSolutionPSDE
    export hdf5, timesteps, nsave, counter, offset, lastentry
    export get_initial_conditions, get_initial_conditions!, set_initial_conditions!,
           get_solution, get_solution!, set_solution!,
           create_hdf5, create_hdf5!


    include("solutions/solution.jl")
    include("solutions/solution_sde.jl")
    include("solutions/solution_psde.jl")
    include("solutions/solution_hdf5.jl")

    include("solutions/atomic_solution_constructors.jl")
    

    export Parameters

    export Integrator, StochasticIntegrator
    export SDEIntegrator, PSDEIntegrator,
           StochasticIntegratorRK, StochasticIntegratorPRK
    
    export nstages, noisedims, integrate!

    export IntegratorSERK, TableauSERK
    export IntegratorSIRK, TableauSIRK
    export IntegratorSIPRK, TableauSIPRK
    export IntegratorSISPRK, TableauSISPRK
    export IntegratorWERK, TableauWERK
    export IntegratorWIRK, TableauWIRK

    include("integrators/abstract_integrators.jl")
    include("integrators/abstract_runge_kutta.jl")
    include("integrators/integrator_cache.jl")

    include("integrators/integrators_serk.jl")
    include("integrators/integrators_sirk.jl")
    include("integrators/integrators_siprk.jl")
    include("integrators/integrators_sisprk.jl")
    include("integrators/integrators_werk.jl")
    include("integrators/integrators_wirk.jl")
    include("integrators/common.jl")
    include("integrators/integrators.jl")


    include("tableaus/tableaus_sirk.jl")

    export TableauStochasticGLRK, TableauStochasticDIRK

    include("tableaus/tableaus_siprk.jl")

    export TableauStochasticStoermerVerlet, TableauStochasticSymplecticEuler
    
    include("tableaus/tableaus_sisprk.jl")

    export TableauStochasticLobattoIIIABD2, TableauModifiedStochasticStoermerVerlet

    include("tableaus/tableaus_serk.jl")

    export TableauPlaten, TableauBurrageR2, TableauBurrageCL
    export TableauBurrageE1, TableauBurrageG5, TableauStochasticHeun
    export TableauStochasticEuler

    include("tableaus/tableaus_werk.jl")

    export TableauRoesslerRS1, TableauRoesslerRS2

    include("tableaus/tableaus_wirk.jl")

    export TableauSRKw1, TableauSRKw2
    
end
