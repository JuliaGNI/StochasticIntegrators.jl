using SafeTestsets

@safetestset "Atomic Solutions                                                                " begin include("atomic_solutions_tests.jl") end
@safetestset "Stochastic Solutions                                                            " begin include("stochastic_solutions_tests.jl") end
@safetestset "Stochastic Tableaus                                                             " begin include("stochastic_tableaus_tests.jl") end
@safetestset "Stochastic Integrators                                                          " begin include("stochastic_integrators_tests.jl") end
