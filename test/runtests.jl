using SafeTestsets

@safetestset "Methods and Tableaus                                                            " begin
    include("methods_tests.jl")
end
@safetestset "Noise Processes                                                                 " begin
    include("processes_tests.jl")
end
@safetestset "Stochastic Integrators                                                          " begin
    include("integrators_tests.jl")
end
@safetestset "Multidimensional Noise                                                          " begin
    include("multidimensional_noise_tests.jl")
end
@safetestset "Package Quality                                                                 " begin
    include("aqua_tests.jl")
end
