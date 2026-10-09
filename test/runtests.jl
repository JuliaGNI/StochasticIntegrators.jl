using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Methods and tableaus" include("methods.jl")
    @safetestset "Noise processes" include("processes.jl")
    @safetestset "Stochastic integrators" include("integrators/integrators.jl")
    @safetestset "Multidimensional noise" include("integrators/multidimensional_noise.jl")
end
if "doctests" in GROUPS
    @safetestset "Doctests" include("quality/doctests.jl")
end
