using GeometricBase: GeometricBase
using InteractiveUtils: subtypes
using StochasticIntegrators
using Test

# The concrete method types that this package owns, found from `GeometricBase.AbstractMethod`
# down. A parametric type comes back as a `UnionAll`, so a type counts as concrete when its body is
# not abstract.
function owned_method_types(T = GeometricBase.AbstractMethod)
    owned = Type[]
    for S in subtypes(T)
        if !isabstracttype(Base.unwrap_unionall(S)) &&
           Base.moduleroot(parentmodule(S)) === StochasticIntegrators
            push!(owned, S)
        end
        append!(owned, owned_method_types(S))
    end
    return owned
end

const OWNED = unique(owned_method_types())

# One instance of each owned method type, built as `test/methods.jl` builds them, together with an
# instance of every owned type that has a constructor of no argument.
const METHODS = [StochasticEuler(),
    StochasticGLRK(1),
    StochasticSymplecticEuler(),
    StochasticLobattoIIIABD2(),
    RoesslerRS1(),
    SRKw1(),
    [T() for T in OWNED if hasmethod(T, Tuple{})]...]

@testset "$(rpad("Interface conformance", 80))" begin
    @testset "$(nameof(typeof(m)))" for m in METHODS
        @test isempty(GeometricBase.interface_violations(m))
    end

    # a new method type without an instance above turns this red
    @testset "$(nameof(T)) is covered" for T in OWNED
        @test any(m -> m isa T, METHODS)
    end
end
