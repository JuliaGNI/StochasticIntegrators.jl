using RungeKutta: Tableau, TableauLobattoIIIA, TableauLobattoIIIB
using StochasticIntegrators
using Test

@testset "$(rpad("Tableaus", 80))" begin
    @test TableauStochasticEuler() isa TableauSERK
    @test TableauStochasticHeun() isa TableauSERK
    @test TableauPlaten() isa TableauSERK
    @test TableauBurrageR2() isa TableauSERK
    @test TableauBurrageCL() isa TableauSERK
    @test TableauBurrageE1() isa TableauSERK
    @test TableauBurrageG5() isa TableauSERK

    @test TableauStochasticGLRK(1) isa TableauSIRK
    @test TableauStochasticGLRK(2) isa TableauSIRK
    @test TableauStochasticDIRK() isa TableauSIRK

    @test TableauStochasticSymplecticEuler() isa TableauSIPRK
    @test TableauStochasticStoermerVerlet() isa TableauSIPRK

    @test TableauStochasticLobattoIIIABD2() isa TableauSISPRK
    @test TableauModifiedStochasticStoermerVerlet() isa TableauSISPRK

    @test TableauRoesslerRS1() isa TableauWERK
    @test TableauRoesslerRS2() isa TableauWERK

    @test TableauSRKw1() isa TableauWIRK
    @test TableauSRKw2() isa TableauWIRK
end

@testset "$(rpad("Second diffusion tableau", 80))" begin
    # only the higher-order Burrage methods carry the ΔZ terms
    @test !hasdiffusion2(TableauStochasticEuler())
    @test !hasdiffusion2(TableauStochasticHeun())
    @test !hasdiffusion2(TableauBurrageR2())
    @test hasdiffusion2(TableauBurrageCL())
    @test hasdiffusion2(TableauBurrageE1())
    @test hasdiffusion2(TableauBurrageG5())
end

@testset "$(rpad("Methods", 80))" begin
    serk = [StochasticEuler(), StochasticHeun(), Platen(),
        BurrageR2(), BurrageCL(), BurrageE1(), BurrageG5()]
    sirk = [StochasticGLRK(1), StochasticGLRK(2), StochasticDIRK()]
    siprk = [StochasticSymplecticEuler(), StochasticStoermerVerlet()]
    sisprk = [StochasticLobattoIIIABD2(), ModifiedStochasticStoermerVerlet()]
    werk = [RoesslerRS1(), RoesslerRS2()]
    wirk = [SRKw1(), SRKw2()]

    for m in serk
        @test m isa SERK
    end
    for m in sirk
        @test m isa SIRK
    end
    for m in siprk
        @test m isa SIPRK
    end
    for m in sisprk
        @test m isa SISPRK
    end
    for m in werk
        @test m isa WERK
    end
    for m in wirk
        @test m isa WIRK
    end

    # the method hierarchy
    for m in [serk; sirk; werk; wirk]
        @test m isa SDEMethod
        @test m isa StochasticMethod
        @test issdemethod(m)
    end
    for m in siprk
        @test m isa PSDEMethod
        @test ispsdemethod(m)
    end
    for m in sisprk
        @test m isa SPSDEMethod
        @test isspsdemethod(m)
    end

    # convergence is a property of the scheme, never of the caller
    for m in [serk; sirk; siprk; sisprk]
        @test convergence(m) == :strong
    end
    for m in [werk; wirk]
        @test convergence(m) == :weak
    end

    # explicit vs implicit
    for m in [serk; werk]
        @test isexplicit(m)
        @test !isimplicit(m)
    end
    for m in [sirk; siprk; sisprk; wirk]
        @test !isexplicit(m)
        @test isimplicit(m)
    end

    # a stochastic problem has no past to extrapolate to
    for m in [serk; sirk; siprk; sisprk; werk; wirk]
        @test default_extrapolation(m) isa NoExtrapolation
    end

    @test nstages(StochasticGLRK(3)) == 3
    @test truncation(StochasticGLRK(1)) == 0
    @test truncation(StochasticGLRK(1; K = 2)) == 2
end

@testset "$(rpad("Shared weights of a split tableau", 80))" begin
    # the Lagrange-d'Alembert conditions force the q equation and the Hamiltonian part of the p
    # equation to share their weight vectors; a tableau that violates this is rejected
    tab = TableauModifiedStochasticStoermerVerlet()
    @test tab.qdrift.b == tab.pdrift1.b
    @test tab.qdiff.b == tab.pdiff1.b

    bad = Tableau(:bad, 1, [0.0 0.0; 0.5 0.5], [0.3, 0.7], [0.0, 1.0])
    @test_throws AssertionError TableauSISPRK(:invalid,
        TableauLobattoIIIA(2), TableauLobattoIIIA(2),
        bad, TableauLobattoIIIB(2), TableauLobattoIIIB(2), TableauLobattoIIIB(2))
end
