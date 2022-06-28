
using StochasticIntegrators
using GeometricProblems.KuboOscillator
using Test


Δt = .1
t0 = 0.
x0 = rand(2)
q0 = rand(1)
p0 = q0.^2
λ0 = rand(1)
λ1 = rand(1)
v0 = rand(2)
y0 = rand(1)
z0 = rand(1)
ΔW = rand(3)
ΔZ = rand(3)


@testset "$(rpad("Atomic Solution Constructors",80))" begin
    sde   = kubo_oscillator_sde_1()
    psde  = kubo_oscillator_psde_1()
    spsde = kubo_oscillator_spsde_1()

    @test typeof(AtomicSolution(sde))   <: AtomicSolutionSDE
    @test typeof(AtomicSolution(psde))  <: AtomicSolutionPSDE
    @test typeof(AtomicSolution(spsde)) <: AtomicSolutionPSDE
end



@testset "$(rpad("Atomic SDE Solution",80))" begin
    asol = AtomicSolutionSDE(t0, x0, ΔW, ΔZ)
    @test get_solution(asol) == (zero(t0), zero(x0))

    set_solution!(asol, (t0, x0))
    @test get_solution(asol) == (t0, x0)
    @test asol.t  == t0
    @test asol.q  == x0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(x0)

    set_increments!(asol, (ΔW, ΔZ))
    δW = zero(ΔW)
    δZ = zero(ΔZ)
    get_increments!(asol, δW, δZ)
    @test get_increments(asol) == (ΔW, ΔZ)
    @test asol.ΔW == ΔW
    @test asol.ΔZ == ΔZ
    @test δW == ΔW
    @test δZ == ΔZ

    reset!(asol, Δt)
    @test asol.t̄ == t0
    @test asol.q̄ == x0

    update!(asol, v0)
    @test asol.t == t0  + Δt
    @test asol.q == x0 .+ v0
end



@testset "$(rpad("Atomic PSDE Solution",80))" begin
    asol = AtomicSolutionPSDE(t0, q0, p0, ΔW, ΔZ)
    @test get_solution(asol) == (zero(t0), zero(q0), zero(p0))
    @test get_increments(asol) == (zero(ΔW), zero(ΔZ))

    set_solution!(asol, (t0, q0, p0))
    @test get_solution(asol) == (t0, q0, p0)
    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.t̄  == zero(t0)
    @test asol.q̄  == zero(q0)
    @test asol.p̄  == zero(p0)

    set_increments!(asol, (ΔW, ΔZ))
    δW = zero(ΔW)
    δZ = zero(ΔZ)
    get_increments!(asol, δW, δZ)
    @test get_increments(asol) == (ΔW, ΔZ)
    @test asol.ΔW == ΔW
    @test asol.ΔZ == ΔZ
    @test δW == ΔW
    @test δZ == ΔZ

    reset!(asol, Δt)
    @test asol.t̄ == t0
    @test asol.q̄ == q0
    @test asol.p̄ == p0

    update!(asol, y0, z0)
    @test asol.t == t0  + Δt
    @test asol.q == q0 .+ y0
    @test asol.p == p0 .+ z0
end
