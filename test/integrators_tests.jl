using StochasticIntegrators
using GeometricIntegrators: Gauss
using GeometricProblems.KuboOscillator
using Random
using Test

import GeometricProblems.KuboOscillator as Kubo

include("utils.jl")

const Δt = Kubo.Δt
const nt = Kubo.nt

@testset "$(rpad("Integrator construction", 80))" begin
    @test GeometricIntegrator(sdeproblem(), BurrageE1()) isa GeometricIntegrator
    @test GeometricIntegrator(sdeproblem(), StochasticGLRK(1)) isa GeometricIntegrator
    @test GeometricIntegrator(sdeproblem(), RoesslerRS1()) isa GeometricIntegrator
    @test GeometricIntegrator(sdeproblem(), SRKw1()) isa GeometricIntegrator
    @test GeometricIntegrator(psdeproblem(), StochasticStoermerVerlet()) isa
          GeometricIntegrator
    @test GeometricIntegrator(spsdeproblem(), ModifiedStochasticStoermerVerlet()) isa
          GeometricIntegrator
end

# The Kubo oscillator conserves energy exactly along every path, so these tolerances measure the
# integrator alone. They are the tolerances the pre-rewrite suite asserted, kept as regression
# anchors for the port.

@testset "$(rpad("SERK integrators", 80))" begin
    @test rel_energy_err(integrate(sdeproblem(), BurrageE1())) < 2E-6
    @test max_rel_energy_err(integrate(sdeensemble(), BurrageE1())) < 2E-6
end

@testset "$(rpad("SIRK integrators", 80))" begin
    @test rel_energy_err(integrate(sdeproblem(), StochasticGLRK(1))) < 1E-14
    @test max_rel_energy_err(integrate(sdeensemble(), StochasticGLRK(1))) < 1E-14
    @test rel_energy_err(integrate(sdeproblem(), StochasticDIRK())) < 1E-14
end

@testset "$(rpad("WERK integrators", 80))" begin
    @test rel_energy_err(integrate(sdeproblem(), RoesslerRS1())) < 1E-5
    @test max_rel_energy_err(integrate(sdeensemble(), RoesslerRS1())) < 1E-5
end

@testset "$(rpad("WIRK integrators", 80))" begin
    @test rel_energy_err(integrate(sdeproblem(), SRKw1())) < 1E-14
    @test max_rel_energy_err(integrate(sdeensemble(), SRKw1())) < 1E-14
end

@testset "$(rpad("SIPRK integrators", 80))" begin
    @test rel_energy_err_pq(integrate(psdeproblem(), StochasticStoermerVerlet())) < 1E-5
    @test max_rel_energy_err_pq(integrate(psdeensemble(), StochasticStoermerVerlet())) <
          2E-5
end

@testset "$(rpad("SISPRK integrators", 80))" begin
    @test rel_energy_err_pq(integrate(spsdeproblem(), ModifiedStochasticStoermerVerlet())) <
          2E-2
    @test max_rel_energy_err_pq(integrate(
        spsdeensemble(), ModifiedStochasticStoermerVerlet())) < 2E-2
    @test rel_energy_err_pq(integrate(spsdeproblem(), StochasticLobattoIIIABD2())) < 2E-2
end

@testset "$(rpad("Zero noise reproduces the deterministic method", 80))" begin
    # Driven by prescribed zero increments an SDE is its own drift, so a stochastic method must
    # reproduce the deterministic method it is built from, up to round-off. This is the sharpest
    # check on the drift half of a scheme. Not bit-identical, and not asserted as such: the
    # stochastic update adds a diffusion term that happens to be zero, so it sums its increment in
    # a different order and the two may differ in the last bit.
    local tspan = (0.0, Δt * nt)
    local nw = grid_length(tspan, Δt)

    zeronoise = GridProcess(zeros(1, nw), zeros(1, nw))
    sde0 = SDEProblem(Kubo.kubo_oscillator_sde_v, Kubo.kubo_oscillator_sde_B, zeronoise,
        tspan, Δt, Kubo.q_init_A; parameters = Kubo.default_parameters())

    sol_sto = integrate(sde0, StochasticGLRK(1))
    sol_det = integrate(odeproblem(), Gauss(1))

    @test sol_sto.q[end] ≈ sol_det.q[end] atol=1E-14
end

@testset "$(rpad("Splitting the forcing does not change the answer", 80))" begin
    # The damped Kubo oscillator, split (SPSDE) and unsplit (PSDE), on one common sample path.
    # The two describe the same dynamics, so the trajectories must agree.
    #
    # This is the regression test for a defect in the pre-rewrite SISPRK residual, where a line
    # continuation beginning with `+` silently dropped the f2 and G2 contributions from the
    # nonlinear system. The undamped Kubo problems have f2 ≡ 0 and G2 ≡ 0, so nothing in the old
    # suite could see it. The damping is deliberately large here — the paper's γ = 0.001 would
    # make the dropped term small enough to hide inside the solver tolerance.
    local nsteps = 200
    local h = 0.05
    local par = (ν = 0.5, γ = 0.5)
    local tspan = (0.0, h * nsteps)

    local nw = grid_length(tspan, h)
    local path = randn(Xoshiro(42), 1, nw) .* sqrt(h)
    local gp = GridProcess(path, zeros(1, nw))

    local q₀ = [2.0]
    local p₀ = [0.0]

    psde = PSDEProblem(Kubo.kubo_oscillator_psde_v, Kubo.kubo_oscillator_damped_f,
        Kubo.kubo_oscillator_psde_B, Kubo.kubo_oscillator_damped_G,
        gp, tspan, h, q₀, p₀; parameters = par)

    spsde = SPSDEProblem(Kubo.kubo_oscillator_spsde_v, Kubo.kubo_oscillator_spsde_f1,
        Kubo.kubo_oscillator_damped_f2, Kubo.kubo_oscillator_spsde_B,
        Kubo.kubo_oscillator_spsde_G1, Kubo.kubo_oscillator_damped_G2,
        gp, tspan, h, q₀, p₀; parameters = par)

    unsplit = integrate(psde, StochasticStoermerVerlet())
    split = integrate(spsde, ModifiedStochasticStoermerVerlet(0.5))

    @test split.q[end] ≈ unsplit.q[end] atol=1E-12
    @test split.p[end] ≈ unsplit.p[end] atol=1E-12

    # and both track the exact solution of the damped oscillator along that same path.
    # `grid_length` may allow one increment more than the run consumes, so the Wiener path is
    # summed over the steps actually taken rather than over the whole prescribed grid.
    local W = cumsum(vec(path))[ntime(split)]
    local qex, pex = exact_solution(tspan[end], W, q₀[begin], p₀[begin], tspan[begin], par)
    @test split.q[end][begin]≈qex atol=1E-2
    @test split.p[end][begin]≈pex atol=1E-2
end

@testset "$(rpad("Prescribed noise is reproducible", 80))" begin
    # A GridProcess fixes the path, so two runs of the same problem must agree exactly even
    # though each method carries its own random number generator.
    local tspan = (0.0, Δt * nt)
    local nw = grid_length(tspan, Δt)
    local path = randn(Xoshiro(7), 1, nw) .* sqrt(Δt)
    local gp = GridProcess(path, zeros(1, nw))
    local prob = SDEProblem(Kubo.kubo_oscillator_sde_v, Kubo.kubo_oscillator_sde_B, gp,
        tspan, Δt, Kubo.q_init_A; parameters = Kubo.default_parameters())

    @test integrate(prob, BurrageE1()).q[end] == integrate(prob, BurrageE1()).q[end]
end
