using GeometricProblems.KuboOscillator
using StochasticIntegrators
using Random
using Statistics
using Test

@testset "$(rpad("Noise processes", 80))" begin
    @test noisedims(WienerProcess(1)) == 1
    @test noisedims(WienerProcess(4)) == 4
    @test_throws AssertionError WienerProcess(0)

    gp = GridProcess(reshape([1.0, 2.0, 3.0], 1, 3))
    @test noisedims(gp) == 1
    @test all(iszero, gp.ΔZ)
    @test GridProcess(ones(2, 5), zeros(2, 5)) == GridProcess(ones(2, 5), zeros(2, 5))
    @test_throws AssertionError GridProcess(ones(2, 5), zeros(2, 4))
end

@testset "$(rpad("Strong increments", 80))" begin
    Δt = 0.01
    n = 200_000
    rng = Xoshiro(1234)
    process = WienerProcess(1)

    ΔW = zeros(1)
    ΔZ = zeros(1)
    W = zeros(n)
    Z = zeros(n)

    for i in 1:n
        sample_noise!(ΔW, ΔZ, process, Val(:strong), Δt, i, rng)
        W[i] = ΔW[1]
        Z[i] = ΔZ[1]
    end

    # ΔW ~ N(0, Δt)
    @test isapprox(mean(W), 0; atol = 5e-4)
    @test isapprox(var(W), Δt; rtol = 2e-2)

    # ΔZ = ∫∫dW dt has mean 0 and variance Δt³/3
    @test isapprox(mean(Z), 0; atol = 5e-6)
    @test isapprox(var(Z), Δt^3 / 3; rtol = 3e-2)

    # and the two are correlated, E[ΔW ΔZ] = Δt²/2
    @test isapprox(mean(W .* Z), Δt^2 / 2; rtol = 3e-2)
end

@testset "$(rpad("Weak increments", 80))" begin
    Δt = 0.01
    n = 200_000
    rng = Xoshiro(1234)
    process = WienerProcess(1)

    ΔW = zeros(1)
    ΔZ = zeros(1)
    W = zeros(n)
    Z = zeros(n)

    for i in 1:n
        sample_noise!(ΔW, ΔZ, process, Val(:weak), Δt, i, rng)
        W[i] = ΔW[1]
        Z[i] = ΔZ[1]
    end

    # Î is supported on exactly three points, Ĩ on exactly two
    @test sort(unique(W)) ≈ [-sqrt(3Δt), 0.0, sqrt(3Δt)]
    @test sort(unique(Z)) ≈ [-sqrt(Δt), sqrt(Δt)]

    # P(Î = 0) = 2/3, P(Î = ±√(3Δt)) = 1/6
    @test isapprox(count(iszero, W) / n, 2 / 3; atol = 5e-3)
    @test isapprox(count(>(0), W) / n, 1 / 6; atol = 5e-3)

    # matching the moments of the true increments is the whole point of a weak scheme
    @test isapprox(mean(W), 0; atol = 5e-4)
    @test isapprox(var(W), Δt; rtol = 2e-2)
    @test isapprox(mean(Z), 0; atol = 5e-4)
    @test isapprox(var(Z), Δt; rtol = 2e-2)
end

@testset "$(rpad("Prescribed increments", 80))" begin
    ΔWs = [0.1 0.2 0.3; -0.1 -0.2 -0.3]
    ΔZs = [1.0 2.0 3.0; -1.0 -2.0 -3.0]
    process = GridProcess(ΔWs, ΔZs)

    ΔW = zeros(2)
    ΔZ = zeros(2)

    for n in 1:3
        sample_noise!(ΔW, ΔZ, process, Val(:strong), 0.1, n, Xoshiro(0))
        @test ΔW == ΔWs[:, n]
        @test ΔZ == ΔZs[:, n]
    end

    # a GridProcess ignores the convergence mode — its increments are given, not drawn
    sample_noise!(ΔW, ΔZ, process, Val(:weak), 0.1, 2, Xoshiro(0))
    @test ΔW == ΔWs[:, 2]

    # and it refuses a step it has no increments for, rather than reading out of bounds
    @test_throws AssertionError sample_noise!(
        ΔW, ΔZ, process, Val(:strong), 0.1, 4, Xoshiro(0))
end

@testset "$(rpad("Increment truncation", 80))" begin
    ΔW = [-3.0, -0.5, 0.0, 0.5, 3.0]

    @test truncate_increments!(copy(ΔW), 0) == ΔW          # A = 0 disables truncation
    @test truncate_increments!(copy(ΔW), 1.0) == [-1.0, -0.5, 0.0, 0.5, 1.0]

    # the Milstein-Tretyakov bound A = √(2KΔt|log Δt|)
    @test StochasticIntegrators.truncation_bound(0, 0.01) == 0
    @test StochasticIntegrators.truncation_bound(1, 0.01) ≈ sqrt(2 * 0.01 * abs(log(0.01)))
end

@testset "$(rpad("Reproducibility", 80))" begin
    # the same seed gives the same trajectory; different seeds do not
    prob = sdeproblem()
    a = integrate(prob, BurrageE1(; rng = Xoshiro(42)))
    b = integrate(prob, BurrageE1(; rng = Xoshiro(42)))
    c = integrate(prob, BurrageE1(; rng = Xoshiro(43)))

    @test a.q[end] == b.q[end]
    @test a.q[end] != c.q[end]
end
