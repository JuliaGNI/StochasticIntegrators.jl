@doc raw"""
    sample_noise!(ΔW, ΔZ, process, convergence, Δt, n, rng)

Draw the increments driving one time step, writing them into `ΔW` and `ΔZ`.

`n` is the index of the step being taken, counted from 1; it is used only by a
`GridProcess`, which prescribes its increments in advance.

## Strong increments

For a `WienerProcess` under `:strong` convergence the increments are
the true ones, drawn from
```math
\Delta W^r = \chi^r \sqrt{\Delta t}, \qquad
\Delta Z^r = \tfrac{1}{2} \Delta t^{3/2} \left( \chi^r + \tfrac{1}{\sqrt 3} \eta^r \right),
```
with ``\chi^r, \eta^r`` independent standard normals. ``\Delta W^r`` is the Wiener increment
``\int_{t_n}^{t_{n+1}} dW^r`` and ``\Delta Z^r`` the iterated integral
``\int_{t_n}^{t_{n+1}} \int_{t_n}^{t} dW^r(\xi) \, dt``; the coefficients are the standard joint
law of the pair. Only schemes carrying a second diffusion tableau use ``\Delta Z``.

## Weak increments

Under `:weak` convergence the true increments are replaced by discrete random variables matching
their low-order moments, which is all a weak scheme is derived against and is much cheaper to
sample:
```math
P\left( \hat I^r = \pm \sqrt{3 \Delta t} \right) = \tfrac{1}{6}, \qquad
P\left( \hat I^r = 0 \right) = \tfrac{2}{3}, \qquad
P\left( \tilde I^r = \pm \sqrt{\Delta t} \right) = \tfrac{1}{2} .
```
``\hat I`` is returned in `ΔW` and ``\tilde I`` in `ΔZ`. Note that these are *not* the same
quantities as in the strong case — they occupy the same arrays because they play the same role in
the update formulas, but ``\tilde I`` is a two-point variable, not an iterated integral.
"""
function sample_noise! end

function sample_noise!(ΔW, ΔZ, ::WienerProcess, ::Val{:strong}, Δt, n, rng)
    for r in eachindex(ΔW, ΔZ)
        χ = randn(rng, typeof(Δt))
        η = randn(rng, typeof(Δt))
        ΔW[r] = χ * sqrt(Δt)
        ΔZ[r] = Δt^(3 // 2) / 2 * (χ + η / sqrt(3))
    end
    (ΔW, ΔZ)
end

function sample_noise!(ΔW, ΔZ, ::WienerProcess, ::Val{:weak}, Δt, n, rng)
    for r in eachindex(ΔW, ΔZ)
        χ = rand(rng, typeof(Δt))
        η = rand(rng, typeof(Δt))

        # Î: P(±√(3Δt)) = 1/6 each, P(0) = 2/3
        ΔW[r] = χ < 1 // 6 ? -sqrt(3Δt) : (χ > 5 // 6 ? +sqrt(3Δt) : zero(Δt))

        # Ĩ: P(±√Δt) = 1/2 each
        ΔZ[r] = η < 1 // 2 ? -sqrt(Δt) : +sqrt(Δt)
    end
    (ΔW, ΔZ)
end

function sample_noise!(ΔW, ΔZ, process::GridProcess, ::Val, Δt, n, rng)
    @assert 1≤n≤size(process.ΔW, 2) "step $n is outside the $(size(process.ΔW, 2)) steps prescribed by the GridProcess"

    for r in eachindex(ΔW, ΔZ)
        ΔW[r] = process.ΔW[r, n]
        ΔZ[r] = process.ΔZ[r, n]
    end
    (ΔW, ΔZ)
end

@doc raw"""
    truncate_increments!(ΔW, A)

Clip every component of `ΔW` to ``[-A, A]`` in place, leaving it untouched when `A == 0`.

See [`truncation`](@ref) for why an implicit scheme wants this.
"""
function truncate_increments!(ΔW, A)
    A > 0 || return ΔW
    for i in eachindex(ΔW)
        ΔW[i] = clamp(ΔW[i], -A, A)
    end
    ΔW
end
