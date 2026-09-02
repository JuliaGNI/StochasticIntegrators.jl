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

# The increments are drawn in the element type of the arrays that receive them, not in the type
# of the time step. Those need not agree — a problem's data type and its time type are separate
# — and sampling in `typeof(Δt)` would then convert on every store. `√Δt` is likewise formed
# once, in that same type, and `Δt^(3/2)` written as `Δt·√Δt` to keep the exponentiation off a
# `Rational` power.

function sample_noise!(ΔW, ΔZ, ::WienerProcess, ::Val{:strong}, Δt, n, rng)
    local DT = eltype(ΔW)
    local h = convert(DT, Δt)
    local sqrth = sqrt(h)
    local sqrt3 = sqrt(convert(DT, 3))

    for r in eachindex(ΔW, ΔZ)
        χ = randn(rng, DT)
        η = randn(rng, DT)
        ΔW[r] = χ * sqrth
        ΔZ[r] = h * sqrth / 2 * (χ + η / sqrt3)
    end
    (ΔW, ΔZ)
end

function sample_noise!(ΔW, ΔZ, ::WienerProcess, ::Val{:weak}, Δt, n, rng)
    local DT = eltype(ΔW)
    local h = convert(DT, Δt)
    local sqrth = sqrt(h)
    local sqrt3h = sqrt(3h)

    # The thresholds are converted once, for the same reason the square roots are hoisted:
    # comparing a `Float32` draw against a `Rational` promotes through a widening path that
    # allocates, where comparing against a `Float32` does not. The comparison is unaffected —
    # `rand(rng, Float32)` returns a multiple of 2^-24 and `Float32(1//6)` has an ulp of 2^-26,
    # so no draw can land on a threshold and the two forms cannot disagree.
    local lo = convert(DT, 1//6)
    local hi = convert(DT, 5//6)
    local half = convert(DT, 1//2)

    for r in eachindex(ΔW, ΔZ)
        χ = rand(rng, DT)
        η = rand(rng, DT)

        # Î: P(±√(3Δt)) = 1/6 each, P(0) = 2/3
        ΔW[r] = χ < lo ? -sqrt3h : (χ > hi ? +sqrt3h : zero(DT))

        # Ĩ: P(±√Δt) = 1/2 each
        ΔZ[r] = η < half ? -sqrth : +sqrth
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
