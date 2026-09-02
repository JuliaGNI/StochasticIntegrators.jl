```@meta
CurrentModule = StochasticIntegrators
```

# Noise Processes

## Where the noise lives

The driving process is part of the **problem**, not of the integrator. Every `SDE`, `PSDE` and
`SPSDE` carries a `noise` field holding an `AbstractStochasticProcess`, and
`noisedims(problem)` reports how many independent Wiener processes drive it — which is what fixes
the number of columns of ``B`` and ``G``, and what an integrator uses to size its increment
vectors.

The process object says *which* noise, not *which realisation*. Increments are drawn by the
integrator, because only the integrator knows whether the scheme needs strong or weak ones.

`GeometricEquations` provides two:

* `WienerProcess(m)` — ``m`` independent standard Wiener processes, sampled afresh each run;
* `GridProcess(ΔW, ΔZ)` — increments prescribed in advance on the time grid.

## Sampling

```@docs
sample_noise!
```

### Strong increments

For a scheme of strong convergence the increments must be the true ones,

```math
\Delta W^r = \chi^r \sqrt{\Delta t} , \qquad
\Delta Z^r = \tfrac{1}{2} \Delta t^{3/2} \left( \chi^r + \tfrac{1}{\sqrt 3} \eta^r \right) ,
```

with ``\chi^r, \eta^r`` independent standard normals. ``\Delta W^r`` is the Wiener increment and
``\Delta Z^r`` the iterated integral ``\int \int dW^r \, dt``; the coefficients above are their
joint law, giving

```math
\mathbb{E}[\Delta W^2] = \Delta t , \qquad
\mathbb{E}[\Delta Z^2] = \tfrac{1}{3}\Delta t^3 , \qquad
\mathbb{E}[\Delta W \, \Delta Z] = \tfrac{1}{2}\Delta t^2 .
```

Only the schemes carrying a second diffusion tableau use ``\Delta Z``.

### Weak increments

A weak scheme is derived against the *moments* of the increments rather than the increments
themselves, so it can be fed a discrete random variable matching those moments at a fraction of
the cost:

```math
P\left( \hat I^r = \pm \sqrt{3 \Delta t} \right) = \tfrac{1}{6} , \qquad
P\left( \hat I^r = 0 \right) = \tfrac{2}{3} , \qquad
P\left( \tilde I^r = \pm \sqrt{\Delta t} \right) = \tfrac{1}{2} .
```

``\hat I`` takes the place of ``\Delta W`` and ``\tilde I`` that of ``\Delta Z``. They are *not*
the same quantities — ``\tilde I`` is a two-point variable, not an iterated integral — but they
occupy the same arrays because they play the same role in the update formulas.

### Truncation

```@docs
truncation
truncation_bound
truncate_increments!
```

A Gaussian increment is unbounded, and the nonlinear system of an implicit scheme need not be
solvable for arbitrarily large ones. Milstein and Tretyakov's remedy is to clip ``\Delta W`` at
``A = \sqrt{2 K \Delta t \, \lvert \log \Delta t \rvert}``, which leaves the mean-square order
intact because the discarded tail is of higher order than the scheme. It is off by default and
enabled per method with `K`:

```@example noise
using StochasticIntegrators
using GeometricProblems.KuboOscillator

integrate(sdeproblem(), StochasticGLRK(1; K = 1)).q[end]
```

## Prescribed paths

A `GridProcess` fixes the realisation instead of drawing it. It has four uses, and all four are
things that cannot be done with a `WienerProcess`.

!!! warning "How many increments a GridProcess needs"
    `GeometricEquations` validates the length of a `GridProcess` when the problem is built, and
    rejects one that is too short. The number it requires is *not* simply `(t₁ - t₀) / Δt`: it is
    `div(t₁ - t₀, Δt, RoundUp)`, which uses exact arithmetic rather than the rounded
    floating-point division, and the two disagree whenever the step is not exactly representable
    in binary.

    The Kubo defaults are exactly such a case. `0.1 / 0.01` evaluates to `10.0`, but `0.1` is a
    shade above one tenth, so the requirement is **eleven** columns for a run that takes ten
    steps. Derive the length rather than assuming it:

    ```@example noise
    tspan, Δt, m = (0.0, 0.1), 0.01, 1

    nw = Int(div(tspan[end] - tspan[begin], Δt, RoundUp))
    gp = GridProcess(randn(m, nw) .* sqrt(Δt))

    nw, Int((tspan[end] - tspan[begin]) / Δt)     # eleven, against the assumed ten
    ```

    A surplus column is harmless — the integrator consumes one per step and ignores the rest.

**Reproducibility.** The same increments give the same trajectory, whatever the state of any
random number generator.

**Comparing methods.** Two schemes can only be compared on a common sample path; run against
independently drawn noise they differ by the noise, not by the method.

```@example noise
import GeometricProblems.KuboOscillator as Kubo

tspan = (0.0, Δt * 100)
nw = Int(div(tspan[end] - tspan[begin], Δt, RoundUp))   # derived, as in the warning above

ΔW = randn(1, nw) .* sqrt(Δt)
prob = SDEProblem(Kubo.kubo_oscillator_sde_v, Kubo.kubo_oscillator_sde_B, GridProcess(ΔW),
                  tspan, Δt, [0.5, 0.0]; parameters = Kubo.default_parameters())

a = integrate(prob, BurrageE1())
b = integrate(prob, StochasticGLRK(1))    # same path, so the difference is the method

maximum(abs, a.q[end] - b.q[end])
```

**Measuring strong convergence.** Coarsening one fine-grid path gives every refinement the same
realisation, which is what makes the error a pathwise error. `scripts/convergence_order.jl` does
this, and computes the matching ``\Delta Z`` by quadrature over the fine subintervals — leaving
them zero would silently cost a scheme that uses them the order those terms buy.

**The deterministic limit.** With zero increments an SDE is its own drift, so a stochastic method
must reproduce the deterministic method it is built from, up to round-off. This is the sharpest
available check on the drift half of a scheme:

```@example noise
using GeometricIntegrators: Gauss

zeronoise = GridProcess(zeros(1, nw), zeros(1, nw))
sde0 = SDEProblem(Kubo.kubo_oscillator_sde_v, Kubo.kubo_oscillator_sde_B, zeronoise,
                  tspan, Δt, [0.5, 0.0]; parameters = Kubo.default_parameters())
ode0 = odeproblem(; timespan = tspan, timestep = Δt)

maximum(abs, integrate(sde0, StochasticGLRK(1)).q[end] - integrate(ode0, Gauss(1)).q[end])
```

The two are not bit-identical, and are not expected to be: the stochastic update adds a diffusion
term that happens to be zero, so it accumulates its increment in a different order from the
deterministic one and the results may differ in the last bit. Agreement at round-off is the claim;
`==` would be asserting something about floating-point summation order rather than about the
scheme.

## Random number generators

Each method carries its own generator, given at construction:

```@example noise
using Random
integrate(sdeproblem(), BurrageE1(; rng = Xoshiro(42))).q[end]
```

The default is `Random.default_rng()`. Two integrators built from the *same* method object share
its generator and therefore its stream; build the method twice for two independent streams, or
pass a seeded generator for a reproducible one. A `GridProcess` ignores the generator entirely.
