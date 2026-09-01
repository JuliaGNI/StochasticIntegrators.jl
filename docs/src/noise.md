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
truncate_increments!
```

A Gaussian increment is unbounded, and the nonlinear system of an implicit scheme need not be
solvable for arbitrarily large ones. Milstein and Tretyakov's remedy is to clip ``\Delta W`` at
``A = \sqrt{2 K \Delta t \, \lvert \log \Delta t \rvert}``, which leaves the mean-square order
intact because the discarded tail is of higher order than the scheme. It is off by default and
enabled per method with `K`:

```julia
integrate(prob, StochasticGLRK(1; K = 1))
```

## Prescribed paths

A `GridProcess` fixes the realisation instead of drawing it. It has four uses, and all four are
things that cannot be done with a `WienerProcess`.

**Reproducibility.** The same increments give the same trajectory, whatever the state of any
random number generator.

**Comparing methods.** Two schemes can only be compared on a common sample path; run against
independently drawn noise they differ by the noise, not by the method.

```julia
using StochasticIntegrators
using GeometricProblems.KuboOscillator
import GeometricProblems.KuboOscillator as Kubo

nt, Δt = 100, 0.01
ΔW  = randn(1, nt) .* sqrt(Δt)
prob = SDEProblem(Kubo.kubo_oscillator_sde_v, Kubo.kubo_oscillator_sde_B, GridProcess(ΔW),
                  (0.0, Δt*nt), Δt, [0.5, 0.0]; parameters = Kubo.default_parameters())

a = integrate(prob, BurrageE1())
b = integrate(prob, StochasticGLRK(1))    # same path, so the difference is the method
```

**Measuring strong convergence.** Coarsening one fine-grid path gives every refinement the same
realisation, which is what makes the error a pathwise error. `scripts/convergence_order.jl` does
this, and computes the matching ``\Delta Z`` by quadrature over the fine subintervals — leaving
them zero would silently cost a scheme that uses them the order those terms buy.

**The deterministic limit.** With zero increments an SDE is its own drift, so a stochastic method
must reproduce the deterministic method it is built from, to the last bit. This is the sharpest
available check on the drift half of a scheme:

```julia
using GeometricIntegrators: Gauss

zeronoise = GridProcess(zeros(1, nt), zeros(1, nt))
sde0 = SDEProblem(Kubo.kubo_oscillator_sde_v, Kubo.kubo_oscillator_sde_B, zeronoise,
                  (0.0, Δt*nt), Δt, [0.5, 0.0]; parameters = Kubo.default_parameters())

integrate(sde0, StochasticGLRK(1)).q[end] == integrate(odeproblem(), Gauss(1)).q[end]
```

## Random number generators

Each method carries its own generator, given at construction:

```julia
using Random
integrate(prob, BurrageE1(; rng = Xoshiro(42)))
```

The default is `Random.default_rng()`. Two integrators built from the *same* method object share
its generator and therefore its stream; build the method twice for two independent streams, or
pass a seeded generator for a reproducible one. A `GridProcess` ignores the generator entirely.
