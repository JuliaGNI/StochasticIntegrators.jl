```@meta
CurrentModule = StochasticIntegrators
```

# StochasticIntegrators.jl

Geometric integrators for stochastic differential equations, built on the
[JuliaGNI](https://github.com/JuliaGNI) ecosystem.

The package implements stochastic Runge-Kutta methods for the `SDE`, `PSDE` and `SPSDE` problem
types of [GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl), plugged into
the integrator interface of
[GeometricIntegratorsBase.jl](https://github.com/JuliaGNI/GeometricIntegratorsBase.jl). Its
distinguishing feature is that most of these methods are **structure preserving**: applied to a
stochastic Hamiltonian system they are symplectic, or — for forced and dissipative systems —
Lagrange-d'Alembert variational integrators. That is what keeps a long integration from
accumulating spurious energy drift, and what makes the results usable for ergodic averages.

Much of the theory these methods rest on is developed in

> M. Kraus and T. M. Tyranowski,
> *Variational integrators for stochastic dissipative Hamiltonian systems*,
> IMA Journal of Numerical Analysis, 2021.
> [`arXiv:1904.06205`](https://arxiv.org/abs/1904.06205)

## Installation

```julia
using Pkg
Pkg.add("StochasticIntegrators")
```

## Quick start

```@example index
using StochasticIntegrators
using GeometricProblems.KuboOscillator

prob = sdeproblem()                       # a Kubo oscillator driven by 1-d noise
sol  = integrate(prob, BurrageE1())       # explicit, strong

sol.q[end]
```

Everything follows the usual JuliaGNI shape: a **problem** describes the system, a **method**
describes the scheme, and `integrate(problem, method)` runs it and returns a `GeometricSolution`.

Two things are specific to the stochastic case.

**The driving noise belongs to the problem.** A `WienerProcess(m)` in the problem's `noise` field
says the system is driven by ``m`` independent Wiener processes; the integrator draws the
increments. Replacing it with a `GridProcess` prescribes them instead, which is how a run is made
reproducible or two methods compared on one sample path — see [Noise Processes](@ref).

**Strong and weak convergence are properties of the method, not options.** A strong method tracks
individual sample paths; a weak one only reproduces expectations, and is fed cheaper increments
accordingly. Choosing the wrong kind silently gives answers of the wrong accuracy, so the package
derives it from the method rather than offering it as a keyword — see [Theory](@ref).

## Choosing a method

| you want | use |
|:--|:--|
| a cheap explicit scheme, sample paths matter | [`BurrageR2`](@ref) — two stages, strong order 1.0 |
| the best strong order per explicit step | [`BurrageCL`](@ref) — four stages, strong order 1.5 |
| accuracy per step, sample paths matter | [`StochasticGLRK`](@ref), [`StochasticDIRK`](@ref) |
| a Hamiltonian system, long integration | [`StochasticStoermerVerlet`](@ref) |
| a *forced* or damped Hamiltonian system | [`ModifiedStochasticStoermerVerlet`](@ref) |
| expectations, ergodic averages, distributions | [`SRKw2`](@ref), [`RoesslerRS1`](@ref) |

The [Methods](@ref) page gives each family in full.

## Contents

```@contents
Pages = ["theory.md", "noise.md", "methods.md", "implementation.md"]
Depth = 2
```
