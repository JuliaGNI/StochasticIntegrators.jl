# Release Notes

All notable changes to StochasticIntegrators.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and deliberately holds no entries. 2 versions were
released before it, the most recent `v0.2.0`, and neither is written up here: the record of
that history is `git log` and the tags. It is named as a gap rather than reconstructed,
because a changelog assembled after the fact loses exactly the reasoning that makes it worth
keeping. The `[Unreleased]` target below is provisional — confirm it when the first entry is
written.

## [Unreleased] — targeting 0.3.0

Rewritten against the current JuliaGNI ecosystem. The package had not loaded for some time: it
imported `GeometricIntegrators.Integrators` and `GeometricIntegrators.Solutions`, submodules that
no longer exist, and its solution layer had been superseded by `GeometricSolutions`. It now targets
`GeometricIntegratorsBase 0.6`, `GeometricEquations 0.21`, `GeometricSolutions 0.6`,
`GeometricBase 0.14`, `RungeKutta 0.6` and `SimpleSolvers 0.13`, and depends on
`GeometricIntegratorsBase` rather than the whole of `GeometricIntegrators`.

The numerics are carried over unchanged. All six schemes reproduce the energy-error tolerances the
previous test suite asserted, and a stochastic method driven by zero noise reproduces the
deterministic method it is built from bit for bit.

### Breaking Changes

- **`Integrator(equation, tableau, Δt)` is replaced by `integrate(problem, method)`.** The time
  step is part of the problem now, and the integrator is built from a problem and a method rather
  than from an equation and a tableau.

  ```julia
  # before
  int = Integrator(sde, TableauBurrageE1(), Δt)
  sol = Solution(sde, Δt, nt, conv = :strong)
  integrate!(int, sol)

  # now
  sol = integrate(sdeproblem(), BurrageE1())
  ```

- **Methods are named objects, and the `Tableau*` types are the layer below them.** `BurrageE1()`
  is a `SERK`, `StochasticGLRK(1)` a `SIRK`, and so on, matching `Gauss(2)` against
  `TableauGauss(2)` in the rest of the ecosystem. The tableau constructors keep their names and
  are still exported, so `SERK(TableauBurrageE1())` is the long form of `BurrageE1()`.

- **The bespoke solution layer is gone**: `SolutionSDE`, `SolutionPSDE`, `AtomicSolutionSDE`,
  `AtomicSolutionPSDE`, their constructors, and the HDF5 reading and writing. `integrate` returns
  a `GeometricSolution` from `GeometricSolutions`, which already covers these problems, and
  compensated summation now comes from `GeometricBase`'s `StateWithError`. The `HDF5` dependency
  is dropped with it; nothing upstream uses HDF5 any more, and the round-trip tests for it had
  been commented out as `# TODO Reactivate!` for some time.

- **`conv = :strong` / `:weak` is no longer a user option.** Which increments a scheme needs is a
  property of the scheme — SERK, SIRK, SIPRK and SISPRK need strong ones, WERK and WIRK weak ones
  — and mismatching them produces plausible numbers of the wrong accuracy, silently. It is derived
  from the method through `convergence(method)`, and there is no way to ask for the wrong one.

- **The noise process comes from the problem.** `WienerProcess(m)` and `GridProcess(ΔW, ΔZ)` live
  in `GeometricEquations` and sit in the equation's `noise` field; `noisedims(problem)` is what
  sizes the increment vectors. The old array-based `WienerProcess`, which held every increment of
  a whole run and was constructed alongside the solution, is gone.

- `K`, the Milstein-Tretyakov truncation parameter, moves from an integrator keyword to a method
  keyword: `StochasticGLRK(1; K = 1)`.

### New Features

- Each method carries its own random number generator, settable at construction —
  `BurrageE1(; rng = Xoshiro(42))` — so a run can be made reproducible without touching global
  state.

- `GridProcess` support throughout, which makes three things possible that were awkward or
  impossible before: driving two different methods along one and the same sample path, measuring
  strong convergence order against a pathwise-exact reference, and reducing a stochastic problem
  to its deterministic drift by prescribing zero increments.

- `scripts/tableau_conditions.jl` checks every tableau against both the Lagrange-d'Alembert
  conditions and the mean-square order conditions of Kraus & Tyranowski, *Variational integrators
  for stochastic dissipative Hamiltonian systems*. It asserts the expected outcome for each
  scheme, the two deliberate failures included: `StochasticLobattoIIIABD2` satisfies the order
  conditions but fails Lagrange-d'Alembert, and `StochasticSymplecticEuler` the reverse.

- `scripts/convergence_order.jl` measures the mean-square order on the Kubo oscillator against its
  closed-form solution along a prescribed path.

- Documentation covering the theory, the noise processes, each method family and the
  implementation.

### Bug Fixes

- **The SISPRK residual dropped its `f2` and `G2` terms.** Four lines in the old
  `function_stages!` continued an expression onto a new line beginning with `+`; Julia read the
  first line as a complete statement and the second as a unary plus whose value was discarded, so
  the forcing contributions were silently missing from the nonlinear system while the final update
  included them — an inconsistent scheme. Nothing could see it, because the only split problem
  available had `f2 ≡ 0` and `G2 ≡ 0`.

  Fixed, and covered by a regression test: the damped Kubo oscillator in split and unsplit form
  describes the same dynamics, so integrating both along one sample path must give the same
  trajectory.

## Open Issues

- **`BurrageE1` measures mean-square order 1.0 where its reference claims 1.5.** It carries a
  second diffusion tableau against the iterated integrals `ΔZ`, which is what should lift it past
  the order-1.0 ceiling. Driven by iterated integrals computed from the sample path,
  `scripts/convergence_order.jl` measures 1.05, and its error is within a few percent of
  `StochasticGLRK(1)` at every refinement — which says both are limited by the same term rather
  than one being half an order better.

  This is not a regression: the arithmetic of the `ΔZ` terms is carried over unchanged and the
  tableau coefficients are untouched. It is either a property of this test problem that the
  published claim does not cover, or a defect predating this work. `BurrageCL` and `BurrageG5`
  carry the same terms and are equally unverified.

- `SDEEnsemble`, `PSDEEnsemble` and `SPSDEEnsemble` have no convenience constructors upstream, so
  an ensemble of stochastic problems has to be assembled by hand from the equation.
