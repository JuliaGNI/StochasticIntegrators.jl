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

### Compatibility

- Floors are the versions that first carry what this package needs, not the current releases:
  `GeometricBase 0.14.10` (`noise`, `noisedims`), `GeometricEquations 0.21.3` (the concrete noise
  processes) and `GeometricIntegratorsBase 0.6.5` (`default_extrapolation(method)`). CI's `min`
  matrix entry resolves exactly these, so a floor left at the major version would have tested a
  combination that cannot work. `GeometricProblems 0.9` is bounded for the same reason — the test
  suite needs the Kubo oscillator on a real noise process.

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

- `test/multidimensional_noise_tests.jl` covers the case of more than one driving Wiener process,
  which nothing else in the suite reaches. The branches that distinguish one noise dimension from
  another — the `qdiff3` block of `WERK` and `WIRK` — were previously never executed.

  Running a weak method on a two-dimensional problem exercises that branch, but on its own would
  not show that it *does* anything: a scheme ignoring `qdiff3` would still produce finite,
  energy-conserving output. So each method is run twice on one prescribed path, once with its own
  tableau and once with `qdiff3` zeroed. With two noise dimensions the two must disagree; with one
  they must agree exactly. Together those pin down both that the branch runs and when.

- Documentation covering the theory, the noise processes, each method family and the
  implementation.

### Performance

- **Each tableau type now carries one type parameter per Butcher tableau it holds.**
  `RungeKutta.Tableau` is `Tableau{T,S,R∞,L}`, so the previous `::Tableau{T}` field annotations
  were abstract: every `tab.qdrift.a[i,j]` in a stage loop inferred as `Any` and materialised the
  whole `SMatrix` on the heap. Measured per step, on the Kubo problems:

  | method | before | after |
  |:--|--:|--:|
  | SIRK | 8880 B | 80 B |
  | WIRK | 18576 B | 80 B |
  | SIPRK | 44528 B | 80 B |
  | SISPRK | 64944 B | 80 B |
  | SERK | 9456 B | 2896 B |
  | WERK | 15584 B | 6096 B |

  The tableaus within one method need not share a concrete type — `R∞` is `Missing` for a tableau
  built from bare coefficients and a number for one of `RungeKutta`'s named families — so a single
  shared parameter would not have been enough. The constructors infer all of them, so nothing
  changes at a call site.

  The two explicit methods keep a residual cost that this does not explain and that is **not** in
  the tableau access: `sol.q` is a `StateWithError`, and indexing one goes through `parent`, which
  allocates 240 B per call in `GeometricBase`. Hoisting it, a function barrier around the stage
  loop, and hoisting the coefficient arrays each left the figure unchanged, so the cause is
  elsewhere and is left open.

- Dropped the redundant `mul!` in the `SIRK` and `SIPRK` initial-guess predictors, which recomputed
  `B1 * Δw` and `G1 * Δw` a few lines after `ΔQ` and `ΔP` already held them, and the four cache
  fields that held the duplicates.

### Bug Fixes

- **Compensated summation was silently dropped.** The final updates added their increment with
  `q .+= Δq`. `sol.q` is a `StateWithError`, which carries a running round-off error beside the
  state — but broadcasting falls through to the generic `AbstractArray` path, writes through
  `setindex!` and never touches the error field, so no compensation happened at all. The
  pre-rewrite `update!(sol, y, k)` did compensate, so this was a regression in the rewrite, not a
  pre-existing gap. The updates now call `GeometricBase.add!`, which is the operation that
  performs it.

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

- **`ntime(problem)` and `ntime(solution)` disagree by one for some time steps.** For the Kubo
  defaults — `timespan = (0.0, 0.1)`, `timestep = 0.01` — `GeometricEquations.ntimesteps` computes
  `div(0.1, 0.01, RoundUp) = 11` while the run takes ten steps and ends correctly at `t = 0.1`,
  so `ntime(problem) == 11` and `ntime(solution) == 10`. `0.1 / 0.01` is exactly `10.0` in
  floating point; `div(…, RoundUp)` uses the exact values, and `0.1` is a shade above one tenth.

  This is upstream and predates the noise work, but `GridProcess` is what exposes it: a process is
  validated against `ntimesteps`, so it must carry eleven increments for a ten-step run. The tests
  and `docs/src/noise.md` derive the length the same way rather than assuming `(t₁ - t₀) / Δt`.

- `SDEEnsemble`, `PSDEEnsemble` and `SPSDEEnsemble` have no convenience constructors upstream, so
  an ensemble of stochastic problems has to be assembled by hand from the equation.

- The two explicit methods still allocate per step (see *Performance*), for a reason that is not
  the tableau field types and that this change does not identify.
