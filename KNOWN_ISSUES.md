# Known issues

### K1 · `ntime(problem)` and `ntime(solution)` disagree by one for some time steps.

- location: `docs/src/noise.md`
- evidence: For the Kubo
  defaults — `timespan = (0.0, 0.1)`, `timestep = 0.01` — `GeometricEquations.ntimesteps` computes
  `div(0.1, 0.01, RoundUp) = 11` while the run takes ten steps and ends correctly at `t = 0.1`,
  so `ntime(problem) == 11` and `ntime(solution) == 10`. `0.1 / 0.01` is exactly `10.0` in
  floating point; `div(…, RoundUp)` uses the exact values, and `0.1` is a shade above one tenth.

  This is upstream and predates the noise work, but `GridProcess` is what exposes it: a process is
  validated against `ntimesteps`, so it must carry eleven increments for a ten-step run. The tests
  and `docs/src/noise.md` derive the length the same way rather than assuming `(t₁ - t₀) / Δt`.
- kind: upstream
- found: 2026-09-02; one issue with GeometricEquations `K2`

### K2 · `SDEEnsemble`, `PSDEEnsemble` and `SPSDEEnsemble` have no convenience constructors upstream, so an ensemble of stochastic problems has to be assembled by hand from the equation.

- location: —
- evidence: —
- kind: upstream
- found: 2026-09-02

### K3 · Aqua's `undefined_exports` fails on two names that this package does not own.

- location: `test/aqua_tests.jl`
- evidence:
  `GeometricEquations` exports `AbstractEquationDELE` and `GeometricIntegratorsBase` exports
  `initialguess!`, neither of which its own module defines; re-exporting those modules inherits
  both. The fix belongs upstream, so `test/aqua_tests.jl` passes `undefined_exports = false` and
  every other Aqua check runs.
- kind: upstream
- found: 2026-09-02

### K4 · `ExplicitImports` reports 8 explicit imports and 4 qualified accesses of names that are not `public` upstream

- location: —
- evidence: — `method`, `solver`, `solverstate`, `components!`, `residual!`, `solversize`,
  `IntegratorCache`, `CacheType`, `AbstractTableau`, `istrilstrict`. These are the ordinary
  interface of the ecosystem and are used as intended; the dependencies simply do not declare
  `public` for them yet. `check_no_stale_explicit_imports` does pass, which is the part that keeps
  Aqua's `stale_deps` honest.
- kind: upstream
- found: 2026-09-02

## Found when this file was split from the CHANGELOG

### K5 · The comment in `test/quality/aqua.jl` points to a CHANGELOG section that does not exist.

- location: `test/quality/aqua.jl:8`
- evidence: the comment reads "See `CHANGELOG.md`, *Open Issues*."; the issue it means is K3 of this file.
- kind: docs
- found: 2026-09-26

### K6 · K3 names `test/aqua_tests.jl`, which does not exist.

- location: `test/quality/aqua.jl:9`
- evidence: commit 0235d3e moved the Aqua test to `test/quality/aqua.jl`; line 9 passes `undefined_exports = false`.
- kind: docs
- found: 2026-09-26

### K7 · K4 describes an ExplicitImports check that the test suite does not run.

- location: `test/`
- evidence: `grep -rn -i explicitimports test Project.toml` prints nothing, so the claim that `check_no_stale_explicit_imports` passes is not checked.
- kind: not verified
- found: 2026-09-26
