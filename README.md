# StochasticIntegrators

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNI.github.io/StochasticIntegrators.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNI.github.io/StochasticIntegrators.jl/latest)
[![Build Status](https://github.com/JuliaGNI/StochasticIntegrators.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/StochasticIntegrators.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/StochasticIntegrators.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/StochasticIntegrators.jl)

StochasticIntegrators.jl is a library of geometric integrators for stochastic differential equations.

Most of these methods are structure preserving: applied to a stochastic Hamiltonian system they
are symplectic, and for forced or dissipative systems they are Lagrange-d'Alembert variational
integrators. That is what keeps a long integration from accumulating spurious energy drift. The
theory is developed in

> M. Kraus and T. M. Tyranowski,
> *Variational integrators for stochastic dissipative Hamiltonian systems*,
> IMA Journal of Numerical Analysis, 2021.
> [arXiv:1904.06205](https://arxiv.org/abs/1904.06205)


## Installation

*StochasticIntegrators.jl* and all of its dependencies can be installed via the Julia REPL by typing 
```
]add StochasticIntegrators
```

## Usage

```julia
using StochasticIntegrators
using GeometricProblems.KuboOscillator

sol = integrate(sdeproblem(), BurrageE1())          # explicit, strong
sol = integrate(psdeproblem(), StochasticStoermerVerlet())   # symplectic
sol = integrate(spsdeproblem(), ModifiedStochasticStoermerVerlet())  # forced systems
sol = integrate(sdeproblem(), SRKw2())              # weak, for expectations
```

A problem carries its own driving noise: `WienerProcess(m)` draws increments as it goes, while
`GridProcess(ΔW, ΔZ)` prescribes them, which is how a run is made reproducible or two methods
compared on one sample path. Whether a scheme needs strong or weak increments follows from the
method and is not a user option.

See the [documentation](https://JuliaGNI.github.io/StochasticIntegrators.jl/latest) for the
theory, the available methods and the implementation.

## Verification

Two scripts in `scripts/` establish the claims the package makes, and are worth running before
trusting a result:

- `tableau_conditions.jl` checks every tableau against the Lagrange-d'Alembert conditions and the
  mean-square order conditions of the paper above, asserting the expected outcome for each scheme
  — including the two that deliberately fail one set.
- `convergence_order.jl` measures the mean-square convergence order on the Kubo oscillator against
  its closed-form solution.

## Development

### Git hooks

Two hooks live in `.githooks`. They are **not active in a fresh clone** — `core.hooksPath` is local
configuration and does not travel with a push — so enable them once per clone:

```sh
git config core.hooksPath .githooks
```

**`pre-commit`** acts on **staged `.jl` files only**, and exits immediately when a commit stages
none, so a documentation- or workflow-only commit is not slowed down by it:

- **JuliaFormatter `--check`**, honouring this repository's own `.JuliaFormatter.toml` — **blocks**
  the commit. Formatting is mechanical and always fixable.
- **`fatou lint`**, when `fatou` is installed — **advisory only**, and deliberately so: its
  `unused-import` rule does not follow `include`, so it flags the load-bearing imports of every
  module file.
- **`using <Package>`**, which catches a syntax error or a broken `include` — **blocks**.

**`pre-push`** runs the full test suite with `--check-bounds=auto`, but **only when pushing to
`main` or `master`**; a topic branch is left to CI. It prints nothing for **10–30 minutes**, which
looks exactly like a network hang and is not one. If you do interrupt it, check for an orphaned
Julia process that the killed hook left behind.

Either hook can be bypassed for a single command with `--no-verify`, for a change you know it does
not apply to:

```sh
git commit --no-verify
git push --no-verify
```

The hooks are generated from one shared copy and are byte-identical across the related
repositories, so edit them there rather than here — a local edit is silently undone by the next
install.
