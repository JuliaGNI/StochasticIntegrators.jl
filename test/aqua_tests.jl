using Aqua
using StochasticIntegrators
using Test

# `undefined_exports` is switched off because it cannot pass here and the cause is not in this
# package: `GeometricEquations` exports `AbstractEquationDELE` and `GeometricIntegratorsBase`
# exports `initialguess!`, neither of which its own module defines, and both are inherited through
# `@reexport`. See `CHANGELOG.md`, *Open Issues*. Everything else Aqua checks is gated here.
Aqua.test_all(StochasticIntegrators; undefined_exports = false)
