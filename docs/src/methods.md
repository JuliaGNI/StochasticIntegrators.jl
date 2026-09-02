```@meta
CurrentModule = StochasticIntegrators
```

# Methods

Six families of stochastic Runge-Kutta methods, distinguished along three axes: which problem type
they apply to, whether they are explicit or implicit, and whether they converge in the strong or
the weak sense.

| family | problem | solve | convergence | structure of the schemes provided |
|:--|:--|:--|:--|:--|
| [`SERK`](@ref) | `SDEProblem` | explicit | strong | — |
| [`SIRK`](@ref) | `SDEProblem` | implicit | strong | symplectic |
| [`WERK`](@ref) | `SDEProblem` | explicit | weak | — |
| [`WIRK`](@ref) | `SDEProblem` | implicit | weak | symplectic |
| [`SIPRK`](@ref) | `PSDEProblem` | implicit | strong | symplectic |
| [`SISPRK`](@ref) | `SPSDEProblem` | implicit | strong | Lagrange-d'Alembert, except `StochasticLobattoIIIABD2` |

The last column is a property of the **coefficients**, not of the family: each of these types
accepts any tableau of its kind, and one that violates the conditions of [Theory](@ref) is not
symplectic whatever its type says. That is why `issymplectic` returns `missing` rather than
`true` — the claim belongs to the named schemes, and `scripts/tableau_conditions.jl` is what
establishes it for each of them.

Each family has a concrete type carrying a tableau and a random number generator, and a set of
named constructors for the individual schemes. `BurrageE1()` *is* a `SERK`; there is no separate
descriptor layer to specialise.

## The method hierarchy

```@docs
StochasticMethod
SDEMethod
PSDEMethod
SPSDEMethod
convergence
```

## Explicit strong methods

```@docs
SERK
TableauSERK
hasdiffusion2
```

Explicit methods cost one evaluation of the drift and diffusion per stage and involve no solve,
which makes them by far the cheapest per step. The price is a bounded stability region: for a
stiff drift or a large noise intensity they force a much smaller time step than the implicit
families, and none of them is structure preserving, so the energy of an oscillator drifts over a
long run.

```@docs
StochasticEuler
StochasticHeun
Platen
BurrageR2
BurrageCL
BurrageE1
BurrageG5
```

[`BurrageCL`](@ref), [`BurrageE1`](@ref) and [`BurrageG5`](@ref) carry a second diffusion tableau
applied to the iterated integrals ``\Delta Z``. Those terms are *necessary* to exceed the
mean-square order 1.0 ceiling that schemes built from ``\Delta W`` alone are subject to, but not
sufficient on their own: `BurrageCL` and `BurrageG5` reach order 1.5, while `BurrageE1` carries
the same terms and is of order 1.0. `scripts/convergence_order.jl` measures both, and the gap
between them is the only check in the package on the ``\Delta Z`` code path.

## Implicit strong methods

```@docs
SIRK
TableauSIRK
```

The stage equations are coupled and solved by Newton's method in ``D \times S`` unknowns. Note
that the unknowns are the stage *increments* ``Y_i``, not stage derivatives: the drift and the
diffusion enter with different coefficients, so unlike the deterministic case there is no single
derivative to solve for.

```@docs
StochasticGLRK
StochasticDIRK
```

## Partitioned methods

```@docs
SIPRK
TableauSIPRK
```

A partitioned method applies different coefficients to the ``q`` and ``p`` equations, which is
what makes it possible to be symplectic without being fully implicit in both. Use these for a
stochastic Hamiltonian system **without** forcing.

```@docs
StochasticSymplecticEuler
StochasticStoermerVerlet
```

## Split partitioned methods

```@docs
SISPRK
TableauSISPRK
```

For a **forced** or dissipative stochastic Hamiltonian system. The ``p`` equation is split so that
the Hamiltonian part ``-\partial H/\partial q`` and the forcing part ``F`` receive different
Runge-Kutta coefficients — the freedom that allows the Lagrange-d'Alembert conditions to be
satisfied in the presence of forcing. See [Theory](@ref).

```@docs
StochasticLobattoIIIABD2
ModifiedStochasticStoermerVerlet
```

## Weak methods

```@docs
WERK
TableauWERK
WIRK
TableauWIRK
```

Weak methods target expectations rather than sample paths. Structurally they are more elaborate
than the strong ones: beyond the drift stages they carry a separate family of internal stages for
each noise dimension, and the diffusion matrix is evaluated at a *different* stage vector for each
of its columns. That is what buys second-order weak accuracy for multi-dimensional noise, and it
is why they cost noticeably more per stage.

```@docs
RoesslerRS1
RoesslerRS2
SRKw1
SRKw2
```

The methods of Wang, Hong & Xu implemented as [`SRKw1`](@ref) and [`SRKw2`](@ref) are symplectic,
which for a stochastic Hamiltonian system is what keeps the energy from drifting over the long
runs a weak method is typically used for. Kraus & Tyranowski show their symplecticity conditions
to be equivalent to the weak Lagrange-d'Alembert conditions, so they remain variational
integrators when forcing is added.

## Usage

Every method is used the same way:

```julia
using StochasticIntegrators
using GeometricProblems.KuboOscillator

integrate(sdeproblem(),   BurrageE1())
integrate(sdeproblem(),   StochasticGLRK(1))
integrate(psdeproblem(),  StochasticStoermerVerlet())
integrate(spsdeproblem(), ModifiedStochasticStoermerVerlet())
```

Methods are matched to problem types by dispatch. Handing a `SIRK` a `PSDEProblem` raises an
error naming both, rather than silently doing something wrong:

```julia
integrate(psdeproblem(), StochasticGLRK(1))
# ERROR: ArgumentError: SIRK does not support problems of type PSDE.
```

Ensembles work without further ceremony, each member drawing its own realisation:

```julia
sol = integrate(sdeensemble(), BurrageE1())
```

Keyword arguments are passed to the method, not to `integrate`:

```julia
using Random

integrate(prob, StochasticGLRK(2; K = 1, rng = Xoshiro(42)))
```

To build an integrator once and step it yourself, or to inspect it:

```julia
int = GeometricIntegrator(sdeproblem(), BurrageE1())
sol = integrate(int)
```

## Tableaus

Each named method wraps a tableau, and the tableaus are available on their own — for inspecting
coefficients, for checking the conditions of [Theory](@ref), or for building a method with
non-default options.

```@docs
TableauStochasticEuler
TableauStochasticHeun
TableauPlaten
TableauBurrageR2
TableauBurrageCL
TableauBurrageE1
TableauBurrageG5
TableauStochasticGLRK
TableauStochasticDIRK
TableauStochasticSymplecticEuler
TableauStochasticStoermerVerlet
TableauStochasticLobattoIIIABD2
TableauModifiedStochasticStoermerVerlet
TableauRoesslerRS1
TableauRoesslerRS2
TableauSRKw1
TableauSRKw2
```

So `BurrageE1()` and `SERK(TableauBurrageE1())` are the same thing, and the second form is how to
pass a tableau built with other parameters:

```julia
SIRK(TableauStochasticDIRK(0.3); K = 1)
```
