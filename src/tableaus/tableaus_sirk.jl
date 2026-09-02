"Tableau for the s-stage Gauss-Lobatto SFIRK method"
function TableauStochasticGLRK(s::Int)
    TableauSIRK(Symbol("StochasticGLRK" * string(s)), TableauGauss(s), TableauGauss(s))
end

"""
Tableau for the 2-stage stochastic symplectic DIRK method
  Tableau for the stochastic symplectic DIRK method
  Satisfies the conditions for Lagrange-d'Alembert integrators.
  Satisfies the conditions for strong convergence of order 1.0 for one Wiener process
"""
function TableauStochasticDIRK(c::Number = 0.5)
    a_drift = [[c/2 0.0]
               [c (1-c)/2]]

    b_drift = [c, 1-c]

    c_drift = [c/2, (1+c)/2]

    TableauSIRK(:StochasticSymplecticDIRK, 2, a_drift, b_drift,
        c_drift, 2, a_drift, b_drift, c_drift)
end

#
# Named methods
#

@doc raw"""
    StochasticGLRK(s; K = 0, rng = Random.default_rng())

`s`-stage stochastic Gauss-Legendre Runge-Kutta method: the Gauss collocation tableau applied to
both the drift and the diffusion.

For `s = 1` this is the stochastic implicit midpoint rule, which is the scheme Kraus & Tyranowski
list first among their examples and the simplest Lagrange-d'Alembert integrator. Because drift and
diffusion share one tableau, the symplecticity conditions reduce to the classical
``\alpha_i a_{ij} + \alpha_j a_{ji} = \alpha_i \alpha_j`` that Gauss satisfies by construction, so
the whole family is symplectic for a stochastic Hamiltonian system.

Setting the noise to zero recovers the deterministic Gauss method of the same number of stages,
which is a useful check on an implementation.
"""
StochasticGLRK(s::Int; kwargs...) = SIRK(TableauStochasticGLRK(s); kwargs...)

@doc raw"""
    StochasticDIRK(c = 0.5; K = 0, rng = Random.default_rng())

Two-stage diagonally implicit stochastic Runge-Kutta method.

This is the family of Kraus & Tyranowski, *Variational integrators for stochastic dissipative
Hamiltonian systems*, §3.3.3, where it is shown to be the **most general** two-stage DIRK
satisfying both the Lagrange-d'Alembert conditions and the conditions for mean-square order 1.0.
Every one of the six coefficient matrices is
```math
a = \begin{pmatrix} c/2 & 0 \\ c & (1-c)/2 \end{pmatrix},
\qquad \alpha = (c, \; 1-c) ,
```
with `c` free. `c = 0` and `c = 1` collapse it to the stochastic midpoint method.

Being diagonally implicit, its stages can be solved one after another rather than all at once,
which is what makes it cheaper than [`StochasticGLRK`](@ref) at equal stage count. The paper
reports it as remaining accurate at time steps where explicit schemes have long since failed.
"""
StochasticDIRK(c::Number = 0.5; kwargs...) = SIRK(TableauStochasticDIRK(c); kwargs...)
