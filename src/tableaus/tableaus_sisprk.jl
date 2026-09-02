"""
Tableau for the 2-stage stochastic LobattoIIIA-IIIB-IIID method
  (based on the deterministic LobattoIIIA-IIIB-IIID due to L. Jay)
  It satisfies the conditions for convergence of order 1.0 for one Wiener process,
  but it doesn't satisfy the conditions for Lagrange-d'Alembert integrators
"""
function TableauStochasticLobattoIIIABD2()
    TableauSISPRK(:StochasticLobattoIIIABD2, TableauLobattoIIIA(2), TableauLobattoIIIA(2),
        TableauLobattoIIIB(2), TableauLobattoIIID(2),
        TableauLobattoIIIB(2), TableauLobattoIIID(2))
end

"""
Tableau for the 2-stage modified stochastic LobattoIIIA-IIIB method
  Satisfies the conditions for Lagrange-d'Alembert integrators
  and the conditions for convergence of order 1.0 for one Wiener process
"""
function TableauModifiedStochasticStoermerVerlet(c::Number = 0.0)
    @assert c ≥ 0.0
    @assert c ≤ 1.0

    a_drift2 = [[c 0.0]
                [c 0.0]]

    b_drift2 = [c, 1-c]

    c_drift2 = [c, c]

    TableauSISPRK(:StochasticModifiedStormerVerlet,
        TableauLobattoIIIA(2), TableauLobattoIIIA(2),
        TableauLobattoIIIB(2), Tableau(:cTableau, 1, a_drift2, b_drift2, c_drift2),
        TableauLobattoIIIB(2), TableauLobattoIIIB(2))
end

#
# Named methods
#

@doc raw"""
    StochasticLobattoIIIABD2(; K = 0, rng = Random.default_rng())

Two-stage stochastic Lobatto IIIA-IIIB-IIID method, with the Lobatto IIID pair of L. Jay in the
forcing slots ``\hat a``, ``\hat b``.

**Satisfies the conditions for mean-square order 1.0 for one Wiener process, but not the
Lagrange-d'Alembert conditions** — so it is not a variational integrator. Concretely, conditions
(6) and (8) of Kraus & Tyranowski fail at ``i = j = 1``: with ``\hat b =`` Lobatto IIID(2)
``= \begin{pmatrix} 1/4 & -1/4 \\ 3/4 & 1/4 \end{pmatrix}``,
```math
\beta_i \hat b_{ij} + \hat\beta_j b_{ji} = \tfrac{1}{2}\cdot\tfrac{1}{4} + \tfrac{1}{2}\cdot 0
= \tfrac{1}{8} \neq \tfrac{1}{4} = \beta_i \hat\beta_j .
```
Conditions (1)–(4), which involve only the IIIA/IIIB pair, do hold. This is checked in
`scripts/tableau_conditions.jl`, which asserts the failure rather than merely tolerating it.

Included for comparison: it shows what is lost when the forcing tableaus are chosen for
deterministic accuracy rather than for the variational structure.
"""
function StochasticLobattoIIIABD2(; kwargs...)
    SISPRK(TableauStochasticLobattoIIIABD2(); kwargs...)
end

@doc raw"""
    ModifiedStochasticStoermerVerlet(c = 0.0; K = 0, rng = Random.default_rng())

Two-stage modified stochastic Störmer-Verlet method: the Lobatto IIIA-IIIB pair with a
one-parameter family
```math
\hat a = \begin{pmatrix} c & 0 \\ c & 0 \end{pmatrix}, \qquad \hat\alpha = (c, \; 1-c)
```
in the forcing drift slot.

Satisfies **both** the Lagrange-d'Alembert conditions and the conditions for mean-square order
1.0. Neither set constrains `c` at all: ``\hat a`` and ``\hat\alpha`` enter the order conditions
only through ``\sum_i \hat\alpha_i = c + (1-c) = 1``, and the symplecticity residuals are
identically zero in `c`. The constructor nevertheless restricts `c` to ``[0,1]``, which keeps the
collocation nodes ``(c, c)`` inside the time step — a restriction on the method rather than one
the conditions impose.

At ``c = 1/2`` it reduces to the plain stochastic Störmer-Verlet method of Kraus & Tyranowski, of
which this is a generalisation to forced systems with a free parameter.

All of the above is checked by `scripts/tableau_conditions.jl`.
"""
function ModifiedStochasticStoermerVerlet(c::Number = 0.0; kwargs...)
    SISPRK(TableauModifiedStochasticStoermerVerlet(c); kwargs...)
end
