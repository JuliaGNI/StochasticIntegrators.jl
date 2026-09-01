"""
Tableau for the stochastic symplectic Euler method
  Tableau for the stochastic symplectic Euler method
  Satisfies the conditions for Lagrange-d'Alembert integrators.
  Satisfies the conditions for strong convergence of order 1.0 for one Wiener process
  for special choices of the stochastic Hamiltonians and forces, e.g., h=h(q), f=0.
"""
function TableauStochasticSymplecticEuler()
    a_q = ones(Float64, 1, 1)
    b_q = [1.0]
    c_q = [1.0]

    a_p = zeros(Float64, 1, 1)
    b_p = [1.0]
    c_p = [0.0]

    TableauSIPRK(:StochasticSymplecticEuler, 1, a_q, b_q, c_q, 1,
        a_q, b_q, c_q, 1, a_p, b_p, c_p, 1, a_p, b_p, c_p)
end

"Tableau for the 2-stage stochastic LobattoIIA-IIB method (Stormer-Verlet)"
function TableauStochasticStoermerVerlet()
    TableauSIPRK(:StochasticStoermerVerlet,
        TableauLobattoIIIA(2), TableauLobattoIIIA(2),
        TableauLobattoIIIB(2), TableauLobattoIIIB(2))
end

#
# Named methods
#

@doc raw"""
    StochasticSymplecticEuler(; K = 0, rng = Random.default_rng())

One-stage stochastic symplectic Euler method — implicit in ``q``, explicit in ``p``.

Satisfies all eight Lagrange-d'Alembert conditions, so it is a variational integrator. Its
mean-square order conditions, however, are **not** satisfied in general:
``\beta^T b e = 1 \neq 1/2``. It attains strong order 1.0 only for special choices of the
stochastic Hamiltonians and forces — for instance ``h = h(q)`` with ``f = 0``, for which every term
of the order conditions carries a vanishing factor and the conditions hold trivially. Outside such
cases expect order 0.5.

Verified against the conditions of Kraus & Tyranowski by `scripts/tableau_conditions.jl`; the
method itself does not appear in that paper.
"""
function StochasticSymplecticEuler(; kwargs...)
    SIPRK(TableauStochasticSymplecticEuler(); kwargs...)
end

@doc raw"""
    StochasticStoermerVerlet(; K = 0, rng = Random.default_rng())

Two-stage stochastic Störmer-Verlet method, built from the Lobatto IIIA-IIIB pair.

This is the scheme of Kraus & Tyranowski §3.3.3, Example 2, and its tableau matches theirs
exactly: ``(a, b) =`` Lobatto IIIA(2), ``(\bar a, \bar b) =`` Lobatto IIIB(2), all weights
``(1/2, 1/2)``. It satisfies both the Lagrange-d'Alembert conditions and the conditions for
mean-square order 1.0.

It is the efficient choice among the strong partitioned methods. The stage equations decouple: the
first is solved on its own, the second becomes explicit when the Hamiltonians are separable, and
the last is an explicit update. If the forcing is independent of ``p`` the whole method is
explicit.
"""
function StochasticStoermerVerlet(; kwargs...)
    SIPRK(TableauStochasticStoermerVerlet(); kwargs...)
end
