"""
Tableau for the 1-stage SRKw1 method due to Wang, Hong & Xu
  Method cited in
  Wang, Hong, Xu, "Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems",
  Commun. Comput. Phys. 21(1), 2017
  According to the paper, the method has weak order 1.0.
"""
function TableauSRKw1(x::Number = 0.0)
    A0 = 0.5*ones(typeof(x), 1, 1)

    A1 = (1.0 - x)*ones(typeof(x), 1, 1)

    B0 = x*ones(typeof(x), 1, 1)

    B1 = 0.5*ones(typeof(x), 1, 1)

    B3 = 0.5*ones(typeof(x), 1, 1)

    α = [1.0]
    β1 = [1.0]

    c0 = [0.5]
    c1 = [1.0 - x]

    TableauWIRK(:SRKw1, A0, A1, B0, B1, B3, α, β1, c0, c1)
end

"""
Tableau for the 4-stage SRKw2 method due to Wang, Hong & Xu
  Method cited in
  Wang, Hong, Xu, "Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems",
  Commun. Comput. Phys. 21(1), 2017
  According to the paper, the method has weak order 2.0 when applied to systems
  driven by one-dimensional noise.
"""
function TableauSRKw2(x1::Number = 0.0, x2::Number = 0.0, x3::Number = 0.0)
    A0 = [[1.0 / 8.0 0.0 0.0 0.0]
          [1.0 / 4.0 1.0 / 8.0 0.0 0.0]
          [1.0 / 4.0 1.0 / 4.0 1.0 / 8.0 0.0]
          [1.0 / 4.0 1.0 / 4.0 1.0 / 4.0 1.0 / 8.0]]

    A1 = [[-1.0 / 6.0 + sqrt(3)/6.0 1.0 / 3.0 - sqrt(3)/6.0 0.0 1.0 / 3.0]
          [1.0 / 2.0 0.0 0.0 0.0]
          [0.0 0.0 0.0 0.0]
          [0.0 0.0 0.0 0.0]]

    B0 = [[5.0 / 6.0 - sqrt(3)/3.0 -1.0 / 2.0 0.0 0.0]
          [-1.0 / 6.0 + sqrt(3)/3.0 1.0 / 2.0 0.0 0.0]
          [1.0 / 2.0 1.0 / 2.0 0.0 0.0]
          [-1.0 / 6.0 1.0 / 2.0 0.0 0.0]]

    B1 = [[1.0 / 4.0 1.0 / 4.0 - sqrt(3)/6.0 0.0 0.0]
          [1.0 / 4.0 + sqrt(3)/6.0 1.0 / 4.0 0.0 0.0]
          [x1 x2 0.0 x3]
          [0.0 -1.0 / 2.0 0.0 0.0]]

    B3 = zeros(typeof(x1), 4, 4)

    α = [0.25, 0.25, 0.25, 0.25]
    β1 = [0.5, 0.5, 0.0, 0.0]

    c0 = [1.0 / 8.0, 3.0 / 8.0, 5.0 / 8.0, 7.0 / 8.0]
    c1 = [0.5, 0.5, 0.0, 0.0]

    TableauWIRK(:SRKw2, A0, A1, B0, B1, B3, α, β1, c0, c1)
end

#
# Named methods
#

@doc raw"""
    SRKw1(λ = 0.0; rng = Random.default_rng())

One-stage weak implicit Runge-Kutta method SRKw1 of Wang, Hong & Xu, *Construction of Symplectic
Runge-Kutta Methods for Stochastic Hamiltonian Systems*, Commun. Comput. Phys. 21(1), 2017. Weak
order 1.0.

`λ` is a free parameter of the family; `λ = 1/2` reduces the method to the implicit midpoint rule
with the Wiener increments replaced by the discrete weak variables.

Reproduced in Kraus & Tyranowski §3.4.3, where it is shown that the symplecticity conditions of
Wang, Hong & Xu are equivalent to the weak Lagrange-d'Alembert conditions — so this is a
variational integrator for forced systems too.
"""
SRKw1(λ::Number = 0.0; kwargs...) = WIRK(TableauSRKw1(λ); kwargs...)

@doc raw"""
    SRKw2(λ₁ = 0.0, λ₂ = 0.0, λ₃ = 0.0; rng = Random.default_rng())

Four-stage weak implicit Runge-Kutta method SRKw2 of Wang, Hong & Xu (2017). Weak order 2.0 for
systems driven by one-dimensional noise.

The three parameters are free and, as the paper notes, have no effect on ``q_{n+1}`` or
``p_{n+1}``; they default to zero for convenience. Of the weak methods here this is the most
accurate.
"""
function SRKw2(λ₁::Number = 0.0, λ₂::Number = 0.0, λ₃::Number = 0.0; kwargs...)
    WIRK(TableauSRKw2(λ₁, λ₂, λ₃); kwargs...)
end
