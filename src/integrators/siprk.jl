@doc raw"""
    TableauSIPRK{T}

Tableau of a stochastic implicit partitioned Runge-Kutta method: `qdrift`, `pdrift` for the drift
parts of the ``q`` and ``p`` equations, `qdiff`, `pdiff` for their diffusion parts.

In the notation of Kraus & Tyranowski, *Variational integrators for stochastic dissipative
Hamiltonian systems*, `qdrift` ``= (a, \alpha)``, `qdiff` ``= (b, \beta)``,
`pdrift` ``= (\bar a, \alpha)`` and `pdiff` ``= (\bar b, \beta)``. Note that the ``q`` and ``p``
tableaus of a pair **share their weight vectors** — only their coefficient matrices differ. That
is a requirement of the Lagrange-d'Alembert conditions, not a convention.

As elsewhere, no overall order is stored: the order of a stochastic Runge-Kutta method depends on
the noise as well as on the coefficients.
"""
struct TableauSIPRK{T} <: AbstractTableau{T}
    name::Symbol
    s::Int
    qdrift::Tableau{T}
    qdiff::Tableau{T}
    pdrift::Tableau{T}
    pdiff::Tableau{T}

    function TableauSIPRK{T}(name, s, qdrift, qdiff, pdrift, pdiff) where {T}
        @assert s == qdrift.s == qdiff.s == pdrift.s == pdiff.s
        new(name, s, qdrift, qdiff, pdrift, pdiff)
    end
end

function TableauSIPRK(name::Symbol, qdrift::Tableau{T}, qdiff::Tableau{T},
        pdrift::Tableau{T}, pdiff::Tableau{T}) where {T}
    TableauSIPRK{T}(name, qdrift.s, qdrift, qdiff, pdrift, pdiff)
end

function TableauSIPRK(name::Symbol, qorder_drift::Int, qa_drift::AbstractMatrix{T},
        qb_drift::AbstractVector{T}, qc_drift::AbstractVector{T}, qorder_diff::Int,
        qa_diff::AbstractMatrix{T}, qb_diff::AbstractVector{T}, qc_diff::AbstractVector{T},
        porder_drift::Int, pa_drift::AbstractMatrix{T}, pb_drift::AbstractVector{T},
        pc_drift::AbstractVector{T}, porder_diff::Int, pa_diff::AbstractMatrix{T},
        pb_diff::AbstractVector{T}, pc_diff::AbstractVector{T}) where {T}
    TableauSIPRK{T}(name, length(qc_drift),
        Tableau(name, qorder_drift, qa_drift, qb_drift, qc_drift),
        Tableau(name, qorder_diff, qa_diff, qb_diff, qc_diff),
        Tableau(name, porder_drift, pa_drift, pb_drift, pc_drift),
        Tableau(name, porder_diff, pa_diff, pb_diff, pc_diff))
end

@doc raw"""
    SIPRK(tableau; K = 0, rng = Random.default_rng())

Stochastic implicit **partitioned** Runge-Kutta method for a `PSDEProblem`, of strong
(mean-square) convergence.

For a stochastic Hamiltonian system
```math
\begin{aligned}
d q &= \frac{\partial H}{\partial p} \, dt + \sum_r \frac{\partial h_r}{\partial p} \circ dW^r, \\
d p &= -\frac{\partial H}{\partial q} \, dt - \sum_r \frac{\partial h_r}{\partial q} \circ dW^r,
\end{aligned}
```
a partitioned method applies different coefficients to the ``q`` and ``p`` equations, which is
what makes it possible to be symplectic without being fully implicit in both. When the
coefficients satisfy the Lagrange-d'Alembert conditions
```math
\alpha_i \bar a_{ij} + \alpha_j a_{ji} = \alpha_i \alpha_j, \qquad
\beta_i \bar b_{ij} + \beta_j b_{ji} = \beta_i \beta_j, \qquad \ldots
```
the scheme is a variational integrator: it arises from a discrete Lagrange-d'Alembert principle
and preserves the symplectic structure pathwise. In practice that is what stops the energy of an
oscillator from drifting over long integrations, which is the whole reason to prefer these methods
to a general-purpose stochastic solver.

The solve is in ``2 D S`` unknowns — the stage increments of both ``q`` and ``p``.

Constructed through [`StochasticSymplecticEuler`](@ref) or [`StochasticStoermerVerlet`](@ref).
"""
struct SIPRK{TT, RNG} <: PSDEMethod
    tableau::TableauSIPRK{TT}
    K::Int
    rng::RNG
end

SIPRK(tableau::TableauSIPRK; K = 0, rng = Random.default_rng()) = SIPRK(tableau, K, rng)

convergence(::SIPRK) = :strong
truncation(method::SIPRK) = method.K
GeometricIntegratorsBase.isexplicit(::SIPRK) = false
GeometricIntegratorsBase.isimplicit(::SIPRK) = true
GeometricIntegratorsBase.issymplectic(::SIPRK) = true

"Cache of a stochastic implicit partitioned Runge-Kutta method."
struct SIPRKCache{DT, D, M, S} <: PSDEIntegratorCache{DT}
    x::Vector{DT}

    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    Δq::Vector{DT}
    Δp::Vector{DT}
    Δy::Vector{DT}
    Δz::Vector{DT}

    Q::Vector{Vector{DT}}
    P::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    F::Vector{Vector{DT}}
    B::Vector{Matrix{DT}}
    G::Vector{Matrix{DT}}
    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}

    ΔQ::Vector{DT}
    ΔP::Vector{DT}
    Y1::Vector{DT}
    Y2::Vector{DT}
    Z1::Vector{DT}
    Z2::Vector{DT}
    V1::Vector{DT}
    V2::Vector{DT}
    F1::Vector{DT}
    F2::Vector{DT}
    B1::Matrix{DT}
    B2::Matrix{DT}
    G1::Matrix{DT}
    G2::Matrix{DT}
    Δw::Vector{DT}

    function SIPRKCache{DT, D, M, S}() where {DT, D, M, S}
        new(zeros(DT, 2 * D * S),
            zeros(DT, M), zeros(DT, M), zeros(DT, D), zeros(DT, D),
            zeros(DT, M), zeros(DT, M),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_matrix(DT, D, M, S),
            create_internal_stage_matrix(DT, D, M, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            zeros(DT, D), zeros(DT, D), zeros(DT, D), zeros(DT, D), zeros(DT, D),
            zeros(DT, D), zeros(DT, D), zeros(DT, D), zeros(DT, D), zeros(DT, D),
            zeros(DT, D, M), zeros(DT, D, M), zeros(DT, D, M), zeros(DT, D, M),
            zeros(DT, M))
    end
end

GeometricIntegratorsBase.nlsolution(cache::SIPRKCache) = cache.x

function GeometricIntegratorsBase.Cache{ST}(problem::AbstractProblemPSDE,
        method::SIPRK; kwargs...) where {ST}
    SIPRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}(;
        kwargs...)
end

@inline function GeometricIntegratorsBase.CacheType(ST, problem::AbstractProblemPSDE,
        method::SIPRK)
    SIPRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}
end

function GeometricIntegratorsBase.solversize(method::SIPRK, problem::AbstractProblemPSDE)
    2 * length(vec(initial_conditions(problem).q)) * nstages(method)
end

GeometricIntegratorsBase.default_solver(::SIPRK) = Newton()

"""
Initial guess for the stage increments of a [`SIPRK`](@ref) method, from the same explicit R2
predictor used for [`SIRK`](@ref), applied to both the `q` and the `p` equation.

When the two tableaus share their collocation nodes the momenta are predicted alongside the
positions in one pass; otherwise a second pass evaluates them at the `pdrift` nodes.
"""
function stage_initial_guess!(sol, history, params,
        int::GeometricIntegrator{<:SIPRK, <:AbstractProblemPSDE})
    local c = cache(int)
    local tab = tableau(int)
    local Δt = timestep(int)
    local t̄ = sol.t - Δt
    local D = length(c.Δq)
    local S = nstages(method(int))
    local equ = equations(int)

    equ.v(c.V1, t̄, sol.q, sol.p, params)
    equ.B(c.B1, t̄, sol.q, sol.p, params)
    equ.f(c.F1, t̄, sol.q, sol.p, params)
    equ.G(c.G1, t̄, sol.q, sol.p, params)

    local samenodes = tab.qdrift.c == tab.pdrift.c

    # positions at the qdrift nodes, and momenta too when the nodes coincide
    for i in eachstage(method(int))
        local Δt_local = tab.qdrift.c[i] * Δt
        c.Δw .= tab.qdrift.c[i] .* c.ΔW

        mul!(c.ΔQ, c.B1, c.Δw)
        mul!(c.ΔP, c.G1, c.Δw)
        @. c.Q[i] = sol.q + 2 // 3 * Δt_local * c.V1 + 2 // 3 * c.ΔQ
        @. c.P[i] = sol.p + 2 // 3 * Δt_local * c.F1 + 2 // 3 * c.ΔP

        local t2 = t̄ + 2 // 3 * Δt_local

        equ.v(c.V2, t2, c.Q[i], c.P[i], params)
        equ.B(c.B2, t2, c.Q[i], c.P[i], params)
        equ.f(c.F2, t2, c.Q[i], c.P[i], params)
        equ.G(c.G2, t2, c.Q[i], c.P[i], params)

        mul!(c.Y1, c.B1, c.Δw)
        mul!(c.Y2, c.B2, c.Δw)
        for j in 1:D
            c.x[(i - 1) * D + j] = Δt_local * (1 // 4 * c.V1[j] + 3 // 4 * c.V2[j]) +
                                   1 // 4 * c.Y1[j] + 3 // 4 * c.Y2[j]
        end

        if samenodes
            mul!(c.Z1, c.G1, c.Δw)
            mul!(c.Z2, c.G2, c.Δw)
            for j in 1:D
                c.x[(S + i - 1) * D + j] = Δt_local *
                                           (1 // 4 * c.F1[j] + 3 // 4 * c.F2[j]) +
                                           1 // 4 * c.Z1[j] + 3 // 4 * c.Z2[j]
            end
        end
    end

    # momenta at the pdrift nodes, when those differ from the qdrift ones
    if !samenodes
        for i in eachstage(method(int))
            local Δt_local = tab.pdrift.c[i] * Δt
            c.Δw .= tab.pdrift.c[i] .* c.ΔW

            mul!(c.ΔQ, c.B1, c.Δw)
            mul!(c.ΔP, c.G1, c.Δw)
            @. c.Q[i] = sol.q + 2 // 3 * Δt_local * c.V1 + 2 // 3 * c.ΔQ
            @. c.P[i] = sol.p + 2 // 3 * Δt_local * c.F1 + 2 // 3 * c.ΔP

            local t2 = t̄ + 2 // 3 * Δt_local

            equ.v(c.V2, t2, c.Q[i], c.P[i], params)
            equ.B(c.B2, t2, c.Q[i], c.P[i], params)
            equ.f(c.F2, t2, c.Q[i], c.P[i], params)
            equ.G(c.G2, t2, c.Q[i], c.P[i], params)

            mul!(c.Z1, c.G1, c.Δw)
            mul!(c.Z2, c.G2, c.Δw)
            for j in 1:D
                c.x[(S + i - 1) * D + j] = Δt_local *
                                           (1 // 4 * c.F1[j] + 3 // 4 * c.F2[j]) +
                                           1 // 4 * c.Z1[j] + 3 // 4 * c.Z2[j]
            end
        end
    end
end

function GeometricIntegratorsBase.components!(x::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:SIPRK, <:AbstractProblemPSDE}) where {ST}
    local c = cache(int, ST)
    local tab = tableau(int)
    local Δt = timestep(int)
    local t̄ = sol.t - Δt
    local D = length(c.Δq)
    local S = nstages(method(int))
    local equ = equations(int)

    # copy x to Y, Z and compute Q, P
    for i in 1:S
        for k in 1:D
            c.Y[i][k] = x[D * (i - 1) + k]
            c.Z[i][k] = x[D * (S + i - 1) + k]
            c.Q[i][k] = sol.q[k] + c.Y[i][k]
            c.P[i][k] = sol.p[k] + c.Z[i][k]
        end
    end

    # compute V, F, B and G at the internal stages
    for i in 1:S
        local tqᵢ = t̄ + Δt * tab.qdrift.c[i]
        local tpᵢ = t̄ + Δt * tab.pdrift.c[i]
        equ.v(c.V[i], tqᵢ, c.Q[i], c.P[i], params)
        equ.f(c.F[i], tpᵢ, c.Q[i], c.P[i], params)
        equ.B(c.B[i], tqᵢ, c.Q[i], c.P[i], params)
        equ.G(c.G[i], tpᵢ, c.Q[i], c.P[i], params)
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:SIPRK, <:AbstractProblemPSDE}) where {ST}
    local c = cache(int, ST)
    local tab = tableau(int)
    local Δt = timestep(int)
    local D = length(c.Δq)
    local M = length(c.Δy)
    local S = nstages(method(int))
    local ΔW = cache(int).ΔW

    for i in 1:S
        for k in 1:D
            local y1 = zero(ST)
            local y2 = zero(ST)
            local z1 = zero(ST)
            local z2 = zero(ST)
            for j in 1:S
                y1 += tab.qdrift.a[i, j] * c.V[j][k] * Δt
                y2 += tab.qdrift.â[i, j] * c.V[j][k] * Δt
                z1 += tab.pdrift.a[i, j] * c.F[j][k] * Δt
                z2 += tab.pdrift.â[i, j] * c.F[j][k] * Δt
                for l in 1:M
                    y1 += tab.qdiff.a[i, j] * c.B[j][k, l] * ΔW[l]
                    y2 += tab.qdiff.â[i, j] * c.B[j][k, l] * ΔW[l]
                    z1 += tab.pdiff.a[i, j] * c.G[j][k, l] * ΔW[l]
                    z2 += tab.pdiff.â[i, j] * c.G[j][k, l] * ΔW[l]
                end
            end
            b[D * (i - 1) + k] = -c.Y[i][k] + (y1 + y2)
            b[D * (S + i - 1) + k] = -c.Z[i][k] + (z1 + z2)
        end
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, x::AbstractVector{ST},
        sol, params,
        int::GeometricIntegrator{<:SIPRK, <:AbstractProblemPSDE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end

function GeometricIntegratorsBase.update!(sol, params, x::AbstractVector{DT},
        int::GeometricIntegrator{<:SIPRK, <:AbstractProblemPSDE}) where {DT}
    local c = cache(int, DT)
    local tab = tableau(int)
    local Δt = timestep(int)

    components!(x, sol, params, int)

    update_solution!(sol.q, sol.p, c.Δq, c.Δp, c.V, c.F, c.B, c.G,
        tab.qdrift.b, tab.qdiff.b, tab.pdrift.b, tab.pdiff.b, Δt, c.ΔW, c.Δy, c.Δz)
    update_solution!(sol.q, sol.p, c.Δq, c.Δp, c.V, c.F, c.B, c.G,
        tab.qdrift.b̂, tab.qdiff.b̂, tab.pdrift.b̂, tab.pdiff.b̂, Δt, c.ΔW, c.Δy, c.Δz)
end

function GeometricIntegratorsBase.integrate_step!(sol, history, params,
        int::GeometricIntegrator{<:SIPRK, <:AbstractProblemPSDE})
    sample_noise!(int, sol)
    stage_initial_guess!(sol, history, params, int)

    local solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int),
        (sol, params, int))
    check_solver_status(solverstatus, int)

    update!(sol, params, nlsolution(int), int)
end
