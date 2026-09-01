@doc raw"""
    TableauSISPRK{T}

Tableau of a stochastic implicit **split** partitioned Runge-Kutta method.

Six Butcher tableaus. `qdrift`, `qdiff` act on the ``q`` equation; the ``p`` equation is split, with
`pdrift1`, `pdiff1` acting on the Hamiltonian part and `pdrift2`, `pdiff2` on the forcing part.
In the notation of Kraus & Tyranowski, *Variational integrators for stochastic dissipative
Hamiltonian systems*, Eq. (SPRK):

| field | paper | acts on |
|:--|:--|:--|
| `qdrift`  | ``(a, \alpha)``            | ``\partial H / \partial p`` |
| `qdiff`   | ``(b, \beta)``             | ``\partial h_r / \partial p`` |
| `pdrift1` | ``(\bar a, \alpha)``       | ``-\partial H / \partial q`` |
| `pdiff1`  | ``(\bar b, \beta)``        | ``-\partial h_r / \partial q`` |
| `pdrift2` | ``(\hat a, \hat \alpha)``  | ``F`` |
| `pdiff2`  | ``(\hat b, \hat \beta)``   | ``f_r`` |

The paper's tableau shares weight vectors: ``\alpha`` is used for both `qdrift` and `pdrift1`, and
``\beta`` for both `qdiff` and `pdiff1`; only the forcing tableaus carry independent weights
``\hat\alpha``, ``\hat\beta``. This is a consequence of the Lagrange-d'Alembert conditions, not a
convention, so the constructor checks it — a tableau that violates it cannot be a variational
integrator whatever else it satisfies.
"""
struct TableauSISPRK{T} <: AbstractTableau{T}
    name::Symbol
    s::Int
    qdrift::Tableau{T}
    qdiff::Tableau{T}
    pdrift1::Tableau{T}
    pdrift2::Tableau{T}
    pdiff1::Tableau{T}
    pdiff2::Tableau{T}

    function TableauSISPRK{T}(name, s, qdrift, qdiff, pdrift1, pdrift2,
            pdiff1, pdiff2) where {T}
        @assert s == qdrift.s == qdiff.s == pdrift1.s == pdrift2.s == pdiff1.s == pdiff2.s
        @assert qdrift.b==pdrift1.b "the drift weights α must be shared by the q equation and " *
                                    "the Hamiltonian part of the p equation"
        @assert qdiff.b==pdiff1.b "the diffusion weights β must be shared by the q equation and " *
                                  "the Hamiltonian part of the p equation"
        new(name, s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2)
    end
end

function TableauSISPRK(name::Symbol, qdrift::Tableau{T}, qdiff::Tableau{T},
        pdrift1::Tableau{T}, pdrift2::Tableau{T},
        pdiff1::Tableau{T}, pdiff2::Tableau{T}) where {T}
    TableauSISPRK{T}(name, qdrift.s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2)
end

function TableauSISPRK(name::Symbol, qorder_drift::Int, qa_drift::AbstractMatrix{T},
        qb_drift::AbstractVector{T}, qc_drift::AbstractVector{T},
        qorder_diff::Int, qa_diff::AbstractMatrix{T}, qb_diff::AbstractVector{T},
        qc_diff::AbstractVector{T},
        porder1_drift::Int, pa1_drift::AbstractMatrix{T}, pb1_drift::AbstractVector{T},
        pc1_drift::AbstractVector{T},
        porder2_drift::Int, pa2_drift::AbstractMatrix{T}, pb2_drift::AbstractVector{T},
        pc2_drift::AbstractVector{T},
        porder1_diff::Int, pa1_diff::AbstractMatrix{T}, pb1_diff::AbstractVector{T},
        pc1_diff::AbstractVector{T},
        porder2_diff::Int, pa2_diff::AbstractMatrix{T}, pb2_diff::AbstractVector{T},
        pc2_diff::AbstractVector{T}) where {T}
    TableauSISPRK{T}(name, length(qc_drift),
        Tableau(name, qorder_drift, qa_drift, qb_drift, qc_drift),
        Tableau(name, qorder_diff, qa_diff, qb_diff, qc_diff),
        Tableau(name, porder1_drift, pa1_drift, pb1_drift, pc1_drift),
        Tableau(name, porder2_drift, pa2_drift, pb2_drift, pc2_drift),
        Tableau(name, porder1_diff, pa1_diff, pb1_diff, pc1_diff),
        Tableau(name, porder2_diff, pa2_diff, pb2_diff, pc2_diff))
end

@doc raw"""
    SISPRK(tableau; K = 0, rng = Random.default_rng())

Stochastic implicit **split** partitioned Runge-Kutta method for a `SPSDEProblem`, of strong
(mean-square) convergence.

This is the method of Kraus & Tyranowski for stochastic **forced** (dissipative) Hamiltonian
systems
```math
\begin{aligned}
d q &= \frac{\partial H}{\partial p} dt + \sum_r \frac{\partial h_r}{\partial p} \circ dW^r , \\
d p &= \left[ -\frac{\partial H}{\partial q} + F \right] dt
     + \sum_r \left[ -\frac{\partial h_r}{\partial q} + f_r \right] \circ dW^r .
\end{aligned}
```
The point of the split is that the Hamiltonian and the forcing terms of the ``p`` equation get
**different** coefficients — ``(\bar a, \bar b)`` against ``(\hat a, \hat b)``. That extra freedom
is what allows the scheme to satisfy the Lagrange-d'Alembert conditions in the presence of
forcing, making it a discrete Lagrange-d'Alembert variational integrator rather than merely a
symplectic method with a dissipative term bolted on. It is the reason such a scheme tracks the
correct energy *decay* of a damped system over long integrations instead of accumulating its own
spurious drift.

Constructed through [`StochasticLobattoIIIABD2`](@ref) or
[`ModifiedStochasticStoermerVerlet`](@ref).
"""
struct SISPRK{TT, RNG} <: SPSDEMethod
    tableau::TableauSISPRK{TT}
    K::Int
    rng::RNG
end

SISPRK(tableau::TableauSISPRK; K = 0, rng = Random.default_rng()) = SISPRK(tableau, K, rng)

convergence(::SISPRK) = :strong
truncation(method::SISPRK) = method.K
GeometricIntegratorsBase.isexplicit(::SISPRK) = false
GeometricIntegratorsBase.isimplicit(::SISPRK) = true

"Cache of a stochastic implicit split partitioned Runge-Kutta method."
struct SISPRKCache{DT, D, M, S} <: PSDEIntegratorCache{DT}
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
    F1::Vector{Vector{DT}}
    F2::Vector{Vector{DT}}
    B::Vector{Matrix{DT}}
    G1::Vector{Matrix{DT}}
    G2::Vector{Matrix{DT}}
    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}

    function SISPRKCache{DT, D, M, S}() where {DT, D, M, S}
        new(zeros(DT, 2 * D * S),
            zeros(DT, M), zeros(DT, M), zeros(DT, D), zeros(DT, D),
            zeros(DT, M), zeros(DT, M),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_matrix(DT, D, M, S),
            create_internal_stage_matrix(DT, D, M, S),
            create_internal_stage_matrix(DT, D, M, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S))
    end
end

GeometricIntegratorsBase.nlsolution(cache::SISPRKCache) = cache.x

function GeometricIntegratorsBase.Cache{ST}(problem::AbstractProblemSPSDE,
        method::SISPRK; kwargs...) where {ST}
    SISPRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}(;
        kwargs...)
end

@inline function GeometricIntegratorsBase.CacheType(ST, problem::AbstractProblemSPSDE,
        method::SISPRK)
    SISPRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}
end

function GeometricIntegratorsBase.solversize(method::SISPRK, problem::AbstractProblemSPSDE)
    2 * length(vec(initial_conditions(problem).q)) * nstages(method)
end

GeometricIntegratorsBase.default_solver(::SISPRK) = Newton()

"""
Initial guess for the stage increments of a [`SISPRK`](@ref) method: zero.

No explicit predictor is used. The split methods are applied to forced systems whose forcing may
be stiff, and the R2 predictor that serves [`SIRK`](@ref) and [`SIPRK`](@ref) well has no such
guarantee here.
"""
function stage_initial_guess!(sol, history, params,
        int::GeometricIntegrator{<:SISPRK, <:AbstractProblemSPSDE})
    nlsolution(int) .= 0
end

function GeometricIntegratorsBase.components!(x::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:SISPRK, <:AbstractProblemSPSDE}) where {ST}
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

    # compute the vector fields at the internal stages
    for i in 1:S
        local tqᵢ = t̄ + Δt * tab.qdrift.c[i]
        local tpᵢ₁ = t̄ + Δt * tab.pdrift1.c[i]
        local tpᵢ₂ = t̄ + Δt * tab.pdrift2.c[i]

        equ.v(c.V[i], tqᵢ, c.Q[i], c.P[i], params)
        equ.f1(c.F1[i], tpᵢ₁, c.Q[i], c.P[i], params)
        equ.f2(c.F2[i], tpᵢ₂, c.Q[i], c.P[i], params)

        equ.B(c.B[i], tqᵢ, c.Q[i], c.P[i], params)
        equ.G1(c.G1[i], tpᵢ₁, c.Q[i], c.P[i], params)
        equ.G2(c.G2[i], tpᵢ₂, c.Q[i], c.P[i], params)
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:SISPRK, <:AbstractProblemSPSDE}) where {ST}
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
                z1 += tab.pdrift1.a[i, j] * c.F1[j][k] * Δt +
                      tab.pdrift2.a[i, j] * c.F2[j][k] * Δt
                z2 += tab.pdrift1.â[i, j] * c.F1[j][k] * Δt +
                      tab.pdrift2.â[i, j] * c.F2[j][k] * Δt
                for l in 1:M
                    y1 += tab.qdiff.a[i, j] * c.B[j][k, l] * ΔW[l]
                    y2 += tab.qdiff.â[i, j] * c.B[j][k, l] * ΔW[l]
                    z1 += tab.pdiff1.a[i, j] * c.G1[j][k, l] * ΔW[l] +
                          tab.pdiff2.a[i, j] * c.G2[j][k, l] * ΔW[l]
                    z2 += tab.pdiff1.â[i, j] * c.G1[j][k, l] * ΔW[l] +
                          tab.pdiff2.â[i, j] * c.G2[j][k, l] * ΔW[l]
                end
            end
            b[D * (i - 1) + k] = -c.Y[i][k] + (y1 + y2)
            b[D * (S + i - 1) + k] = -c.Z[i][k] + (z1 + z2)
        end
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, x::AbstractVector{ST},
        sol, params,
        int::GeometricIntegrator{<:SISPRK, <:AbstractProblemSPSDE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end

function GeometricIntegratorsBase.update!(sol, params, x::AbstractVector{DT},
        int::GeometricIntegrator{<:SISPRK, <:AbstractProblemSPSDE}) where {DT}
    local c = cache(int, DT)
    local tab = tableau(int)
    local Δt = timestep(int)

    components!(x, sol, params, int)

    update_solution!(sol.q, sol.p, c.Δq, c.Δp, c.V, c.F1, c.F2, c.B, c.G1, c.G2,
        tab.qdrift.b, tab.qdiff.b, tab.pdrift1.b, tab.pdrift2.b,
        tab.pdiff1.b, tab.pdiff2.b, Δt, c.ΔW, c.Δy, c.Δz)
    update_solution!(sol.q, sol.p, c.Δq, c.Δp, c.V, c.F1, c.F2, c.B, c.G1, c.G2,
        tab.qdrift.b̂, tab.qdiff.b̂, tab.pdrift1.b̂, tab.pdrift2.b̂,
        tab.pdiff1.b̂, tab.pdiff2.b̂, Δt, c.ΔW, c.Δy, c.Δz)
end

function GeometricIntegratorsBase.integrate_step!(sol, history, params,
        int::GeometricIntegrator{<:SISPRK, <:AbstractProblemSPSDE})
    sample_noise!(int, sol)
    stage_initial_guess!(sol, history, params, int)

    local solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int),
        (sol, params, int))
    check_solver_status(solverstatus, int)

    update!(sol, params, nlsolution(int), int)
end
