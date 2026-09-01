@doc raw"""
    TableauWIRK{T}

Tableau of a weak implicit Runge-Kutta method, in the layout of Wang, Hong & Xu, *Construction of
Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems*, Commun. Comput. Phys. 21(1),
2017.

In the notation of that paper `qdrift0` ``= A^{(0)}``, `qdrift1` ``= A^{(1)}``,
`qdiff0` ``= B^{(0)}``, `qdiff1` ``= B^{(1)}``, `qdiff3` ``= B^{(3)}``, with weights
`qdrift0.b` ``= \alpha`` and `qdiff0.b` ``= \beta``, and nodes `qdrift0.c`, `qdrift1.c`.

Kraus & Tyranowski, §3.4, show that the symplecticity conditions of that paper are equivalent to
their weak Lagrange-d'Alembert conditions, so these methods remain variational integrators when a
forcing term is added.
"""
struct TableauWIRK{T} <: AbstractTableau{T}
    name::Symbol
    s::Int
    qdrift0::Tableau{T}
    qdrift1::Tableau{T}
    qdiff0::Tableau{T}
    qdiff1::Tableau{T}
    qdiff3::Tableau{T}

    function TableauWIRK{T}(name, s, qdrift0, qdrift1, qdiff0, qdiff1, qdiff3) where {T}
        @assert s == qdrift0.s == qdrift1.s == qdiff0.s == qdiff1.s == qdiff3.s
        new(name, s, qdrift0, qdrift1, qdiff0, qdiff1, qdiff3)
    end
end

function TableauWIRK(name::Symbol, qdrift0::Tableau{T}, qdrift1::Tableau{T},
        qdiff0::Tableau{T}, qdiff1::Tableau{T}, qdiff3::Tableau{T}) where {T}
    TableauWIRK{T}(name, qdrift0.s, qdrift0, qdrift1, qdiff0, qdiff1, qdiff3)
end

function TableauWIRK(name::Symbol, A0::AbstractMatrix{T}, A1::AbstractMatrix{T},
        B0::AbstractMatrix{T}, B1::AbstractMatrix{T}, B3::AbstractMatrix{T},
        α::AbstractVector{T}, β1::AbstractVector{T},
        c0::AbstractVector{T}, c1::AbstractVector{T}) where {T}
    TableauWIRK(name, Tableau(name, 0, A0, α, c0),
        Tableau(name, 0, A1, α, c1),
        Tableau(name, 0, B0, β1, c0),
        Tableau(name, 0, B1, β1, c1),
        Tableau(name, 0, B3, β1, c1))
end

@doc raw"""
    WIRK(tableau; rng = Random.default_rng())

**Weak** implicit Runge-Kutta method for a `SDEProblem`.

Like [`WERK`](@ref) it targets expectations rather than sample paths, and like [`SIRK`](@ref) it
solves for its stage increments. It carries ``m+1`` families of internal stages — one, ``Q^{(0)}``,
against which the drift is evaluated, and one ``Q^{(l)}`` per noise dimension, against which
column ``l`` of the diffusion is evaluated. The solver system is correspondingly larger:
``D \times S \times (1 + m)`` unknowns rather than ``D \times S``.

The methods of Wang, Hong & Xu implemented here are symplectic, which for a stochastic Hamiltonian
system is what keeps the energy from drifting over long runs — the property that makes a weak
method usable for ergodic averages in the first place.

Constructed through [`SRKw1`](@ref) or [`SRKw2`](@ref).
"""
struct WIRK{TT, RNG} <: SDEMethod
    tableau::TableauWIRK{TT}
    rng::RNG
end

WIRK(tableau::TableauWIRK; rng = Random.default_rng()) = WIRK(tableau, rng)

convergence(::WIRK) = :weak
GeometricIntegratorsBase.isexplicit(::WIRK) = false
GeometricIntegratorsBase.isimplicit(::WIRK) = true

@doc raw"""
Cache of a weak implicit Runge-Kutta method.

`Q` and `Y` are indexed `[l, i]` with `l` running from `0` to `M`: `l = 0` is the family of stages
carrying the drift, `l ≥ 1` the family against which column `l` of the diffusion is evaluated.
"""
struct WIRKCache{DT, D, M, S} <: SDEIntegratorCache{DT}
    x::Vector{DT}

    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    Δq::Vector{DT}
    Δy::Vector{DT}

    Q::Matrix{Vector{DT}}
    Y::Matrix{Vector{DT}}
    V::Vector{Vector{DT}}
    B::Vector{Matrix{DT}}

    tB::Matrix{DT}

    function WIRKCache{DT, D, M, S}() where {DT, D, M, S}
        # Q and Y are indexed [l, i] with l ∈ 0:M, stored as (M+1) × S with a one-based first
        # index; `stageindex` maps l to it.
        Q = [zeros(DT, D) for _ in 0:M, _ in 1:S]
        Y = [zeros(DT, D) for _ in 0:M, _ in 1:S]

        new(zeros(DT, D * S * (1 + M)),
            zeros(DT, M), zeros(DT, M), zeros(DT, D), zeros(DT, M),
            Q, Y,
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_matrix(DT, D, M, S),
            zeros(DT, D, M))
    end
end

"Row index in the `Q`/`Y` arrays of a `WIRKCache` for the stage family `l ∈ 0:M`."
@inline stageindex(l::Integer) = l + 1

GeometricIntegratorsBase.nlsolution(cache::WIRKCache) = cache.x

function GeometricIntegratorsBase.Cache{ST}(problem::AbstractProblemSDE,
        method::WIRK; kwargs...) where {ST}
    WIRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}(;
        kwargs...)
end

@inline function GeometricIntegratorsBase.CacheType(ST, problem::AbstractProblemSDE,
        method::WIRK)
    WIRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}
end

function GeometricIntegratorsBase.solversize(method::WIRK, problem::AbstractProblemSDE)
    length(vec(initial_conditions(problem).q)) * nstages(method) * (1 + noisedims(problem))
end

GeometricIntegratorsBase.default_solver(::WIRK) = Newton()

@doc raw"""
Initial guess for the stage increments of a [`WIRK`](@ref) method: zero.

The explicit predictor used for [`SIRK`](@ref) is deliberately *not* used here. These methods
converge only in the weak sense, so an explicit strong integrator has no reason to produce
anything close to the stage values they are solving for, and starting Newton from a confidently
wrong point is worse than starting it from zero.
"""
function stage_initial_guess!(sol, history, params,
        int::GeometricIntegrator{<:WIRK, <:AbstractProblemSDE})
    nlsolution(int) .= 0
end

function GeometricIntegratorsBase.components!(x::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:WIRK, <:AbstractProblemSDE}) where {ST}
    local c = cache(int, ST)
    local tab = tableau(int)
    local Δt = timestep(int)
    local t̄ = sol.t - Δt
    local equ = equations(int)
    local D = length(c.Δq)
    local M = length(c.Δy)
    local S = nstages(method(int))

    # copy x to Y and compute Q
    for i in 1:S
        for k in 1:D
            c.Y[stageindex(0), i][k] = x[D * (i - 1) + k]
        end
        for l in 1:M
            for k in 1:D
                c.Y[stageindex(l), i][k] = x[D * S + D * M * (i - 1) + D * (l - 1) + k]
            end
        end
        for l in 0:M
            for k in 1:D
                c.Q[stageindex(l), i][k] = sol.q[k] + c.Y[stageindex(l), i][k]
            end
        end
    end

    # compute V = v(Q⁰) and B = B(Q^l), column by column
    for i in 1:S
        equ.v(c.V[i], t̄ + Δt * tab.qdrift0.c[i], c.Q[stageindex(0), i], params)

        local t1 = t̄ + Δt * tab.qdrift1.c[i]
        for l in 1:M
            diffusion_column!(c.B[i], c.tB, l, equ, t1, c.Q[stageindex(l), i], params)
        end
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:WIRK, <:AbstractProblemSDE}) where {ST}
    local c = cache(int, ST)
    local tab = tableau(int)
    local Δt = timestep(int)
    local D = length(c.Δq)
    local M = length(c.Δy)
    local S = nstages(method(int))
    local ΔW = cache(int).ΔW

    # the terms corresponding to Y⁰
    for i in 1:S
        for k in 1:D
            local y1 = zero(ST)
            local y2 = zero(ST)
            for j in 1:S
                y1 += tab.qdrift0.a[i, j] * c.V[j][k] * Δt
                y2 += tab.qdrift0.â[i, j] * c.V[j][k] * Δt
                for l in 1:M
                    y1 += tab.qdiff0.a[i, j] * c.B[j][k, l] * ΔW[l]
                    y2 += tab.qdiff0.â[i, j] * c.B[j][k, l] * ΔW[l]
                end
            end
            b[D * (i - 1) + k] = -c.Y[stageindex(0), i][k] + (y1 + y2)
        end
    end

    # the terms corresponding to Y^l
    for i in 1:S
        for l in 1:M
            for k in 1:D
                local y1 = zero(ST)
                local y2 = zero(ST)
                for j in 1:S
                    y1 += tab.qdrift1.a[i, j] * c.V[j][k] * Δt
                    y2 += tab.qdrift1.â[i, j] * c.V[j][k] * Δt

                    # the noise terms, taken from the B1 tableau on the diagonal and B3 off it
                    for r in 1:M
                        if r == l
                            y1 += tab.qdiff1.a[i, j] * c.B[j][k, r] * ΔW[r]
                            y2 += tab.qdiff1.â[i, j] * c.B[j][k, r] * ΔW[r]
                        else
                            y1 += tab.qdiff3.a[i, j] * c.B[j][k, r] * ΔW[r]
                            y2 += tab.qdiff3.â[i, j] * c.B[j][k, r] * ΔW[r]
                        end
                    end
                end
                b[D * S + D * M * (i - 1) + D * (l - 1) + k] = -c.Y[stageindex(l), i][k] +
                                                               (y1 + y2)
            end
        end
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, x::AbstractVector{ST},
        sol, params,
        int::GeometricIntegrator{<:WIRK, <:AbstractProblemSDE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end

function GeometricIntegratorsBase.update!(sol, params, x::AbstractVector{DT},
        int::GeometricIntegrator{<:WIRK, <:AbstractProblemSDE}) where {DT}
    local c = cache(int, DT)
    local tab = tableau(int)

    components!(x, sol, params, int)

    # same final update as for SIRK
    update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift0.b, tab.qdiff0.b,
        timestep(int), c.ΔW, c.Δy)
    update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift0.b̂, tab.qdiff0.b̂,
        timestep(int), c.ΔW, c.Δy)
end

function GeometricIntegratorsBase.integrate_step!(sol, history, params,
        int::GeometricIntegrator{<:WIRK, <:AbstractProblemSDE})
    sample_noise!(int, sol)
    stage_initial_guess!(sol, history, params, int)

    local solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int),
        (sol, params, int))
    check_solver_status(solverstatus, int)

    update!(sol, params, nlsolution(int), int)
end
