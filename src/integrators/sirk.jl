@doc raw"""
    TableauSIRK{T}

Tableau of a stochastic implicit Runge-Kutta method: `qdrift` for the drift, `qdiff` for the
diffusion.

As for [`TableauSERK`](@ref) no overall order is stored, because the order of a stochastic
Runge-Kutta method depends on the noise as well as on the coefficients. The `o` fields are the
classical orders of the underlying deterministic schemes.
"""
struct TableauSIRK{T} <: AbstractTableau{T}
    name::Symbol
    s::Int
    qdrift::Tableau{T}
    qdiff::Tableau{T}

    function TableauSIRK{T}(name, s, qdrift, qdiff) where {T}
        @assert s == qdrift.s == qdiff.s
        new(name, s, qdrift, qdiff)
    end
end

function TableauSIRK(name::Symbol, qdrift::Tableau{T}, qdiff::Tableau{T}) where {T}
    TableauSIRK{T}(name, qdrift.s, qdrift, qdiff)
end

function TableauSIRK(name::Symbol, order_drift::Int, a_drift::AbstractMatrix{T},
        b_drift::AbstractVector{T}, c_drift::AbstractVector{T}, order_diff::Int,
        a_diff::AbstractMatrix{T}, b_diff::AbstractVector{T},
        c_diff::AbstractVector{T}) where {T}
    TableauSIRK{T}(name, length(c_drift),
        Tableau(name, order_drift, a_drift, b_drift, c_drift),
        Tableau(name, order_diff, a_diff, b_diff, c_diff))
end

@doc raw"""
    SIRK(tableau; K = 0, rng = Random.default_rng())

Stochastic **implicit** Runge-Kutta method for a `SDEProblem`, of strong (mean-square)
convergence.

The internal stages are coupled, so each step solves
```math
Y_i = \Delta t \sum_j a^{v}_{ij} \, v(Q_j) + \sum_r \Delta W^r \sum_j a^{B}_{ij} \, B^{\cdot r}(Q_j),
\qquad Q_i = q_n + Y_i
```
for the stage increments ``Y_i`` by Newton's method, in ``D \times S`` unknowns. Note that the
unknowns are the increments themselves rather than the stage derivatives: unlike the deterministic
case the drift and the diffusion enter with different coefficients, so there is no single
derivative to solve for.

`K` sets the Milstein-Tretyakov truncation of the Wiener increments — see [`truncation`](@ref).
It defaults to `0`, no truncation.

Constructed through [`StochasticGLRK`](@ref) or [`StochasticDIRK`](@ref).
"""
struct SIRK{TT, RNG} <: SDEMethod
    tableau::TableauSIRK{TT}
    K::Int
    rng::RNG
end

SIRK(tableau::TableauSIRK; K = 0, rng = Random.default_rng()) = SIRK(tableau, K, rng)

convergence(::SIRK) = :strong
truncation(method::SIRK) = method.K
GeometricIntegratorsBase.isexplicit(::SIRK) = false
GeometricIntegratorsBase.isimplicit(::SIRK) = true

@doc raw"""
Cache of a stochastic implicit Runge-Kutta method.

### Fields

* `x`: nonlinear solver solution vector, holding the stage increments ``Y``
* `ΔW`, `ΔZ`: increments of the driving process for the current step
* `Δq`, `Δy`: scratch for the final update
* `Q`, `V`, `B`, `Y`: internal stages, drift, diffusion and stage increments
* `ΔQ`, `Y1`, `Y2`, `V1`, `V2`, `B1`, `B2`, `Δw`: scratch of the initial-guess predictor
"""
struct SIRKCache{DT, D, M, S} <: SDEIntegratorCache{DT}
    x::Vector{DT}

    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    Δq::Vector{DT}
    Δy::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    B::Vector{Matrix{DT}}
    Y::Vector{Vector{DT}}

    ΔQ::Vector{DT}
    Y1::Vector{DT}
    Y2::Vector{DT}
    V1::Vector{DT}
    V2::Vector{DT}
    B1::Matrix{DT}
    B2::Matrix{DT}
    Δw::Vector{DT}

    function SIRKCache{DT, D, M, S}() where {DT, D, M, S}
        new(zeros(DT, D * S),
            zeros(DT, M), zeros(DT, M), zeros(DT, D), zeros(DT, M),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_matrix(DT, D, M, S),
            create_internal_stage_vector(DT, D, S),
            zeros(DT, D), zeros(DT, D), zeros(DT, D), zeros(DT, D), zeros(DT, D),
            zeros(DT, D, M), zeros(DT, D, M), zeros(DT, M))
    end
end

GeometricIntegratorsBase.nlsolution(cache::SIRKCache) = cache.x

function GeometricIntegratorsBase.Cache{ST}(problem::AbstractProblemSDE,
        method::SIRK; kwargs...) where {ST}
    SIRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}(;
        kwargs...)
end

@inline function GeometricIntegratorsBase.CacheType(ST, problem::AbstractProblemSDE,
        method::SIRK)
    SIRKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}
end

function GeometricIntegratorsBase.solversize(method::SIRK, problem::AbstractProblemSDE)
    length(vec(initial_conditions(problem).q)) * nstages(method)
end

GeometricIntegratorsBase.default_solver(::SIRK) = Newton()

@doc raw"""
Initial guess for the stage increments of a [`SIRK`](@ref) method.

Rather than extrapolating from the solution history, this runs the two-stage explicit R2 method of
Burrage & Burrage from ``(t_n, q_n)`` out to each collocation node ``c_i``, and takes the resulting
increment as the guess for ``Y_i``. The predictor is cheap — two evaluations of the drift and
diffusion per stage — and for time steps and noise intensities that are not too large it lands
close enough that Newton converges in very few iterations.

If the noise is strong enough that this predictor is poor, the alternatives are a smaller time
step, a higher-order explicit predictor, or simply `nlsolution(int) .= 0`.

This is deliberately *not* the framework's `initial_guess!` hook. That hook runs before
`integrate_step!`, and the predictor needs the increments of the step it is predicting, which are
not drawn until `integrate_step!` begins. It is called from there instead, right after
[`sample_noise!`](@ref).
"""
function stage_initial_guess!(sol, history, params,
        int::GeometricIntegrator{<:SIRK, <:AbstractProblemSDE})
    local c = cache(int)
    local tab = tableau(int)
    local Δt = timestep(int)
    local t̄ = sol.t - Δt
    local D = length(c.Δq)

    # evaluate the drift and the diffusion at (t̄, q̄) — the same for every stage
    equations(int).v(c.V1, t̄, sol.q, params)
    equations(int).B(c.B1, t̄, sol.q, params)

    for i in eachstage(method(int))
        local Δt_local = tab.qdrift.c[i] * Δt
        c.Δw .= tab.qdrift.c[i] .* c.ΔW

        mul!(c.ΔQ, c.B1, c.Δw)
        @. c.Q[i] = sol.q + 2 // 3 * Δt_local * c.V1 + 2 // 3 * c.ΔQ

        local t2 = t̄ + 2 // 3 * Δt_local

        equations(int).v(c.V2, t2, c.Q[i], params)
        equations(int).B(c.B2, t2, c.Q[i], params)

        mul!(c.Y1, c.B1, c.Δw)
        mul!(c.Y2, c.B2, c.Δw)

        for j in 1:D
            c.x[(i - 1) * D + j] = Δt_local * (1 // 4 * c.V1[j] + 3 // 4 * c.V2[j]) +
                                   1 // 4 * c.Y1[j] + 3 // 4 * c.Y2[j]
        end
    end
end

@doc raw"""
Unpack the solver vector `x = (Y_1^1, …, Y_1^D, Y_2^1, …)` into the stage increments `Y`, form the
internal stages `Q_i = q_n + Y_i`, and evaluate the drift and diffusion there.
"""
function GeometricIntegratorsBase.components!(x::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:SIRK, <:AbstractProblemSDE}) where {ST}
    local c = cache(int, ST)
    local tab = tableau(int)
    local Δt = timestep(int)
    local t̄ = sol.t - Δt
    local D = length(c.Δq)

    # copy x to Y and compute Q
    for i in eachindex(c.Q, c.Y)
        for k in eachindex(c.Q[i])
            c.Y[i][k] = x[D * (i - 1) + k]
            c.Q[i][k] = sol.q[k] + c.Y[i][k]
        end
    end

    # compute V = v(Q) and B = B(Q)
    for i in eachindex(c.Q, c.V, c.B)
        local tᵢ = t̄ + Δt * tab.qdrift.c[i]
        equations(int).v(c.V[i], tᵢ, c.Q[i], params)
        equations(int).B(c.B[i], tᵢ, c.Q[i], params)
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:SIRK, <:AbstractProblemSDE}) where {ST}
    local c = cache(int, ST)
    local tab = tableau(int)
    local Δt = timestep(int)
    local D = length(c.Δq)
    local S = nstages(method(int))

    # the increments are real numbers and shared by every arithmetic path in this step
    local ΔW = cache(int).ΔW

    # compute b = - (Y - A V - Â B ΔW)
    for i in 1:S
        for k in 1:D
            local y1 = zero(ST)
            local y2 = zero(ST)
            local y3 = zero(ST)
            local y4 = zero(ST)
            for j in 1:S
                y1 += tab.qdrift.a[i, j] * c.V[j][k] * Δt
                y2 += tab.qdrift.â[i, j] * c.V[j][k] * Δt
                for l in 1:length(ΔW)
                    y3 += tab.qdiff.a[i, j] * c.B[j][k, l] * ΔW[l]
                    y4 += tab.qdiff.â[i, j] * c.B[j][k, l] * ΔW[l]
                end
            end
            b[D * (i - 1) + k] = -c.Y[i][k] + (y1 + y2) + (y3 + y4)
        end
    end
end

function GeometricIntegratorsBase.residual!(b::AbstractVector{ST}, x::AbstractVector{ST},
        sol, params,
        int::GeometricIntegrator{<:SIRK, <:AbstractProblemSDE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end

function GeometricIntegratorsBase.update!(sol, params, x::AbstractVector{DT},
        int::GeometricIntegrator{<:SIRK, <:AbstractProblemSDE}) where {DT}
    local c = cache(int, DT)
    local tab = tableau(int)

    components!(x, sol, params, int)

    update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift.b, tab.qdiff.b,
        timestep(int), c.ΔW, c.Δy)
    update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift.b̂, tab.qdiff.b̂,
        timestep(int), c.ΔW, c.Δy)
end

function GeometricIntegratorsBase.integrate_step!(sol, history, params,
        int::GeometricIntegrator{<:SIRK, <:AbstractProblemSDE})
    # draw the increments of the driving process for this step
    sample_noise!(int, sol)

    # compute the initial guess, which needs the increments already in place
    stage_initial_guess!(sol, history, params, int)

    # call the nonlinear solver
    local solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int),
        (sol, params, int))
    check_solver_status(solverstatus, int)

    # compute the final update
    update!(sol, params, nlsolution(int), int)
end
