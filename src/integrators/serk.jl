@doc raw"""
    TableauSERK{T}

Tableau of a stochastic explicit Runge-Kutta method.

Three Butcher tableaus, all strictly lower triangular: `qdrift` for the drift, `qdiff` for the
diffusion terms carrying the Wiener increments ``\Delta W``, and `qdiff2` for the terms carrying
the iterated integrals ``\Delta Z``. A method that uses only ``\Delta W`` leaves `qdiff2` empty,
which is recorded by giving it the name `:NULL`; [`hasdiffusion2`](@ref) is the predicate.

No overall order is stored. Unlike in the deterministic case the order of a stochastic
Runge-Kutta method is not a property of the coefficients alone — it depends on the dimension of
the driving Wiener process and on whether the diffusion matrix satisfies the commutativity
conditions. The `o` fields of the three tableaus are the *classical* orders of the underlying
deterministic schemes and nothing more.
"""
struct TableauSERK{T} <: AbstractTableau{T}
    name::Symbol
    s::Int

    qdrift::Tableau{T}
    qdiff::Tableau{T}
    qdiff2::Tableau{T}

    function TableauSERK{T}(name, qdrift, qdiff, qdiff2) where {T}
        @assert qdrift.s == qdiff.s == qdiff2.s
        @assert qdrift.c[1] == qdiff.c[1] == qdiff2.c[1] == 0
        @assert istrilstrict(qdrift.a)
        @assert istrilstrict(qdiff.a)
        @assert istrilstrict(qdiff2.a)
        @assert !(qdrift.s == 1 && qdrift.a[1, 1] ≠ 0)
        @assert !(qdiff.s == 1 && qdiff.a[1, 1] ≠ 0)
        @assert !(qdiff2.s == 1 && qdiff2.a[1, 1] ≠ 0)
        new(name, qdrift.s, qdrift, qdiff, qdiff2)
    end
end

function TableauSERK(name, qdrift::Tableau{T}, qdiff::Tableau{T}, qdiff2::Tableau{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff, qdiff2)
end

function TableauSERK(name, qdrift::Tableau{T}, qdiff::Tableau{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff,
        Tableau{T}(:NULL, 0, zero(qdrift.a), zero(qdrift.b), zero(qdrift.c)))
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::AbstractMatrix{T},
        b_drift::AbstractVector{T}, c_drift::AbstractVector{T},
        order_diff::Int, a_diff::AbstractMatrix{T}, b_diff::AbstractVector{T},
        c_diff::AbstractVector{T}, order_diff2::Int, a_diff2::AbstractMatrix{T},
        b_diff2::AbstractVector{T}, c_diff2::AbstractVector{T}) where {T}
    TableauSERK{T}(name, Tableau(name, order_drift, a_drift, b_drift, c_drift),
        Tableau(name, order_diff, a_diff, b_diff, c_diff),
        Tableau(name, order_diff2, a_diff2, b_diff2, c_diff2))
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::AbstractMatrix{T},
        b_drift::AbstractVector{T}, c_drift::AbstractVector{T},
        order_diff::Int, a_diff::AbstractMatrix{T}, b_diff::AbstractVector{T},
        c_diff::AbstractVector{T}) where {T}
    TableauSERK{T}(name, Tableau(name, order_drift, a_drift, b_drift, c_drift),
        Tableau(name, order_diff, a_diff, b_diff, c_diff),
        Tableau{T}(:NULL, 0, zero(a_drift), zero(b_drift), zero(c_drift)))
end

"Whether a tableau carries a second diffusion tableau, applied to the iterated integrals `ΔZ`."
hasdiffusion2(tab::TableauSERK) = tab.qdiff2.name !== :NULL

@doc raw"""
    SERK(tableau; rng = Random.default_rng())

Stochastic **explicit** Runge-Kutta method for a `SDEProblem`, of strong (mean-square)
convergence.

Each internal stage is built from the stages before it, so a step costs one evaluation of the
drift and the diffusion per stage and involves no solve:
```math
Q_i = q_n + \Delta t \sum_{j<i} a^{v}_{ij} V_j
          + \sum_r \Delta W^r \sum_{j<i} a^{B}_{ij} B_j^{\cdot r}
          + \frac{1}{\Delta t} \sum_r \Delta Z^r \sum_{j<i} a^{B_2}_{ij} B_j^{\cdot r} .
```

Being explicit it is cheap per step but has a bounded stability region, which for a stiff drift or
a large noise intensity forces a much smaller time step than the implicit methods.

Constructed through one of the named schemes: [`StochasticEuler`](@ref),
[`StochasticHeun`](@ref), [`Platen`](@ref), [`BurrageR2`](@ref), [`BurrageCL`](@ref),
[`BurrageE1`](@ref), [`BurrageG5`](@ref).
"""
struct SERK{TT, RNG} <: SDEMethod
    tableau::TableauSERK{TT}
    rng::RNG
end

SERK(tableau::TableauSERK; rng = Random.default_rng()) = SERK(tableau, rng)

convergence(::SERK) = :strong
GeometricIntegratorsBase.isexplicit(::SERK) = true
GeometricIntegratorsBase.isimplicit(::SERK) = false

@doc raw"""
Cache of a stochastic explicit Runge-Kutta method.

### Fields

* `ΔW`, `ΔZ`: increments of the driving process for the current step
* `Δq`, `Δy`: scratch for the final update
* `Q`: internal stages, `Q[i][k]`
* `V`: drift evaluated at the internal stages
* `B`: diffusion matrix evaluated at the internal stages
"""
struct SERKCache{DT, D, M, S} <: SDEIntegratorCache{DT}
    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    Δq::Vector{DT}
    Δy::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    B::Vector{Matrix{DT}}

    function SERKCache{DT, D, M, S}() where {DT, D, M, S}
        new(zeros(DT, M), zeros(DT, M), zeros(DT, D), zeros(DT, M),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_matrix(DT, D, M, S))
    end
end

function GeometricIntegratorsBase.Cache{ST}(problem::AbstractProblemSDE,
        method::SERK; kwargs...) where {ST}
    SERKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}(;
        kwargs...)
end

@inline function GeometricIntegratorsBase.CacheType(ST, problem::AbstractProblemSDE,
        method::SERK)
    SERKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}
end

function GeometricIntegratorsBase.integrate_step!(sol, history, params,
        int::GeometricIntegrator{<:SERK, <:AbstractProblemSDE})
    # draw the increments of the driving process for this step
    sample_noise!(int, sol)

    local c = cache(int)
    local tab = tableau(int)
    local Δt = timestep(int)

    # on entry sol.t holds tₙ₊₁ while sol.q still holds qₙ
    local t̄ = sol.t - Δt

    for i in eachstage(method(int))
        for k in eachindex(c.Q[i])
            # contribution from the drift part
            local ydrift = zero(eltype(c.Δq))
            for j in 1:(i - 1)
                ydrift += tab.qdrift.a[i, j] * c.V[j][k]
            end

            # ΔW contribution from the diffusion part
            c.Δy .= 0
            for j in 1:(i - 1)
                for l in eachindex(c.Δy)
                    c.Δy[l] += tab.qdiff.a[i, j] * c.B[j][k, l]
                end
            end

            c.Q[i][k] = sol.q[k] + Δt * ydrift + dot(c.Δy, c.ΔW)

            # ΔZ contribution from the diffusion part
            if hasdiffusion2(tab)
                c.Δy .= 0
                for j in 1:(i - 1)
                    for l in eachindex(c.Δy)
                        c.Δy[l] += tab.qdiff2.a[i, j] * c.B[j][k, l]
                    end
                end

                c.Q[i][k] += dot(c.Δy, c.ΔZ) / Δt
            end
        end

        local tᵢ = t̄ + Δt * tab.qdrift.c[i]
        equations(int).v(c.V[i], tᵢ, c.Q[i], params)
        equations(int).B(c.B[i], tᵢ, c.Q[i], params)
    end

    # compute the final update
    if hasdiffusion2(tab)
        update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift.b, tab.qdiff.b, tab.qdiff2.b,
            Δt, c.ΔW, c.ΔZ, c.Δy)
        update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift.b̂, tab.qdiff.b̂, tab.qdiff2.b̂,
            Δt, c.ΔW, c.ΔZ, c.Δy)
    else
        update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift.b, tab.qdiff.b, Δt, c.ΔW, c.Δy)
        update_solution!(sol.q, c.Δq, c.V, c.B, tab.qdrift.b̂, tab.qdiff.b̂, Δt, c.ΔW, c.Δy)
    end
end
