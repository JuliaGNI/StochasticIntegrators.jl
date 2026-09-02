@doc raw"""
    TableauWERK{T}

Tableau of a weak explicit Runge-Kutta method, in the layout of Rößler, *Second order Runge-Kutta
methods for Stratonovich stochastic differential equations*, BIT 47 (2007), Eq. (5.1).

Seven strictly lower triangular Butcher tableaus. In the notation of the paper,

| field | paper | field | paper |
|:--|:--|:--|:--|
| `qdrift0` | ``A^{(0)}`` | `qdiff0` | ``B^{(0)}`` |
| `qdrift1` | ``A^{(1)}`` | `qdiff1` | ``B^{(1)}`` |
| `qdrift2` | ``A^{(2)}`` | `qdiff2` | ``B^{(2)}`` |
|           |             | `qdiff3` | ``B^{(3)}`` |

with weights `qdrift0.b` ``= \alpha``, `qdiff0.b` ``= \beta_1``, `qdiff3.b` ``= \beta_2`` and
nodes `qdrift0.c`, `qdrift1.c`, `qdrift2.c` ``= c^{(0)}, c^{(1)}, c^{(2)}``. The `o` fields are
meaningless here and set to zero.
"""
struct TableauWERK{T} <: AbstractTableau{T}
    name::Symbol
    s::Int

    qdrift0::Tableau{T}
    qdrift1::Tableau{T}
    qdrift2::Tableau{T}

    qdiff0::Tableau{T}
    qdiff1::Tableau{T}
    qdiff2::Tableau{T}
    qdiff3::Tableau{T}

    function TableauWERK{T}(name, qdrift0, qdrift1, qdrift2,
            qdiff0, qdiff1, qdiff2, qdiff3) where {T}
        @assert qdrift0.s == qdrift1.s == qdrift2.s ==
                qdiff0.s == qdiff1.s == qdiff2.s == qdiff3.s
        @assert qdrift0.c[1] == qdrift1.c[1] == qdrift2.c[1] ==
                qdiff0.c[1] == qdiff1.c[1] == qdiff2.c[1] == qdiff3.c[1] == 0
        for tab in (qdrift0, qdrift1, qdrift2, qdiff0, qdiff1, qdiff2, qdiff3)
            @assert istrilstrict(tab.a)
            @assert !(tab.s == 1 && tab.a[1, 1] ≠ 0)
        end
        new(name, qdrift0.s, qdrift0, qdrift1, qdrift2, qdiff0, qdiff1, qdiff2, qdiff3)
    end
end

function TableauWERK(name, qdrift0::Tableau{T}, qdrift1::Tableau{T}, qdrift2::Tableau{T},
        qdiff0::Tableau{T}, qdiff1::Tableau{T}, qdiff2::Tableau{T},
        qdiff3::Tableau{T}) where {T}
    TableauWERK{T}(name, qdrift0, qdrift1, qdrift2, qdiff0, qdiff1, qdiff2, qdiff3)
end

function TableauWERK(name::Symbol, A0::AbstractMatrix{T}, A1::AbstractMatrix{T},
        A2::AbstractMatrix{T}, B0::AbstractMatrix{T}, B1::AbstractMatrix{T},
        B2::AbstractMatrix{T}, B3::AbstractMatrix{T}, α::AbstractVector{T},
        β1::AbstractVector{T}, β2::AbstractVector{T}, c0::AbstractVector{T},
        c1::AbstractVector{T}, c2::AbstractVector{T}) where {T}
    TableauWERK{T}(name, Tableau(name, 0, A0, α, c0),
        Tableau(name, 0, A1, α, c1), Tableau(name, 0, A2, α, c2),
        Tableau(name, 0, B0, β1, c0), Tableau(name, 0, B1, β1, c1),
        Tableau(name, 0, B2, β1, c2), Tableau(name, 0, B3, β2, c1))
end

@doc raw"""
    WERK(tableau; rng = Random.default_rng())

**Weak** explicit Runge-Kutta method for a `SDEProblem`.

A weak method approximates expectations ``E[\varphi(q(T))]`` rather than individual sample paths,
which is what one actually wants for computing statistics — mean energy, ergodic limits,
distributions. That freedom is what makes it cheap: the true Gaussian increments are replaced by
the discrete random variables of [`sample_noise!`](@ref), which match the moments the scheme was
derived against and cost a single uniform draw apiece.

Structurally it is more elaborate than a strong method. Beyond the drift stages ``H^{(0)}_i`` it
carries a separate family of stages ``H^{(l)}_i`` and ``\hat H^{(l)}_i`` for each noise dimension
``l``, and the diffusion matrix is evaluated at a *different* stage vector for each of its
columns. That is what buys second-order weak accuracy for multi-dimensional noise.

Constructed through [`RoesslerRS1`](@ref) or [`RoesslerRS2`](@ref).
"""
struct WERK{TT, RNG} <: SDEMethod
    tableau::TableauWERK{TT}
    rng::RNG
end

WERK(tableau::TableauWERK; rng = Random.default_rng()) = WERK(tableau, rng)

convergence(::WERK) = :weak
GeometricIntegratorsBase.isexplicit(::WERK) = true
GeometricIntegratorsBase.isimplicit(::WERK) = false

@doc raw"""
Cache of a weak explicit Runge-Kutta method.

### Fields

* `ΔW`, `ΔZ`: the discrete random variables ``\hat I`` and ``\tilde I`` of the current step
* `Δq`, `Δy`: scratch for the final update
* `Q0`: the internal stage ``H^{(0)}_i`` of the current stage index
* `Q1`, `Q2`: the internal stages ``H^{(l)}_i`` and ``\hat H^{(l)}_i``, one per noise dimension
* `V`: drift evaluated at the ``H^{(0)}_i``
* `B1`, `B2`: diffusion, column `l` evaluated at ``H^{(l)}_i`` and ``\hat H^{(l)}_i`` respectively
* `tB`: scratch holding a full diffusion matrix from which one column is taken
"""
struct WERKCache{DT, D, M, S} <: SDEIntegratorCache{DT}
    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    Δq::Vector{DT}
    Δy::Vector{DT}

    Q0::Vector{DT}
    Q1::Vector{Vector{DT}}
    Q2::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    B1::Vector{Matrix{DT}}
    B2::Vector{Matrix{DT}}

    tB::Matrix{DT}

    function WERKCache{DT, D, M, S}() where {DT, D, M, S}
        new(zeros(DT, M), zeros(DT, M), zeros(DT, D), zeros(DT, M),
            zeros(DT, D),
            create_internal_stage_vector(DT, D, M),
            create_internal_stage_vector(DT, D, M),
            create_internal_stage_vector(DT, D, S),
            create_internal_stage_matrix(DT, D, M, S),
            create_internal_stage_matrix(DT, D, M, S),
            zeros(DT, D, M))
    end
end

function GeometricIntegratorsBase.Cache{ST}(problem::AbstractProblemSDE,
        method::WERK; kwargs...) where {ST}
    WERKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}(;
        kwargs...)
end

@inline function GeometricIntegratorsBase.CacheType(ST, problem::AbstractProblemSDE,
        method::WERK)
    WERKCache{ST, length(vec(initial_conditions(problem).q)),
        noisedims(problem), nstages(method)}
end

@doc raw"""
Evaluate column `l` of the diffusion matrix at the stage vector `q` and store it in column `l` of
`B`.

The scheme needs each column of ``B`` at a *different* argument, but the equation interface only
offers the whole matrix at one argument, so the full matrix is formed in scratch and one column
kept. That costs ``m`` times more diffusion evaluations than a column-wise interface would; it is
the price of not widening `GeometricEquations`, and it is why the weak methods are noticeably more
expensive per stage than the strong ones for multi-dimensional noise.
"""
function diffusion_column!(B, tB, l, equ, t, q, params)
    equ.B(tB, t, q, params)
    for k in axes(B, 1)
        B[k, l] = tB[k, l]
    end
    B
end

function GeometricIntegratorsBase.integrate_step!(sol, history, params,
        int::GeometricIntegrator{<:WERK, <:AbstractProblemSDE})
    # draw the increments of the driving process for this step
    sample_noise!(int, sol)

    local c = cache(int)
    local tab = tableau(int)
    local Δt = timestep(int)
    local t̄ = sol.t - Δt
    local equ = equations(int)

    # the first internal stages all equal the solution at the previous step
    equ.v(c.V[1], t̄, sol.q, params)
    equ.B(c.B1[1], t̄, sol.q, params)
    c.B2[1] .= c.B1[1]

    for i in 2:tab.s
        # the internal stage H^(0)_i
        for k in eachindex(c.Q0)
            local ydrift = zero(eltype(c.Q0))
            for j in 1:(i - 1)
                ydrift += tab.qdrift0.a[i, j] * c.V[j][k]
            end

            c.Δy .= 0
            for j in 1:(i - 1)
                for l in eachindex(c.Δy)
                    c.Δy[l] += tab.qdiff0.a[i, j] * c.B1[j][k, l]
                end
            end

            c.Q0[k] = sol.q[k] + Δt * ydrift + dot(c.Δy, c.ΔW)
        end

        # the internal stages H^(l)_i, one per noise dimension
        for k in eachindex(c.Q1[1])
            local ydrift = zero(eltype(c.Q0))
            for j in 1:(i - 1)
                ydrift += tab.qdrift1.a[i, j] * c.V[j][k]
            end

            for r in eachindex(c.Q1)
                c.Δy .= 0
                for j in 1:(i - 1)
                    for l in eachindex(c.Δy)
                        if l == r
                            c.Δy[l] += tab.qdiff1.a[i, j] * c.B1[j][k, l]
                        else
                            c.Δy[l] += tab.qdiff3.a[i, j] * c.B1[j][k, l]
                        end
                    end
                end

                c.Q1[r][k] = sol.q[k] + Δt * ydrift + dot(c.Δy, c.ΔW)
            end
        end

        # the internal stages Ĥ^(l)_i, one per noise dimension
        for k in eachindex(c.Q2[1])
            local ydrift = zero(eltype(c.Q0))
            for j in 1:(i - 1)
                ydrift += tab.qdrift2.a[i, j] * c.V[j][k]
            end

            for r in eachindex(c.Q2)
                local ydiff2 = zero(eltype(c.Q0))
                for j in 1:(i - 1)
                    for l in eachindex(c.ΔW, c.ΔZ)
                        # the terms carrying the random variables Î_{r,l}
                        if l < r
                            ydiff2 += tab.qdiff2.a[i, j] * c.B1[j][k, l] *
                                      c.ΔW[r] * c.ΔZ[l]
                        elseif l > r
                            ydiff2 += -tab.qdiff2.a[i, j] * c.B1[j][k, l] *
                                      c.ΔW[l] * c.ΔZ[r]
                        end
                    end
                end

                c.Q2[r][k] = sol.q[k] + Δt * ydrift + ydiff2 / sqrt(Δt)
            end
        end

        # the new values of V
        equ.v(c.V[i], t̄ + Δt * tab.qdrift0.c[i], c.Q0, params)

        # the new values of B1 and B2, each column at its own internal stage
        local t1 = t̄ + Δt * tab.qdrift1.c[i]
        local t2 = t̄ + Δt * tab.qdrift2.c[i]

        for l in eachindex(c.Q1)
            diffusion_column!(c.B1[i], c.tB, l, equ, t1, c.Q1[l], params)
        end

        for l in eachindex(c.Q2)
            diffusion_column!(c.B2[i], c.tB, l, equ, t2, c.Q2[l], params)
        end
    end

    # compute the final update
    update_solution_weak!(sol.q, c.Δq, c.V, c.B1, c.B2, tab.qdrift0.b, tab.qdiff0.b,
        tab.qdiff3.b, Δt, c.ΔW, c.Δy)
end
