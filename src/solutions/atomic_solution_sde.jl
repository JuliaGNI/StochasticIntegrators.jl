"""
Atomic solution for an SDE.

### Parameters

* `DT`: data type
* `TT`: time step type
* `AT`: array type
* `IT`: internal variable types

### Fields

* `t`: time of current time step
* `t̄`: time of previous time step
* `q`: current solution of q
* `q̄`: previous solution of q
* `q̃`: compensated summation error of q
* `ΔW`: Wiener process driving the stochastic process q
* `ΔZ`: Wiener process driving the stochastic process q
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions)
"""
mutable struct SolutionStepSDE{DT<:Number,TT<:Real,AT<:AbstractArray{DT},IT<:NamedTuple} <: SolutionStep{DT,TT,AT}
    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    ΔW::AT
    ΔZ::AT

    K::Int

    internal::IT

    function SolutionStepSDE{DT,TT,AT,IT}(nd, nm, internal::IT) where {DT,TT,AT,IT}
        new(zero(TT), zero(TT),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nm)), AT(zeros(DT, nm)), 0,
            internal)
    end

    function SolutionStepSDE{DT,TT,AT,IT}(t::TT, q::AT, ΔW::AT, ΔZ::AT, internal::IT) where {DT,TT,AT,IT}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(ΔW), zero(ΔZ), 0,
            internal)
    end
end

SolutionStepSDE(::Type{DT}, ::Type{TT}, ::Type{AT}, nd, nm, internal::IT=NamedTuple()) where {DT,TT,AT,IT} = SolutionStepSDE{DT,TT,AT,IT}(nd, nm, internal)
SolutionStepSDE(t::TT, q::AT, ΔW::AT, ΔZ::AT, internal::IT=NamedTuple()) where {DT,TT,AT<:AbstractArray{DT},IT} = SolutionStepSDE{DT,TT,AT,IT}(t, q, ΔW, ΔZ, internal)

function set_solution!(asol::SolutionStepSDE, sol)
    t, q = sol
    asol.t = t
    asol.q .= q
end

function set_increments!(asol::SolutionStepSDE, incs)
    ΔW, ΔZ = incs
    asol.ΔW .= ΔW
    asol.ΔZ .= ΔZ
end

function get_solution(asol::SolutionStepSDE)
    (asol.t, asol.q)
end

function get_increments(asol::SolutionStepSDE)
    (asol.ΔW, asol.ΔZ)
end

function get_increments!(asol::SolutionStepSDE, ΔW, ΔZ)
    ΔW .= asol.ΔW
    ΔZ .= asol.ΔZ
end

function reset!(asol::SolutionStepSDE, Δt)
    asol.t̄ = asol.t
    asol.q̄ .= asol.q
    asol.t += Δt
end

function update!(asol::SolutionStepSDE{DT}, y::Vector{DT}) where {DT}
    for k in eachindex(y)
        update!(asol, y[k], k)
    end
end

function update!(asol::SolutionStepSDE{DT}, y::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
end
