
# Final updates of the stochastic Runge-Kutta methods.
#
# Each of these accumulates the whole increment of a state variable into a scratch vector and adds
# it in one go, rather than adding term by term.
#
# The addition goes through `GeometricBase.add!`, not `q .+= Δq`. `sol.q` is a `StateWithError`,
# which carries a running round-off error alongside the state, and `add!` is what performs the
# compensated summation against it — broadcasting falls through to the generic `AbstractArray`
# path, writes the state through `setindex!` and never touches the error field, so it drops the
# compensation entirely.
#
# The tableaus are applied twice by the callers, once with the weights `b` and once with the
# correction weights `b̂` that `RungeKutta.Tableau` carries — that pair is a higher-precision
# splitting of the same weights, and keeping the two additions separate is deliberate.

@doc raw"""
Final update of a stochastic Runge-Kutta method (SIRK, WIRK):
```math
q_{n+1} = q_{n} + \Delta t \sum_{i} b^{v}_{i} V_{i}
                + \sum_{r} \Delta W^r \sum_{i} b^{B}_{i} B_{i}^{\cdot r} .
```

* `V`: drift ``v`` evaluated at the internal stages
* `B`: diffusion matrix ``B`` evaluated at the internal stages
* `bdrift`, `bdiff`: weights of the drift and diffusion tableaus
* `ΔW`: increments of the driving process
* `Δq`, `Δy`: scratch of length ``D`` and ``M``
"""
function update_solution!(q, Δq, V, B, bdrift, bdiff, Δt::Number, ΔW, Δy)
    @assert length(bdrift) == length(bdiff) == length(V) == length(B)

    for i in eachindex(V, B)
        @assert length(Δq) == length(V[i]) == size(B[i], 1)
        @assert length(ΔW) == size(B[i], 2)
    end

    # contribution from the drift part
    for k in eachindex(Δq)
        local x = zero(eltype(Δq))
        for i in eachindex(bdrift, V)
            x += bdrift[i] * V[i][k]
        end
        Δq[k] = Δt * x
    end

    # contribution from the diffusion part
    for k in eachindex(Δq)
        Δy .= 0
        for i in eachindex(bdiff, B)
            for l in eachindex(Δy)
                Δy[l] += bdiff[i] * B[i][k, l]
            end
        end
        Δq[k] += dot(Δy, ΔW)
    end

    add!(q, Δq)
end

@doc raw"""
Final update of a stochastic explicit Runge-Kutta method (SERK), which carries a second diffusion
tableau against the iterated integrals ``\Delta Z``:
```math
q_{n+1} = q_{n} + \Delta t \sum_{i} b^{v}_{i} V_{i}
                + \sum_{r} \Delta W^r \sum_{i} b^{B}_{i} B_{i}^{\cdot r}
                + \frac{1}{\Delta t} \sum_{r} \Delta Z^r \sum_{i} b^{B_2}_{i} B_{i}^{\cdot r} .
```
"""
function update_solution!(q, Δq, V, B, bdrift, bdiff, bdiff2, Δt::Number, ΔW, ΔZ, Δy)
    @assert length(bdiff2) == length(B)

    for i in eachindex(B)
        @assert length(ΔZ) == size(B[i], 2)
    end

    # contributions from the drift and the ΔW part of the diffusion
    update_solution!(q, Δq, V, B, bdrift, bdiff, Δt, ΔW, Δy)

    # contribution from the ΔZ part of the diffusion
    for k in eachindex(Δq)
        Δy .= 0
        for i in eachindex(bdiff2, B)
            for l in eachindex(Δy)
                Δy[l] += bdiff2[i] * B[i][k, l]
            end
        end
        Δq[k] = dot(Δy, ΔZ) / Δt
    end

    add!(q, Δq)
end

@doc raw"""
Final update of a stochastic partitioned Runge-Kutta method (SIPRK):
```math
\begin{aligned}
q_{n+1} &= q_{n} + \Delta t \sum_{i} b^{v}_{i} V_{i}
                 + \sum_{r} \Delta W^r \sum_{i} b^{B}_{i} B_{i}^{\cdot r} , \\
p_{n+1} &= p_{n} + \Delta t \sum_{i} b^{f}_{i} F_{i}
                 + \sum_{r} \Delta W^r \sum_{i} b^{G}_{i} G_{i}^{\cdot r} .
\end{aligned}
```
"""
function update_solution!(q, p, Δq, Δp, V, F, B, G,
        bqdrift, bqdiff, bpdrift, bpdiff, Δt::Number, ΔW, Δy, Δz)
    @assert length(bqdrift) == length(bqdiff) == length(bpdrift) == length(bpdiff) ==
            length(V) == length(F) == length(B) == length(G)

    for i in eachindex(V, F, B, G)
        @assert length(Δq) == length(Δp) == length(V[i]) == length(F[i]) ==
                size(B[i], 1) == size(G[i], 1)
        @assert length(ΔW) == size(B[i], 2) == size(G[i], 2)
    end

    # contribution from the drift part
    for k in eachindex(Δq, Δp)
        local x = zero(eltype(Δq))
        local y = zero(eltype(Δp))
        for i in eachindex(bqdrift, bpdrift, V, F)
            x += bqdrift[i] * V[i][k]
            y += bpdrift[i] * F[i][k]
        end
        Δq[k] = Δt * x
        Δp[k] = Δt * y
    end

    # contribution from the diffusion part
    for k in eachindex(Δq, Δp)
        Δy .= 0
        Δz .= 0
        for i in eachindex(bqdiff, bpdiff, B, G)
            for l in eachindex(Δy, Δz)
                Δy[l] += bqdiff[i] * B[i][k, l]
                Δz[l] += bpdiff[i] * G[i][k, l]
            end
        end
        Δq[k] += dot(Δy, ΔW)
        Δp[k] += dot(Δz, ΔW)
    end

    add!(q, Δq)
    add!(p, Δp)
end

@doc raw"""
Final update of a stochastic split partitioned Runge-Kutta method (SISPRK). The ``p`` equation
carries two drift and two diffusion tableaus, applied to the Hamiltonian part ``F_1``, ``G_1`` and
the forcing part ``F_2``, ``G_2`` respectively:
```math
p_{n+1} = p_{n} + \Delta t \sum_{i} \left( b^{f_1}_{i} F^{1}_{i} + b^{f_2}_{i} F^{2}_{i} \right)
                + \sum_{r} \Delta W^r \sum_{i}
                  \left( b^{G_1}_{i} G^{1 \cdot r}_{i} + b^{G_2}_{i} G^{2 \cdot r}_{i} \right) .
```
"""
function update_solution!(q, p, Δq, Δp, V, F1, F2, B, G1, G2,
        bqdrift, bqdiff, bpdrift1, bpdrift2, bpdiff1, bpdiff2, Δt::Number, ΔW, Δy, Δz)
    @assert length(bqdrift) == length(bqdiff) == length(bpdrift1) == length(bpdrift2) ==
            length(bpdiff1) == length(bpdiff2) == length(V) == length(F1) == length(F2) ==
            length(B) == length(G1) == length(G2)

    for i in eachindex(V, F1, F2, B, G1, G2)
        @assert length(Δq) == length(Δp) == length(V[i]) == length(F1[i]) ==
                length(F2[i]) == size(B[i], 1) == size(G1[i], 1) == size(G2[i], 1)
        @assert length(ΔW) == size(B[i], 2) == size(G1[i], 2) == size(G2[i], 2)
    end

    # contribution from the drift part
    for k in eachindex(Δq, Δp)
        local x = zero(eltype(Δq))
        local y = zero(eltype(Δp))
        for i in eachindex(bqdrift, bpdrift1, bpdrift2, V, F1, F2)
            x += bqdrift[i] * V[i][k]
            y += bpdrift1[i] * F1[i][k] + bpdrift2[i] * F2[i][k]
        end
        Δq[k] = Δt * x
        Δp[k] = Δt * y
    end

    # contribution from the diffusion part
    for k in eachindex(Δq, Δp)
        Δy .= 0
        Δz .= 0
        for i in eachindex(bqdiff, bpdiff1, bpdiff2, B, G1, G2)
            for l in eachindex(Δy, Δz)
                Δy[l] += bqdiff[i] * B[i][k, l]
                Δz[l] += bpdiff1[i] * G1[i][k, l] + bpdiff2[i] * G2[i][k, l]
            end
        end
        Δq[k] += dot(Δy, ΔW)
        Δp[k] += dot(Δz, ΔW)
    end

    add!(q, Δq)
    add!(p, Δp)
end

@doc raw"""
Final update of a weak explicit Runge-Kutta method (WERK):
```math
q_{n+1} = q_{n} + \Delta t \sum_{i} \alpha_{i} V_{i}
                + \sum_{r} \hat I^r \sum_{i} \beta^{1}_{i} B_{1,i}^{\cdot r}
                + \sqrt{\Delta t} \sum_{i} \beta^{2}_{i} \sum_{r} B_{2,i}^{\cdot r} .
```

* `B1`: diffusion evaluated at the internal stages ``H^{(r)}_i``
* `B2`: diffusion evaluated at the internal stages ``\hat H^{(r)}_i``

This carries its own name rather than being another `update_solution!` method. It takes two
families of diffusion stages and a ``\sqrt{\Delta t}`` term that no strong update has, so it is a
different formula rather than an overload of the same one — and as an overload it was ambiguous
with the SERK update above, which has the same arity and differs only in which positional slot
holds the time step.
"""
function update_solution_weak!(q, Δq, V, B1, B2, α, β1, β2, Δt::Number, ΔW, Δy)
    @assert length(α) == length(β1) == length(β2) == length(V) == length(B1) == length(B2)

    for i in eachindex(V, B1, B2)
        @assert length(Δq) == length(V[i]) == size(B1[i], 1) == size(B2[i], 1)
        @assert length(ΔW) == size(B1[i], 2) == size(B2[i], 2)
    end

    # contribution from the drift part
    for k in eachindex(Δq)
        local x = zero(eltype(Δq))
        for i in eachindex(α, V)
            x += α[i] * V[i][k]
        end
        Δq[k] = Δt * x
    end

    # contribution from the diffusion term carrying the random variables Î
    for k in eachindex(Δq)
        Δy .= 0
        for i in eachindex(β1, B1)
            for l in eachindex(Δy)
                Δy[l] += β1[i] * B1[i][k, l]
            end
        end
        Δq[k] += dot(Δy, ΔW)
    end

    # contribution from the second diffusion term
    for k in eachindex(Δq)
        local x = zero(eltype(Δq))
        for i in eachindex(β2, B2)
            for l in axes(B2[i], 2)
                x += β2[i] * B2[i][k, l]
            end
        end
        Δq[k] += sqrt(Δt) * x
    end

    add!(q, Δq)
end
