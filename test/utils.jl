using LinearAlgebra: norm

@doc raw"""
    grid_length(timespan, timestep)

Number of increment columns a `GridProcess` must carry to drive a run over `timespan` with
`timestep`.

This is **not** simply `(t₁ - t₀) / Δt`. `GeometricEquations` validates a `GridProcess` against
`ntimesteps`, which rounds the quotient up using exact arithmetic rather than the rounded
floating-point division — and for a step that is not exactly representable in binary the two
disagree. The Kubo default is the case in point: `0.1 / 0.01` evaluates to exactly `10.0`, but
`0.1` is a shade above one tenth, so `div(0.1, 0.01, RoundUp)` is `11` and a process carrying ten
columns is rejected.

Matching the upstream computation here keeps these tests honest about what the API actually
demands. See the note in CHANGELOG.md: the run itself takes ten steps, so `ntime(problem)` and
`ntime(solution)` disagree by one for such a time step.
"""
grid_length(timespan, timestep) = Int(abs(div(timespan[end] - timespan[begin], timestep,
    RoundUp)))

# The Kubo oscillator conserves its energy exactly along every sample path, because its diffusion
# is proportional to its drift and the solution is therefore the deterministic one evaluated at
# the random time θ(t) = t + νW(t). Any drift in these quantities is the integrator's.

"Relative error of the energy `‖q‖²/2` between the first and last step of an SDE solution."
function rel_energy_err(sol::GeometricSolution)
    en_ref = 0.5 * norm(sol.q[begin])^2
    en_end = 0.5 * norm(sol.q[end])^2
    abs((en_end - en_ref) / en_ref)
end

"Relative error of the energy `(‖q‖² + ‖p‖²)/2` between the first and last step of a PSDE solution."
function rel_energy_err_pq(sol::GeometricSolution)
    en_ref = 0.5 * (norm(sol.q[begin])^2 + norm(sol.p[begin])^2)
    en_end = 0.5 * (norm(sol.q[end])^2 + norm(sol.p[end])^2)
    abs((en_end - en_ref) / en_ref)
end

"Largest relative energy error over the members of an ensemble solution."
function max_rel_energy_err(sol::EnsembleSolution)
    maximum(rel_energy_err, sol)
end

function max_rel_energy_err_pq(sol::EnsembleSolution)
    maximum(rel_energy_err_pq, sol)
end
