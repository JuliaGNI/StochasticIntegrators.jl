using LinearAlgebra: norm

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
