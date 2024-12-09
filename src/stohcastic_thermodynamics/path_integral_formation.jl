using Optim

"""
    calculate_action(path, drift, diffusion)

Calculates the action functional for a given path in a stochastic system.

# Arguments
- path::Vector{Float64}: The trajectory of the system.
- drift::Function: Drift function, `drift(x)`.
- diffusion::Float64: Diffusion coefficient.

# Returns
- action::Float64: The action value.
"""
function calculate_action(path, drift, diffusion)
    action = 0.0
    for i in 1:(length(path) - 1)
        dx = path[i+1] - path[i]
        drift_term = (drift(path[i])^2 / (2 * diffusion)) - drift(path[i]) * dx / diffusion
        action += (dx^2 / (2 * diffusion)) + drift_term
    end
    return action
end
