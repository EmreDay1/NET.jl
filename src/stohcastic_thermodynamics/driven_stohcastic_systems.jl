"""
    simulate_driven_system(drift, external_force, t_range)

Simulates a system with a time-dependent external driving force.

# Arguments
- drift::Function: Drift function, `drift(x)`.
- external_force::Function: Time-dependent force, `force(t)`.
- t_range::Vector{Float64}: Time range.

# Returns
- solution::Matrix{Float64}: Time series of positions.
"""
function simulate_driven_system(drift, external_force, t_range)
    positions = zeros(length(t_range))
    for t in 2:length(t_range)
        dt = t_range[t] - t_range[t-1]
        positions[t] = positions[t-1] + drift(positions[t-1]) * dt + external_force(t_range[t-1]) * dt
    end
    return positions
end
