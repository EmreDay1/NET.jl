"""
    calculate_escape_rate(barrier_height, diffusion)

Calculates the Kramers escape rate from a metastable state.

# Arguments
- barrier_height::Float64: Potential barrier height.
- diffusion::Float64: Diffusion coefficient.

# Returns
- rate::Float64: Escape rate.
"""
function calculate_escape_rate(barrier_height, diffusion)
    return exp(-barrier_height / diffusion)
end
