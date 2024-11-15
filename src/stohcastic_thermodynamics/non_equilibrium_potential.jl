using LinearAlgebra

"""
    calculate_non_equilibrium_potential(grid, drift, diffusion)

Calculates the non-equilibrium potential landscape for a given drift and diffusion.

# Arguments
- grid::Vector{Float64}: Spatial grid for the potential.
- drift::Function: Drift function, `drift(x)`.
- diffusion::Float64: Diffusion coefficient.

# Returns
- potential::Vector{Float64}: Non-equilibrium potential landscape.
"""
function calculate_non_equilibrium_potential(grid, drift, diffusion)
    potential = zeros(length(grid))
    for i in 2:length(grid)
        dx = grid[i] - grid[i-1]
        potential[i] = potential[i-1] + drift(grid[i]) * dx / diffusion
    end
    return potential
end
