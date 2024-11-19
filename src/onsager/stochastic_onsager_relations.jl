using LinearAlgebra
using Logging
using Random
using Plots

"""
    struct StochasticOnsagerMatrix

Represents a phenomenological matrix for stochastic Onsager relations.
- L: 3D array where L[:,:,k] is the matrix of coefficients at grid point or time step k.
"""
mutable struct StochasticOnsagerMatrix
    L::Array{Float64, 3}
end

"""
    compute_stochastic_fluxes(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2}) -> Array{Float64, 2}

Computes the fluxes (J) with stochastic contributions at each point based on the forces (F) using the phenomenological matrix (L).
- L: A `StochasticOnsagerMatrix` type representing the stochastic phenomenological coefficients.
- F: 2D array representing forces at each point.
- noise: 2D array representing noise contributions at each point.
Returns: A 2D array of computed fluxes with stochastic contributions.
"""
function compute_stochastic_fluxes(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2})
    validate_dimensions_stochastic(L, F, noise)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)

    for k in 1:num_points
        J[:, k] = L.L[:, :, k] * F[:, k] + noise[:, k]
    end

    log_info("Stochastic flux computation complete for $num_points points with external noise.")
    return J
end

"""
    validate_dimensions_stochastic(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2})

Ensures that the input dimensions of the matrix L, the force array F, and the noise array match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_stochastic(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2})
    num_vars, num_points = size(F)
    if size(L.L, 1) != num_vars || size(L.L, 3) != num_points
        throw(DimensionMismatch("L must match the number of variables and points in F."))
    end
    if size(noise) != (num_vars, num_points)
        throw(DimensionMismatch("Noise dimensions must match the dimensions of F (variables × points)."))
    end
end

"""
    visualize_stochastic_fluxes(J::Array{Float64, 2}, points::Vector, title_name::String)

Creates a heatmap to visualize stochastic fluxes over points (e.g., space or time).
- J: 2D array of fluxes at each point.
- points: Vector representing the x-axis (e.g., spatial or time points).
"""
function visualize_stochastic_fluxes(J::Array{Float64, 2}, points::Vector, title_name::String)
    p = heatmap(
        points,
        1:size(J, 1),
        J,
        xlabel="Grid Points",
        ylabel="Flux Variables",
        title=title_name,
        color=:inferno,
        colorbar=true
    )
    display(p)
end
