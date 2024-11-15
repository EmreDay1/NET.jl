using LinearAlgebra
using Logging
using Plots

"""
    struct SpatialOnsagerMatrix

Represents a spatially dependent phenomenological coefficient matrix for Onsager relations.
- L: 3D array where L[:,:,k] is the matrix of coefficients at grid point k.
"""
mutable struct SpatialOnsagerMatrix
    L::Array{Float64, 3}
end

"""
    compute_spatial_fluxes(L::SpatialOnsagerMatrix, F::Array{Float64, 2}) -> Array{Float64, 2}

Computes the fluxes (J) at each spatial point based on the forces (F) using the spatially varying Onsager matrix (L).
- L: A `SpatialOnsagerMatrix` type representing spatially dependent phenomenological coefficients.
- F: 2D array representing forces at each spatial point, where each column is a force vector at a grid point.
Returns: A 2D array of computed fluxes at each spatial point.
"""
function compute_spatial_fluxes(L::SpatialOnsagerMatrix, F::Array{Float64, 2})
    validate_dimensions_spatial(L, F)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)
    
    for k in 1:num_points
        J[:, k] = L.L[:, :, k] * F[:, k]
    end
    
    log_info("Spatial flux computation complete for $num_points points.")
    return J
end

"""
    validate_dimensions_spatial(L::SpatialOnsagerMatrix, F::Array{Float64, 2})

Ensures that the input dimensions of the matrix L and the force array F match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_spatial(L::SpatialOnsagerMatrix, F::Array{Float64, 2})
    num_vars, num_points = size(F)
    if size(L.L, 1) != num_vars || size(L.L, 3) != num_points
        throw(DimensionMismatch("L must match the number of variables and points in F."))
    end
end

"""
    visualize_spatial_fluxes(J::Array{Float64, 2}, title_name::String)

Creates a heatmap visualization of the spatial fluxes for easy interpretation.
- J: 2D array of fluxes at each spatial point.
"""
function visualize_spatial_fluxes(J::Array{Float64, 2}, title_name::String)
    heatmap(J, xlabel="Grid Points", ylabel="Flux Variables",
            title=title_name, color=:inferno, colorbar=true)
end

