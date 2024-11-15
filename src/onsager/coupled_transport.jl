using LinearAlgebra
using Logging
using Plots

"""
    struct CoupledTransportMatrix

Represents a phenomenological matrix for coupled transport systems with multiple interacting variables.
- L: 3D array for spatially varying transport coefficients.
"""
mutable struct CoupledTransportMatrix
    L::Array{Float64, 3}
end

"""
    compute_coupled_fluxes(L::CoupledTransportMatrix, F::Array{Float64, 2}) -> Array{Float64, 2}

Computes the fluxes (J) based on interacting forces (F) using the coupled transport matrix (L).
- L: A `CoupledTransportMatrix` type representing the coupled transport coefficients.
- F: 2D array representing forces at each spatial point.
Returns: A 2D array of computed fluxes.
"""
function compute_coupled_fluxes(L::CoupledTransportMatrix, F::Array{Float64, 2})
    validate_dimensions_coupled(L, F)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)
    
    for k in 1:num_points
        J[:, k] = L.L[:, :, k] * F[:, k]
    end
    
    log_info("Coupled transport flux computation complete for $num_points points.")
    return J
end

"""
    validate_dimensions_coupled(L::CoupledTransportMatrix, F::Array{Float64, 2})

Ensures that the input dimensions of the matrix L and the force array F match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_coupled(L::CoupledTransportMatrix, F::Array{Float64, 2})
    num_vars, num_points = size(F)
    if size(L.L, 1) != num_vars || size(L.L, 3) != num_points
        throw(DimensionMismatch("L must match the number of variables and points in F."))
    end
end

"""
    visualize_coupled_fluxes(J::Array{Float64, 2}, title_name::String)

Creates a heatmap visualization of the coupled fluxes for easy interpretation.
- J: 2D array of fluxes at each spatial point.
"""
function visualize_coupled_fluxes(J::Array{Float64, 2}, title_name::String)
    heatmap(J, xlabel="Grid Points", ylabel="Flux Variables",
            title=title_name, color=:cividis, colorbar=true)
end

