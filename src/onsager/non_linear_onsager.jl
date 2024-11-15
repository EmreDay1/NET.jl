using LinearAlgebra
using Logging
using Plots

"""
    struct NonLinearOnsagerMatrix

Represents a phenomenological coefficient matrix for non-linear Onsager relations.
- L: 3D array for spatially varying phenomenological coefficients.
"""
mutable struct NonLinearOnsagerMatrix
    L::Array{Float64, 3}
end

"""
    compute_nonlinear_fluxes(L::NonLinearOnsagerMatrix, F::Array{Float64, 2}, F_nl::Array{Float64, 2}) -> Array{Float64, 2}

Computes the fluxes (J) from forces (F) with non-linear contributions using the phenomenological coefficient matrix (L).
- L: A `NonLinearOnsagerMatrix` type representing the phenomenological coefficients.
- F: 2D array of linear forces at each spatial point.
- F_nl: 2D array of non-linear corrections to the forces at each spatial point.
Returns: A 2D array of computed fluxes with non-linear contributions.
"""
function compute_nonlinear_fluxes(L::NonLinearOnsagerMatrix, F::Array{Float64, 2}, F_nl::Array{Float64, 2})
    validate_dimensions_nonlinear(L, F, F_nl)
    num_vars, num_points = size(F)
    J_linear = zeros(num_vars, num_points)
    J_nonlinear = zeros(num_vars, num_points)
    
    for k in 1:num_points
        J_linear[:, k] = L.L[:, :, k] * F[:, k]
        J_nonlinear[:, k] = L.L[:, :, k] * F_nl[:, k]
    end
    
    J_total = J_linear + J_nonlinear
    log_info("Non-linear flux computation complete for $num_points points.")
    return J_total
end

"""
    validate_dimensions_nonlinear(L::NonLinearOnsagerMatrix, F::Array{Float64, 2}, F_nl::Array{Float64, 2})

Ensures that the input dimensions of the matrix L, the linear force array F, and the non-linear force array F_nl match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_nonlinear(L::NonLinearOnsagerMatrix, F::Array{Float64, 2}, F_nl::Array{Float64, 2})
    num_vars, num_points = size(F)
    if size(F) != size(F_nl)
        throw(DimensionMismatch("F and F_nl must have the same dimensions."))
    end
    if size(L.L, 1) != num_vars || size(L.L, 3) != num_points
        throw(DimensionMismatch("L must match the number of variables and points in F and F_nl."))
    end
end

"""
    visualize_nonlinear_fluxes(J::Array{Float64, 2}, title_name::String)

Creates a heatmap visualization of the non-linear fluxes for easy interpretation.
- J: 2D array of fluxes at each spatial point.
"""
function visualize_nonlinear_fluxes(J::Array{Float64, 2}, title_name::String)
    heatmap(J, xlabel="Grid Points", ylabel="Flux Variables",
            title=title_name, color=:plasma, colorbar=true)
end

