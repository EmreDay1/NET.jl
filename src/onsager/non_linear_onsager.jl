using LinearAlgebra
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
    struct NonLinearFluxSystem

Represents a system of non-linear Onsager relations with multiple non-linear corrections.
- matrix: A `NonLinearOnsagerMatrix` type for the phenomenological coefficients.
- linear_forces: 2D array for linear forces across spatial points.
- nonlinear_forces: Vector of 2D arrays for multiple non-linear force corrections.
"""
mutable struct NonLinearFluxSystem
    matrix::NonLinearOnsagerMatrix
    linear_forces::Array{Float64, 2}
    nonlinear_forces::Vector{Array{Float64, 2}}
end

"""
    compute_total_fluxes(system::NonLinearFluxSystem) -> Array{Float64, 2}

Computes the total fluxes for a system with multiple non-linear contributions.
- system: A `NonLinearFluxSystem` object.
Returns: A 2D array of computed fluxes with both linear and non-linear contributions.
"""
function compute_total_fluxes(system::NonLinearFluxSystem)

    L = system.matrix.L
    F = system.linear_forces
    F_nl_list = system.nonlinear_forces

  
    validate_dimensions_multinonlinear(system)

    num_vars, num_points = size(F)
    J_linear = zeros(num_vars, num_points)
    J_total = zeros(num_vars, num_points)


    for k in 1:num_points
        J_linear[:, k] = L[:, :, k] * F[:, k]
    end
    J_total .= J_linear 


    for F_nl in F_nl_list
        for k in 1:num_points
            J_total[:, k] .+= L[:, :, k] * F_nl[:, k]
        end
    end


    return J_total
end

"""
    validate_dimensions_multinonlinear(system::NonLinearFluxSystem)

Ensures that the input dimensions of the matrix L, linear forces F, and all non-linear force arrays match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_multinonlinear(system::NonLinearFluxSystem)
    L = system.matrix.L
    F = system.linear_forces
    F_nl_list = system.nonlinear_forces

    num_vars, num_points = size(F)

   
    if size(L, 1) != num_vars || size(L, 3) != num_points
        throw(DimensionMismatch("L must match the number of variables and points in F."))
    end

   
    for F_nl in F_nl_list
        if size(F) != size(F_nl)
            throw(DimensionMismatch("All forces (F and F_nl) must have the same dimensions."))
        end
    end
end

"""
    visualize_total_fluxes(J::Array{Float64, 2}, title_name::String)

Creates a heatmap visualization of the total fluxes (including non-linear contributions).
- J: 2D array of fluxes at each spatial point.
"""
function visualize_total_fluxes(J::Array{Float64, 2}, title_name::String)
    heatmap(J, xlabel="Grid Points", ylabel="Flux Variables",
            title=title_name, color=:plasma, colorbar=true)
end
