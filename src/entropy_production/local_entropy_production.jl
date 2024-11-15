using LinearAlgebra
using Logging
using Plots

"""
    compute_local_entropy_production(J::Array{Float64, 2}, X::Array{Float64, 2}) -> Array{Float64, 1}

Calculates the local entropy production at each grid point based on the fluxes (J) and corresponding thermodynamic forces (X).
- J: 2D array of fluxes at each point, where each column represents a flux vector at a grid point.
- X: 2D array of thermodynamic forces at each point, where each column represents a force vector at a grid point.
Returns: A 1D array of entropy production values for each grid point.
"""
function compute_local_entropy_prozduction(J::Array{Float64, 2}, X::Array{Float64, 2})
    validate_dimensions_entropy(J, X)
    entropy_production = sum(J .* X, dims=1)  
    log_info("Local entropy production calculation complete for $(size(J, 2)) grid points.")
    return vec(entropy_production)  
end

"""
    validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})

Ensures that the input dimensions of the fluxes (J) and forces (X) match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})
    if size(J) != size(X)
        throw(DimensionMismatch("J and X must have the same dimensions."))
    end
end

"""
    visualize_local_entropy_production_heatmap(σ::Array{Float64, 1}, grid_points::Vector, title_name::String)

Creates a heatmap to visualize local entropy production across grid points.
- σ: 1D array of entropy production values at each grid point.
- grid_points: Vector representing the grid points for the x-axis.
"""
function visualize_local_entropy_production_heatmap(σ::Array{Float64, 1}, grid_points::Vector, title_name::String)

    σ_matrix = reshape(σ, 1, length(σ))
    
    p = heatmap(
        grid_points,
        1:1,  
        σ_matrix,
        xlabel="Grid Points",
        ylabel="Entropy Production",
        title=title_name,
        color=:plasma,
        colorbar=true
    )
    display(p) 
end
