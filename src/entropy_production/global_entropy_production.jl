using LinearAlgebra

"""
    compute_global_entropy_production(J::Array{Float64, 2}, X::Array{Float64, 2}) -> Float64

Calculates the global entropy production for the entire system by summing the local entropy productions across all grid points.
- J: 2D array of fluxes where each column represents a flux vector at a grid point.
- X: 2D array of thermodynamic forces where each column represents a force vector at a grid point.
Returns: A scalar value of total entropy production in the system.
"""
function compute_global_entropy_production(J::Array{Float64, 2}, X::Array{Float64, 2})
    validate_dimensions_entropy(J, X)  # Reuse the existing validation function
    local_entropy_production = compute_local_entropy_production(J, X)
    global_entropy_production = sum(local_entropy_production)
    log_info("Global entropy production calculation complete with a total of $(global_entropy_production).")
    return global_entropy_production
end

# Reusing the existing helper function to validate dimensions
function validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})
    if size(J) != size(X)
        throw(DimensionMismatch("J and X must have the same dimensions."))
    end
end

# Reusing the existing local entropy production computation function
function compute_local_entropy_production(J::Array{Float64, 2}, X::Array{Float64, 2})
    validate_dimensions_entropy(J, X)
    entropy_production = sum(J .* X, dims=1)
    return vec(entropy_production)  # Convert matrix output to vector
end
