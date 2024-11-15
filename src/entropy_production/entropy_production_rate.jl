using LinearAlgebra
using Logging

"""
    compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                         J2::Array{Float64, 2}, X2::Array{Float64, 2}, Δt::Float64) -> Float64

Calculates the rate of irreversible entropy production by comparing entropy production at two different states over time Δt.
- J1, X1: 2D arrays of fluxes and thermodynamic forces at the initial state.
- J2, X2: 2D arrays of fluxes and thermodynamic forces at the final state.
- Δt: Time interval between the two states.
Returns: The rate of entropy production in the system per unit time.
"""
function compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                              J2::Array{Float64, 2}, X2::Array{Float64, 2}, Δt::Float64)
    validate_dimensions_entropy(J1, X1)
    validate_dimensions_entropy(J2, X2)

    # Calculate local entropy production at initial and final states
    entropy_initial = sum(J1 .* X1, dims=1)
    entropy_final = sum(J2 .* X2, dims=1)

    # Calculate the change in entropy production and normalize by time
    Δentropy = sum(entropy_final) - sum(entropy_initial)
    entropy_rate = Δentropy / Δt

    log_info("Entropy rate calculation complete: rate = $(entropy_rate) per unit time.")
    return entropy_rate
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
