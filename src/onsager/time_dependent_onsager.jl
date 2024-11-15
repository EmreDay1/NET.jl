using LinearAlgebra
using Logging
using Plots

"""
    struct TimeDependentOnsagerMatrix

Represents a time-dependent phenomenological coefficient matrix for Onsager relations.
- L: 3D array where L[:,:,t] is the matrix of coefficients at time point t.
"""
mutable struct TimeDependentOnsagerMatrix
    L::Array{Float64, 3}
end

"""
    compute_time_dependent_fluxes(L::TimeDependentOnsagerMatrix, F::Array{Float64, 2}) -> Array{Float64, 2}

Computes the fluxes (J) at each time point based on the time-dependent forces (F) using the time-dependent Onsager matrix (L).
- L: A `TimeDependentOnsagerMatrix` type representing time-dependent phenomenological coefficients.
- F: 2D array representing forces at each time step, where each column is a force vector at a time point.
Returns: A 2D array of computed fluxes at each time point.
"""
function compute_time_dependent_fluxes(L::TimeDependentOnsagerMatrix, F::Array{Float64, 2})
    validate_dimensions_time_dependent(L, F)
    num_vars, num_times = size(F)
    J = zeros(num_vars, num_times)

    for t in 1:num_times
        J[:, t] = L.L[:, :, t] * F[:, t]
    end

    log_info("Time-dependent flux computation complete for $num_times time points.")
    return J
end

"""
    validate_dimensions_time_dependent(L::TimeDependentOnsagerMatrix, F::Array{Float64, 2})

Ensures that the input dimensions of the matrix L and the force array F match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_time_dependent(L::TimeDependentOnsagerMatrix, F::Array{Float64, 2})
    num_vars, num_times = size(F)
    if size(L.L, 1) != num_vars || size(L.L, 3) != num_times
        throw(DimensionMismatch("L must match the number of variables and time points in F."))
    end
end

"""
    visualize_time_dependent_fluxes(J::Array{Float64, 2}, time_points::Vector, title_name::String)

Creates a time series plot of the fluxes for easy interpretation.
- J: 2D array of fluxes at each time point.
- time_points: Vector representing the time points for the x-axis.
"""
function visualize_time_dependent_fluxes(J::Array{Float64, 2}, time_points::Vector, title_name::String)
    num_vars = size(J, 1)
    p = plot(xlabel="Time", ylabel="Flux Magnitude", title=title_name, lw=2, legend=:topright)

    for i in 1:num_vars
        plot!(p, time_points, J[i, :], label="Flux $i")
    end

    display(p)
end

