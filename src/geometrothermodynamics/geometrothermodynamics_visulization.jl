using Plots
using LinearAlgebra

"""
    plot_thermodynamic_potential(time, potential_values)

Visualize the evolution of a thermodynamic potential over time.

# Parameters
- `time`: A vector of time points.
- `potential_values`: A vector of potential values corresponding to time.
"""
function plot_thermodynamic_potential(time::Vector{Float64}, potential_values::Vector{Float64})
    plot(time, potential_values, xlabel="Time", ylabel="Potential", lw=2, label="Thermodynamic Potential")
end

"""
    plot_legendre_transformation(time, original_values, transformed_values)

Visualize the original and Legendre-transformed potentials over time.

# Parameters
- `time`: A vector of time points.
- `original_values`: A vector of original potential values.
- `transformed_values`: A vector of Legendre-transformed potential values.
"""
function plot_legendre_transformation(time::Vector{Float64}, original_values::Vector{Float64}, transformed_values::Vector{Float64})
    plot(time, original_values, xlabel="Time", ylabel="Potential", lw=2, label="Original Potential")
    plot!(time, transformed_values, lw=2, label="Legendre Transformed Potential")
end

"""
    plot_ricci_curvature(time, curvature_values)

Visualize the Ricci curvature of the thermodynamic manifold over time.

# Parameters
- `time`: A vector of time points.
- `curvature_values`: A vector of Ricci curvature values corresponding to time.
"""
function plot_ricci_curvature(time::Vector{Float64}, curvature_values::Vector{Float64})
    plot(time, curvature_values, xlabel="Time", ylabel="Ricci Curvature", lw=2, label="Ricci Curvature")
end

"""
    plot_thermodynamic_metrics(metrics)

Visualize the thermodynamic metric tensor as a heatmap.

# Parameters
- `metrics`: A matrix representing the thermodynamic metric tensor.
"""
function plot_thermodynamic_metrics(metrics::Matrix{Float64})
    heatmap(metrics, xlabel="Variables", ylabel="Variables", color=:viridis, title="Thermodynamic Metric")
end

"""
    plot_3d_potential_surface(X, Y, Z)

Visualize the potential surface as a 3D plot.

# Parameters
- `X`: A matrix of x-coordinates.
- `Y`: A matrix of y-coordinates.
- `Z`: A matrix of potential values corresponding to (X, Y).
"""
function plot_3d_potential_surface(X::Matrix{Float64}, Y::Matrix{Float64}, Z::Matrix{Float64})
    surface(X, Y, Z, xlabel="X1", ylabel="X2", zlabel="Potential", title="Potential Surface")
end
