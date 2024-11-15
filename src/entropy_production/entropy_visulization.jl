using Plots

"""
    visualize_entropy_production(entropy_values::Array{Float64, 1}, title::String)

Creates a plot to visualize entropy production values.
- entropy_values: 1D array of entropy production values.
- title: Title for the plot.
"""
function visualize_entropy_production(entropy_values::Array{Float64, 1}, title::String)
    p = plot(entropy_values, title=title, xlabel="Index", ylabel="Entropy Production",
             legend=false, linewidth=2, color=:blue)
    display(p)
end
