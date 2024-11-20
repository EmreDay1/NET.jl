using Random
using Plots

"""
    compute_msd(positions, time)

Computes the mean squared displacement (MSD) for a particle given its positions over time.

# Arguments
- `positions::Matrix{Float64}`: Positions of the particle over time for each dimension.
- `time::Vector{Float64}`: Time vector corresponding to the positions.

# Returns
- `msd::Vector{Float64}`: Mean squared displacement at each time step.
"""
function compute_msd(positions::Matrix{Float64}, time::Vector{Float64})
    num_steps, dims = size(positions)
    msd = zeros(Float64, num_steps)
    for t in 1:num_steps
        displacements = positions[t:num_steps, :] .- positions[1:num_steps - t + 1, :]
        msd[t] = mean(sum(displacements.^2, dims=2))
    end
    return msd
end
