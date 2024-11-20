using Random
using Plots

k_B = 1.380649e-23  # Boltzmann constant (J/K)

"""
    langevin_dynamics(num_steps, dt, gamma, T, mass, potential, grad_potential; dims=1)

Simulates Langevin dynamics for a particle in a potential field with thermal noise, damping, and deterministic forces.

# Arguments
- `num_steps::Int`: Number of simulation steps.
- `dt::Float64`: Time step for the simulation (seconds).
- `gamma::Float64`: Damping coefficient (kg/s).
- `T::Float64`: Temperature (Kelvin).
- `mass::Float64`: Mass of the particle (kg).
- `potential::Function`: Function defining the potential energy field.
- `grad_potential::Function`: Function defining the gradient (force) of the potential.
- `dims::Int`: Dimensionality of the simulation (default: 1).

# Returns
- `time::Vector{Float64}`: Time vector.
- `positions::Matrix{Float64}`: Position of the particle at each time step for each dimension.
- `velocities::Matrix{Float64}`: Velocity of the particle at each time step for each dimension.
"""
function langevin_dynamics(num_steps::Int, dt::Float64, gamma::Float64, T::Float64, mass::Float64,
                           potential::Function, grad_potential::Function; dims::Int = 1)

    time = collect(0:dt:(num_steps - 1) * dt)
    positions = zeros(Float64, num_steps, dims)
    velocities = zeros(Float64, num_steps, dims)

    sigma = sqrt(2 * k_B * T * gamma / mass)
    friction_factor = exp(-gamma * dt / mass)

    positions[1, :] .= randn(dims)
    velocities[1, :] .= sigma * randn(dims)

    for i in 2:num_steps
        x = positions[i - 1, :]
        v = velocities[i - 1, :]

        random_noise = sigma * sqrt(1 - friction_factor^2) * randn(dims)
        v_half = v .* friction_factor .+ 0.5 * dt .* (-grad_potential(x) ./ mass) .+ random_noise
        x_new = x .+ dt .* v_half
        v_new = v_half .+ 0.5 * dt .* (-grad_potential(x_new) ./ mass)

        positions[i, :] .= x_new
        velocities[i, :] .= v_new
    end

    return time, positions, velocities
end

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
