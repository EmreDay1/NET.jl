using Random
using Plots

# Physical constants
k_B = 1.380649e-23  # Boltzmann constant (J/K)

"""
    langevin_dynamics(num_steps, dt, gamma, T, mass, potential, grad_potential)

Simulates Langevin dynamics for a particle in a potential field with thermal noise, damping, and deterministic forces.

# Arguments
- `num_steps::Int`: Number of simulation steps.
- `dt::Float64`: Time step for the simulation (seconds).
- `gamma::Float64`: Damping coefficient (kg/s).
- `T::Float64`: Temperature (Kelvin).
- `mass::Float64`: Mass of the particle (kg).
- `potential::Function`: Function defining the potential energy field.
- `grad_potential::Function`: Function defining the gradient (force) of the potential.

# Returns
- `time::Vector{Float64}`: Time vector.
- `positions::Vector{Float64}`: Position of the particle at each time step.
- `velocities::Vector{Float64}`: Velocity of the particle at each time step.
"""
function langevin_dynamics(num_steps::Int, dt::Float64, gamma::Float64, T::Float64, mass::Float64,
                           potential::Function, grad_potential::Function)

    time = collect(0:dt:(num_steps - 1) * dt)
    positions = zeros(Float64, num_steps)
    velocities = zeros(Float64, num_steps)


    sqrt_2kBT_gamma_dt = sqrt(2 * k_B * T * gamma / mass * dt)

  
    positions[1] = 0.5  
    velocities[1] = 0.0 


    for i in 2:num_steps
        x = positions[i - 1]
        v = velocities[i - 1]

        
        force = -grad_potential(x)
        random_force = sqrt_2kBT_gamma_dt * randn()


        v_new = v + dt * (force / mass - gamma * v / mass) + random_force
        x_new = x + dt * v_new

        velocities[i] = v_new
        positions[i] = x_new
    end

    return time, positions, velocities
end


