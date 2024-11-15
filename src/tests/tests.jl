""" 


coupled transport test


"""

L_data = rand(3, 3, 50) 
L = CoupledTransportMatrix(L_data)
F = rand(3, 50)  


J = compute_coupled_fluxes(L, F)


visualize_coupled_fluxes(J, "Coupled Transport Fluxes")


"""


spatial onsager test


"""


L_data = rand(4, 4, 100) 
L = SpatialOnsagerMatrix(L_data)
F = rand(4, 100) 


J = compute_spatial_fluxes(L, F)


visualize_spatial_fluxes(J, "Spatial Onsager Fluxes")


L_data = rand(4, 4, 100) 
L = SpatialOnsagerMatrix(L_data)
F = rand(4, 100)


J = compute_spatial_fluxes(L, F)

visualize_spatial_fluxes(J, "Spatial Onsager Fluxes")

"""


time dependent onsager test


"""

L_data = rand(3, 3, 50) 
L = TimeDependentOnsagerMatrix(L_data)


F = rand(3, 50)


time_points = collect(1:50)


J = compute_time_dependent_fluxes(L, F)


visualize_time_dependent_fluxes(J, time_points, "Time-Dependent Onsager Fluxes")


log_info("Visualization of time-dependent Onsager fluxes completed.")

"""


stochastic onsager test


"""


num_vars = 4  
num_points = 100  

L_data = rand(num_vars, num_vars, num_points)  
L = StochasticOnsagerMatrix(L_data)


F = rand(num_vars, num_points)


noise_level = 0.1


points = collect(1:num_points)


J = compute_stochastic_fluxes(L, F, noise_level)


visualize_stochastic_fluxes(J, points, "Stochastic Onsager Fluxes")


log_info("Visualization of stochastic Onsager fluxes completed.")

"""


local entropy production test


"""

num_vars = 3 
num_points = 100  


J = rand(num_vars, num_points) 
X = rand(num_vars, num_points) 


grid_points = collect(1:num_points)

   
σ = compute_local_entropy_production(J, X)

  
visualize_local_entropy_production_heatmap(σ, grid_points, "Local Entropy Production Heatmap")

log_info("Visualization of local entropy production completed.")

"""


langevin_dynamics


"""

# Double-well potential and gradient
function double_well_potential(x::Float64)
    return x^4 - 2 * x^2  # Simple double-well potential
end

function grad_double_well_potential(x::Float64)
    return 4 * x^3 - 4 * x  # Gradient of the double-well potential
end

"""
    test_langevin_dynamics()

Runs a test of Langevin dynamics with a double-well potential and generates improved plots for
visualizing the particle's position and velocity over time.
"""
function test_langevin_dynamics()
    # Simulation parameters
    num_steps = 5000        # Number of time steps
    dt = 0.01               # Time step (s)
    gamma = 1.0             # Damping coefficient (kg/s)
    T = 300.0               # Temperature (K)
    mass = 1e-21            # Mass of the particle (kg)

    # Run the simulation
    time, positions, velocities = langevin_dynamics(
        num_steps, dt, gamma, T, mass, double_well_potential, grad_double_well_potential
    )

    # Plot position over time
    plot(time, positions, xlabel="Time (s)", ylabel="Position (m)", label="Position",
         title="Langevin Dynamics in Double-Well Potential", lw=2)
    scatter!(time[1:500:end], positions[1:500:end], label="Sampled Positions", alpha=0.6)

    # Phase-space plot (Position vs. Velocity)
    plot(positions, velocities, xlabel="Position (m)", ylabel="Velocity (m/s)", label="Phase Space",
         title="Phase Space of Langevin Dynamics", lw=2)
end

# Run the test example
test_langevin_dynamics()

"""


fluctuation_theorems


"""

function test_fluctuation_theorems()
    # Set parameters for the simulation
    T = 300.0              # Temperature in Kelvin (room temperature)
    ΔF_true = -1e-19       # Increased free energy difference in Joules for stability
    num_trials = 10000     # Number of trials to simulate for averaging

    # 1. Test Jarzynski's Equality
    ΔF_estimated = jarzynski_equality(num_trials, ΔF_true, T)
    println("Jarzynski Equality Test:")
    println("True ΔF (Joules): ", ΔF_true)
    println("Estimated ΔF from Jarzynski Equality (Joules): ", ΔF_estimated)

    # 2. Test Crooks' Theorem
    mean_ratio = crooks_theorem(num_trials, ΔF_true, T)
    println("\nCrooks Theorem Test:")
    println("Mean Ratio (should be close to 1 if theorem holds): ", mean_ratio)
end

# Run the test example when the script is executed
test_fluctuation_theorems()



