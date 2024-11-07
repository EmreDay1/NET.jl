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



