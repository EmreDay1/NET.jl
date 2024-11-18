using Random
using LinearAlgebra
using Plots
using Thermodynamics  


function meshgrid(x, y)
    xg = repeat(x, 1, length(y))  
    yg = repeat(y', length(x), 1) 
    return xg, yg
end


"""
    onsager_relations_spatial(L::Array, X::Array)

Computes the fluxes (J) from the forces (X) using the phenomenological coefficients matrix (L),
while considering spatial dependence. The spatial domain is discretized, and forces and fluxes are computed
at each grid point.
- L: 3D array, with L[i,j,k] representing the (i,j) element of the phenomenological matrix at grid point k.
- X: 2D array, with X[j,k] representing the force on variable j at grid point k.
Returns: A 2D array of fluxes (J) at each grid point.
"""
function onsager_relations_spatial(L::Array, X::Array)
    num_variables, num_points = size(X) 
    J = zeros(num_variables, num_points) 

    for k in 1:num_points
        J[:, k] = L[:, :, k] * X[:, k]
    end

    return J  
end


"""
    thermoelectric_transport_spatial(L::Array, T::Array, μ::Array)

Computes the heat flux and charge flux for a thermoelectric system with spatial dependence.
- L: 3D array, with L[i,j,k] representing the phenomenological matrix at grid point k.
- T: 1D array, temperature gradient at each spatial point.
- μ: 1D array, chemical potential gradient at each spatial point.
Returns: A 2D array of fluxes (heat and charge flux) at each spatial point.
"""
function thermoelectric_transport_spatial(L::Array, T::Array, μ::Array)
    num_points = length(T)  #


    additional_force_1 = zeros(num_points)
    additional_force_2 = zeros(num_points) 


    forces = Matrix(hcat(T, μ, additional_force_1, additional_force_2)')  

    return onsager_relations_spatial(L, forces)
end


"""
    nonlinear_onsager_spatial(L::Array, X::Array, X2::Array)

Computes the fluxes from forces with nonlinear corrections, including spatial dependence.
- L: 3D array for spatial phenomenological coefficients.
- X: 2D array of forces (linear) at each spatial point.
- X2: 2D array of nonlinear force corrections at each spatial point.
Returns: A 2D array of fluxes (linear + nonlinear contributions) at each spatial point.
"""
function nonlinear_onsager_spatial(L::Array, X::Array, X2::Array)
    J_linear = onsager_relations_spatial(L, X) 
    J_nonlinear = onsager_relations_spatial(L, X2) 
    return J_linear + J_nonlinear  
end


domain_min = (-1.0, -1.0)
domain_max = (1.0, 1.0)
num_particles = 20
positions = [rand(2) .* (domain_max .- domain_min) .+ domain_min for _ in 1:num_particles]
positions = reduce(hcat, positions)


velocities = [randn(2) for _ in 1:num_particles]
velocities = reduce(hcat, velocities)


target_system_temperature = 300.0  
particle_temperatures = [target_system_temperature + randn() * 10.0 for _ in 1:num_particles]

collision_distance = 0.05
dt = 0.01
total_time = 10.0
num_steps = Int(total_time / dt)
constant_zone_temperature = 300.0  



particle_mass = 0.02  
specific_heat_capacity = 1005.0  
thermal_conductivity = 0.1 
diffusion_coefficient = 0.3  
soret_coefficient = 0.01  
heat_loss_coefficient = 0.05 


grid_size = (50, 50)
grid_x = LinRange(domain_min[1], domain_max[1], grid_size[1])
grid_y = LinRange(domain_min[2], domain_max[2], grid_size[2])
grid_temperatures = [target_system_temperature + 50.0 * sin(x * 2 * π) * cos(y * 2 * π) for x in grid_x, y in grid_y]


function calculate_temperature_gradient(grid_temperatures, x_idx, y_idx, grid_x, grid_y)
    dx = grid_x[2] - grid_x[1]
    dy = grid_y[2] - grid_y[1]
    
    x_idx = clamp(x_idx, 2, length(grid_x) - 1)
    y_idx = clamp(y_idx, 2, length(grid_y) - 1)

    dTdx = (grid_temperatures[x_idx+1, y_idx] - grid_temperatures[x_idx-1, y_idx]) / (2 * dx)
    dTdy = (grid_temperatures[x_idx, y_idx+1] - grid_temperatures[x_idx, y_idx-1]) / (2 * dy)

    return [dTdx, dTdy]
end


cold_zones = [0 for x in 1:grid_size[1], y in 1:grid_size[2]] 
cold_zone_temperature = 280.0 


for i in 1:div(grid_size[1], 2), j in 1:div(grid_size[2], 2)
    cold_zones[i, j] = 1  
end


function apply_heat_loss!(particle_temperatures, positions, grid_x, grid_y, cold_zones, cold_zone_temperature, heat_loss_coefficient, dt)
    for i in 1:num_particles
        
        x_idx = findfirst(x -> x >= positions[1, i], grid_x)
        y_idx = findfirst(y -> y >= positions[2, i], grid_y)

       
        if x_idx === nothing || y_idx === nothing
            continue
        end


        x_idx = clamp(x_idx, 1, length(grid_x))
        y_idx = clamp(y_idx, 1, length(grid_y))


        if cold_zones[x_idx, y_idx] == 1
            heat_loss = heat_loss_coefficient * (particle_temperatures[i] - cold_zone_temperature) * dt
            particle_temperatures[i] -= heat_loss
        end
    end
end


function apply_fluxes!(positions, velocities, particle_temperatures, grid_temperatures, grid_x, grid_y, dt, particle_mass, specific_heat_capacity)
    for i in 1:num_particles

        x_idx = findfirst(x -> x >= positions[1, i], grid_x)
        y_idx = findfirst(y -> y >= positions[2, i], grid_y)

        if x_idx === nothing || y_idx === nothing
            continue 
        end

        x_idx = clamp(x_idx, 2, length(grid_x) - 1)
        y_idx = clamp(y_idx, 2, length(grid_y) - 1)

       
        temp_gradient = calculate_temperature_gradient(grid_temperatures, x_idx, y_idx, grid_x, grid_y)


        heat_flux = -thermal_conductivity * temp_gradient
        heat_transfer = norm(heat_flux) * dt


        delta_T = heat_transfer / (specific_heat_capacity * particle_mass)
        particle_temperatures[i] += delta_T

        velocities[:, i] .+= temp_gradient * dt * 0.01

      
        velocities[:, i] .+= randn(2) * 0.01 * sqrt(0.3)
    end
end



function handle_collisions!(positions, velocities, particle_temperatures)
    for i in 1:num_particles
        for j in i+1:num_particles
            distance = norm(positions[:, i] .- positions[:, j])
            if distance < collision_distance

                repulsion_force = (positions[:, i] .- positions[:, j]) / distance^2
                velocities[:, i] .+= repulsion_force * 0.01
                velocities[:, j] .-= repulsion_force * 0.01


                overlap_distance = collision_distance - distance
                separation_vector = normalize(positions[:, i] .- positions[:, j]) * overlap_distance * 0.5
                positions[:, i] .+= separation_vector
                positions[:, j] .-= separation_vector

 
                vi, vj = velocities[:, i], velocities[:, j]
                velocities[:, i], velocities[:, j] = vj, vi

              
                temp_diff = particle_temperatures[i] - particle_temperatures[j]
                heat_transfer = 0.05 * temp_diff
                particle_temperatures[i] -= heat_transfer
                particle_temperatures[j] += heat_transfer
            end
        end
    end
end

function apply_constant_temperature_zones!(particle_temperatures, positions, grid_x, grid_y, constant_temperature_zones, constant_zone_temperature, specific_heat_capacity, particle_mass, dt)
    for i in 1:length(particle_temperatures)
     
        x_idx = findfirst(x -> x >= positions[1, i], grid_x)
        y_idx = findfirst(y -> y >= positions[2, i], grid_y)

        if x_idx === nothing || y_idx === nothing
            continue  # Skip if no valid grid index
        end


        x_idx = clamp(x_idx, 1, length(grid_x))
        y_idx = clamp(y_idx, 1, length(grid_y))


        if constant_temperature_zones[x_idx, y_idx] == 1
  
            delta_T = (constant_zone_temperature - particle_temperatures[i]) * 0.1 * dt
            particle_temperatures[i] += delta_T
        end
    end
end

function handle_boundaries!(positions, velocities, domain_min, domain_max)
    num_particles = size(positions, 2)
    for i in 1:num_particles
        for dim in 1:2  
            if positions[dim, i] < domain_min[dim]
                positions[dim, i] = domain_min[dim]  
                velocities[dim, i] = -velocities[dim, i]  
            elseif positions[dim, i] > domain_max[dim]
                positions[dim, i] = domain_max[dim] 
                velocities[dim, i] = -velocities[dim, i]  
            end
        end
    end
end




function update_positions_and_velocities!(positions, velocities, particle_temperatures, dt, domain_min, domain_max, grid_temperatures, grid_x, grid_y, constant_temperature_zones, constant_zone_temperature, specific_heat_capacity, particle_mass)
    apply_fluxes!(positions, velocities, particle_temperatures, grid_temperatures, grid_x, grid_y, dt, particle_mass, specific_heat_capacity)

    apply_heat_loss!(particle_temperatures, positions, grid_x, grid_y, cold_zones, cold_zone_temperature, heat_loss_coefficient, dt)


  
    apply_constant_temperature_zones!(particle_temperatures, positions, grid_x, grid_y, constant_temperature_zones, constant_zone_temperature, specific_heat_capacity, particle_mass, dt)

    for i in 1:num_particles
 
        positions[:, i] .= positions[:, i] .+ velocities[:, i] .* dt
    end

    
    handle_boundaries!(positions, velocities, domain_min, domain_max)
end


function kinetic_energy(velocities)
    total_kinetic_energy = 0.0
    for i in 1:num_particles
        speed = norm(velocities[:, i])
        total_kinetic_energy += 0.5 * particle_mass * speed^2
    end
    return total_kinetic_energy
end


function calculate_heat_transfer(temp1, temp2, specific_heat, mass)
    heat_transfer = specific_heat * mass * (temp1 - temp2)
    return heat_transfer
end



anim = @animate for step in 1:num_steps
   
    handle_collisions!(positions, velocities, particle_temperatures)
    
  
    update_positions_and_velocities!(positions, velocities, particle_temperatures, dt, domain_min, domain_max, grid_temperatures, grid_x, grid_y, constant_temperature_zones, constant_zone_temperature, specific_heat_capacity, particle_mass)


    heatmap(grid_x, grid_y, grid_temperatures, color=cgrad(:inferno), alpha=0.5, xlabel="x", ylabel="y", title="Grid Temperature Field")

    scatter!(positions[1, :], positions[2, :], markersize=5, marker_z=particle_temperatures, color=cgrad(:plasma), legend=false)

    
    scatter!([positions[1, 1]], [positions[2, 1]], markersize=5, marker_z=[particle_temperatures[1]], 
             color=:plasma, markerstrokewidth=2, markerstroke=:orange, legend=false)

    x_idx = findfirst(x -> x >= positions[1, 1], grid_x)
    y_idx = findfirst(y -> y >= positions[2, 1], grid_y)
    if x_idx !== nothing && y_idx !== nothing
        temp_gradient = calculate_temperature_gradient(grid_temperatures, x_idx, y_idx, grid_x, grid_y)
        heat_flux = -thermal_conductivity * temp_gradient
        heat_flux_value = norm(heat_flux)
    else
        heat_flux_value = 0.0
    end


    arrow_scale = 0.005  
    quiver!([positions[1, 1]], [positions[2, 1]], quiver=([arrow_scale * heat_flux[1]], [arrow_scale * heat_flux[2]]), 
            arrow=:true, linewidth=1, color=:black, legend=false)

    # Set axis limits
    xlims!(domain_min[1], domain_max[1])
    ylims!(domain_min[2], domain_max[2])

    
    total_kinetic_energy = kinetic_energy(velocities)
    title!("Kinetic Energy: $total_kinetic_energy J | Particle 1 Heat Flux: $(round(heat_flux_value, digits=3)) W")
end


function apply_constant_temperature_zones!(particle_temperatures, positions, grid_x, grid_y, constant_temperature_zones, constant_zone_temperature, specific_heat_capacity, particle_mass, dt)
    for i in 1:num_particles
     
        x_idx = findfirst(x -> x >= positions[1, i], grid_x)
        y_idx = findfirst(y -> y >= positions[2, i], grid_y)

        if x_idx === nothing || y_idx === nothing
            continue 
        end

        x_idx = clamp(x_idx, 1, length(grid_x))
        y_idx = clamp(y_idx, 1, length(grid_y))

      
        if constant_temperature_zones[x_idx, y_idx] == 1
            delta_T = (constant_zone_temperature - particle_temperatures[i]) * 0.1 * dt
            heat_transfer = calculate_heat_transfer(particle_temperatures[i], constant_zone_temperature, specific_heat_capacity, particle_mass)
            particle_temperatures[i] += delta_T
        end
    end
end



gif(anim, "sim.gif")


function apply_constant_temperature_zones!(particle_temperatures, positions, grid_x, grid_y, constant_temperature_zones, constant_zone_temperature, specific_heat_capacity, particle_mass, dt)
    for i in 1:num_particles

        x_idx = findfirst(x -> x >= positions[1, i], grid_x)
        y_idx = findfirst(y -> y >= positions[2, i], grid_y)

        if x_idx === nothing || y_idx === nothing
            continue  
        end

        x_idx = clamp(x_idx, 1, length(grid_x))
        y_idx = clamp(y_idx, 1, length(grid_y))

        # If the particle is in a constant temperature zone, adjust its temperature
        if constant_temperature_zones[x_idx, y_idx] == 1
            delta_T = (constant_zone_temperature - particle_temperatures[i]) * 0.1 * dt
            heat_transfer = calculate_heat_transfer(particle_temperatures[i], constant_zone_temperature, specific_heat_capacity, particle_mass)
            particle_temperatures[i] += delta_T
        end
    end
end
