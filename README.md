# NET.jl


## Table of Contents

### * Introduction
### * Onsager Relations
### * Entropy Production
### * Stohcastic Thermodynamics
### * Example Usage 1: Non Equilibrium Particle Simulation


## Introduction
NET.jl is a Julia library built for the sole purpose of working with elements of Non Equilibrium Thermodynamics(NET) principles in Julia. This library allows for the visulizaiton of certain phenomena, creation of Non Equilibrium thermodynamical systems, computation of Non Equilibrium Thermodynamics problems. The library mainly consists of 3 subsets, for the time being: Stohcastic Thermodynamics, Onsager Relations, Entropy Prodcution which are fundamental componenets of NET.

## Onsager Relations

Onsager relations are relations which express the relationship between fluxes and forces in a Thermodynamics system. Forces such as temperature gradients and concetration difference drive the fluxes, which are the determinanats of behaviors in the system, such as heat flux or diffusion (in the sense of movement of particles). 

There are 6 main files in the onsager relations subsection of the library being the spatial, linear, non Linear, time dependent, stohcastic onsager relations and coupled transport. In the following section the physics and programming of these classes will be explained.

#### Linear Onsager Relations

In the class firstly the Onsager Matrix is initialized which is a necessity. Two types of Onsager Matrix exist in this context: 2D and 3D Onsager matrices. The main difference for these 2 types of Onsager Matrices is that a 3D Onsager Matrix accounts for local thermal propetries as well. 
```julia
mutable struct OnsagerMatrix
    L::Union{Array{Float64, 2}, Array{Float64, 3}}
end

```

Next the onsager relation function is initialized where the matrices are multiplied by the current forces of the system to get a flux vector (generally heat flux). The sole difference between the 2D and 3D matrice based onsager is established via the slicing technique, in the 3D matrice the each point has a individual submatrice as well since their properties (coefficients etc) are different as well.

```julia
function compute_fluxes(L::OnsagerMatrix, F::Array{Float64, 2})
    validate_dimensions(L, F)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)
    for k in 1:num_points
        if ndims(L.L) == 2
            J[:, k] = L.L * F[:, k]
        else
            J[:, k] = L.L[:, :, k] * F[:, k]
        end
    end
    log_info("Flux computation complete for $num_points points.")
    return J
end
```
There are additional functions in the library in order to validate dimensions and do a quick visulization; they won't be explained here since they aren't completely related with the physics of the library- this will keep through for all files containing classes.

#### Non-Linear Onsager Relations


