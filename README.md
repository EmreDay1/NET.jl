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

There are 5 main files in the onsager relations subsection of the library being the linear, non Linear, time dependent, stohcastic onsager relations and coupled transport. In the following section the physics and programming of these classes will be explained.

### Linear Onsager Relations

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
There are additional functions in the library in order to validate dimensions and do a quick visulization; they won't be explained here since they aren't completely related with the physics of the library- this will keep through for all files containing classes. Also note that the 2D matrix is used to express simple linear onsager relationships, meanwhile the 3D Matrix linear onsager expresses another phenomena known as spatial onsager relations.

### Non-Linear Onsager Relations

Non-linear onsager relations are another type of onsager relations which are frequently seen in thermodynamical system which are not in equilibrium. It's difference from linear onsager is that higher order dimensions are in play in the calcuations on top of the linear relations (e.g. includes variables such as second order coupling). In this type of onsager relations there is only one type of Onsager relations because since this involves higher order dimensions as well the local thermal properties are a must for the calculations due to those being some of the higher properties used when a Non-linear onsager calculation is involved. Below is the matricie

```julia
mutable struct NonLinearOnsagerMatrix
    L::Array{Float64, 3}
end
```
The next structure in the code is the NonLinearFluxSystem structure which is a basic stroage structure to represent the forces in the system 
```julia
mutable struct NonLinearFluxSystem
    matrix::NonLinearOnsagerMatrix
    linear_forces::Array{Float64, 2}
    nonlinear_forces::Vector{Array{Float64, 2}}
end
```

Then after these structures are created again the flux computation function is yet again created, but only a 3D version of the array exist here; in the function firstly the linear fluxes are computed and broadcasted upon the total Onsager flux array and later the same process is repeated for the non-linear fluxes, below is the structure.

```julia
function compute_total_fluxes(system::NonLinearFluxSystem)
    L = system.matrix.L
    F = system.linear_forces
    F_nl_list = system.nonlinear_forces

    validate_dimensions_multinonlinear(system)

    num_vars, num_points = size(F)
    J_linear = zeros(num_vars, num_points)
    J_total = zeros(num_vars, num_points)

    for k in 1:num_points
        J_linear[:, k] = L[:, :, k] * F[:, k]  # Linear contribution
    end
    J_total .= J_linear  # Initialize total flux

    for F_nl in F_nl_list  # Add non-linear contributions
        for k in 1:num_points
            J_total[:, k] .+= L[:, :, k] * F_nl[:, k]
        end
    end

    return J_total
end

```
A similar dimension validation and visulization flux- like the one for the linear onsager- exists for this function as well





