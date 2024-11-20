# NET.jl


## Table of Contents

### * Introduction
### * Onsager Relations
### * Entropy Production
### * Stohcastic Thermodynamics


## Introduction
NET.jl is a Julia library built for the sole purpose of working with elements of Non Equilibrium Thermodynamics(NET) principles in Julia. This library allows for the visulizaiton of certain phenomena, creation of Non Equilibrium thermodynamical systems, computation of Non Equilibrium Thermodynamics problems. The library mainly consists of 3 subsets, for the time being: Stohcastic Thermodynamics, Onsager Relations, Entropy Prodcution which are fundamental componenets of NET.

## Onsager Relations

Onsager relations are relations which express the relationship between fluxes and forces in a Thermodynamics system. Forces such as temperature gradients and concetration difference drive the fluxes, which are the determinanats of behaviors in the system, such as heat flux or diffusion (in the sense of movement of particles). 

There are 4 main files in the onsager relations subsection of the library being the linear, non Linear, time dependent and stohcastic onsager relations. In the following section the physics and programming of these classes will be explained.

### Linear Onsager Relations

In the class firstly the Onsager Matrix is initialized which is a necessity. Two types of Onsager Matrix exist in this context: 2D and 3D Onsager matrices. The main difference for these 2 types of Onsager Matrices is that a 3D Onsager Matrix accounts for local thermal propetries as well. 
```julia
mutable struct OnsagerMatrix
    L::Union{Array{Float64, 2}, Array{Float64, 3}}
end

```

Next the onsager relation function is initialized where the matrices are multiplied by the current forces of the system to get a flux vector (generally heat flux). The sole difference between the 2D and 3D matrice based onsager is established via the slicing technique, in the 3D array the each point has a individual "subarray" as well since their properties (coefficients etc) are different as well.

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

Non-linear onsager relations are another type of onsager relations which are frequently seen in thermodynamical system which are not in equilibrium. It's difference from linear onsager is that higher order dimensions are in play in the calcuations on top of the linear relations (e.g. includes variables such as second order coupling). In this type of onsager relations there is only one type of Onsager relations because since this involves higher order dimensions as well the local thermal properties are a must for the calculations due to those being some of the higher properties used when a Non-linear onsager calculation is involved. Below is the array

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
        J_linear[:, k] = L[:, :, k] * F[:, k] 
    end
    J_total .= J_linear 

    for F_nl in F_nl_list  
        for k in 1:num_points
            J_total[:, k] .+= L[:, :, k] * F_nl[:, k]
        end
    end

    return J_total
end

```
A similar dimension validation and visulization flux- like the one for the linear onsager- exists for this function as well

### Stohcastic Onsager Relations

Stohcastic onsager relations are spatial onsager relations which have a componenet of noise in the thermodynamical system as well. Again a 3D array is initalized for the computations of the stohcastic onsager relations

```julia
mutable struct StochasticOnsagerMatrix
    L::Array{Float64, 3}
end
```

In the computation function when the column vector J is generated the sole difference with the other computation functions is that in the computation process external noise is involved in the process as well. Again for the stroing of the data a 3D Array is used since stochastic onsager relations are a type of non linear onsager relations.

```julia
mutable struct StochasticOnsagerMatrix
    L::Array{Float64, 3}
end

```
After the storage array is defined the Stochastic flux computation function is defined. The main difference between the spatial linear onsager annd Stochastic onsager relations flux computation is know the external noise in the system is spesficially considered. 

```julia

function compute_stochastic_fluxes(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2})
    validate_dimensions_stochastic(L, F, noise)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)

    for k in 1:num_points
        J[:, k] = L.L[:, :, k] * F[:, k] + noise[:, k]
    end

    log_info("Stochastic flux computation complete for $num_points points with external noise.")
    return J
end


```
Here after the matrix multiplication noise is added onto the resultant array as well

### Time Dependent Onsager Relations

This is the last type of Onsager relations. It is another spesific type of Onsager Relations which considers a time based gradient rather than a distance based gradient.


```julia
function compute_time_dependent_fluxes(L::TimeDependentOnsagerMatrix, F::Array{Float64, 2})
    validate_dimensions_time_dependent(L, F)
    num_vars, num_times = size(F)
    J = zeros(num_vars, num_times)

    for t in 1:num_times
        J[:, t] = L.L[:, :, t] * F[:, t]
    end

    log_info("Time-dependent flux computation complete for $num_times time points.")
    return J
end

```

###### Important Note
The 2D matrices can be multiplied with a 3D matricie if the first 2 dimensions are the same because in that case the first two dimensions of the 3D matricie are treated as a 2D matricie and the third dim is depth so each depth is considered as a 2D matricie and the amount of depth is the amount of 2D matricie multiplication done. Below is the time spesific 3D matricie based function for time dependent Onsager calculation.


## Entropy Production

Entropy, per unit thermal energy which can't be used for work, production is a fundamental concept in non equilibrium thermodynamics; entropy production describes the rate at which entropy is produced. Like Onsager relations it is driven by thermodynamical forces. This concept describes the dissapation of energy, energy which has become unable to be used in work. In this subsection there are 4 main files

### Local Entropy Production

This is the entropy production for a local particle, it is calculated as the dot producted of the flux and force vectors

```julia

function compute_local_entropy_prozduction(J::Array{Float64, 2}, X::Array{Float64, 2})
    validate_dimensions_entropy(J, X)
    entropy_production = sum(J .* X, dims=1)  
    log_info("Local entropy production calculation complete for $(size(J, 2)) grid points.")
    return vec(entropy_production)  
end
```

In addition to these the class includes a dimension validation and a visulization function in the form of a heatmap

Dimension Validation:

```julia
function validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})
    if size(J) != size(X)
        throw(DimensionMismatch("J and X must have the same dimensions."))
    end
end
```

Visulization:

```julia
function visualize_local_entropy_production_heatmap(σ::Array{Float64, 1}, grid_points::Vector, title_name::String)

    σ_matrix = reshape(σ, 1, length(σ))
    
    p = heatmap(
        grid_points,
        1:1,  
        σ_matrix,
        xlabel="Grid Points",
        ylabel="Entropy Production",
        title=title_name,
        color=:plasma,
        colorbar=true
    )
    display(p) 
end
```

### Time Dependent Entropy Production

This class is for the computing of entropies based on fluxes and forces which are dependent on time. It consists of 2 functions one being the outer function containinig the ordinary differential equation and the solver and the inner one being the one which calculates the derivative of entropy production. There are three params for the outer functions. The outer function takes the time span, flux and force functions meanwhile the inner functions takes 2 parameters p which is a parameter to add coefficients to fluxes and t which is the time, used as input to the time dependent functions.

```julia

function analyze_time_dependent_entropy_production(J_func::Function, X_func::Function, t_span::Tuple{Float64, Float64})
    function dσ_dt(p,t)
        J = p * J_func(t)
        X = X_func(t)
        return sum(J .* X)
    end

    σ_0 = 0.0  
    prob = ODEProblem(dσ_dt, σ_0, t_span)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    return sol
end

```

### Entropy Production Rate

A function which computes the derivative of the entropy production based on entropy production in two points

```julia
function compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                              J2::Array{Float64, 2}, X2::Array{Float64, 2}, Δt::Float64)
    validate_dimensions_entropy(J1, X1)
    validate_dimensions_entropy(J2, X2)


    entropy_initial = sum(J1 .* X1, dims=1)
    entropy_final = sum(J2 .* X2, dims=1)


    Δentropy = sum(entropy_final) - sum(entropy_initial)
    entropy_rate = Δentropy / Δt

    log_info("Entropy rate calculation complete: rate = $(entropy_rate) per unit time.")
    return entropy_rate
end
```

### Entropy Optimization 

This is a custom built function in order to minimize the system's entropy based on Gradient Descent, firstly the derivatives are computed and then the results are inputed into gradient descent to find the most optimal entropy.

```julia

function entropy_gradient(J::Array{Float64, 2}, X::Array{Float64, 2})
    ∇J = X  
    ∇X = J  
    return (∇J, ∇X)
end

function gradient_descent(J::Array{Float64, 2}, X::Array{Float64, 2}, learning_rate::Float64, max_iterations::Int)
    for i in 1:max_iterations
        ∇J, ∇X = entropy_gradient(J, X)
        J -= learning_rate * ∇J
        X -= learning_rate * ∇X
        

        @show compute_entropy(J, X)
    end
    return (J, X)
end


```


### Stohcastic Thermodynamics

Stochastic thermodynamics is a framework that generalizes classical thermodynamics to systems influenced by randomness and fluctuations. It provides tools to study the interplay between energy, entropy, and noise, enabling the analysis of small-scale systems or processes driven by stochastic dynamics.

#### Stohcastic Driven Systems

This class allows for certain conditions of a stohcastic enviroment to be simulated which is necessary for an physically accurate. Namely this function allows for the evolution of properties such as the diffusion process, randmom movement due to thermal fluctuations; drift, systematic and non-random motion caused by one gradient- force- overweighing the other foces and pushing the particle to one direction; Stohcastic Growth, random growth and decay rates in particle properties; Nonlinear Dynamics of Noise, higher order interactions which affect behavior in the system.

The function created for this functionality is the "simulate_stochastic_system" in this function firstly the dimensions for realizations(particle trajectories) and number of points are done. Later a empty array is initialized. After that the first raw is set to x0 since it is the initial position. Later a nested loop is initalized to itirate over the data to calculate the trajectories., with outer being number of trajectories and inner being num of time step. 

In the loop for calculation there are two types SDE's in the loop the Euler Maruyama and the Milstein method to calculate the effect of stohcastic components on the system.

###### Euler-Maruyama

<img width="720" alt="Ekran Resmi 2024-11-20 12 22 12" src="https://github.com/user-attachments/assets/54f3a953-dcd6-4ce9-86d6-14cf9c124f8a">

The Euler-Maruyama eqution is the rate of change of the force plus the positon at initla time plus the amplitue of the noise times the root of the time times random noise, Wiener process to be spesific, which is a driver in the stohcastic system as well. This equation is more well suited to stohcastic systems with weak noise since non-linear interactions are better captured in Milstein and since in small noises the comptutationally expensive calculation of Milstein aren't necessary the usage of Eular-Maruyama is more pragmatic. Below is the code for Euler Maruyama

```julia

            if method == "Euler-Maruyama"

                trajectories[t, r] = trajectories[t - 1, r] +
                                     drift(trajectories[t - 1, r]) * dt +
                                     noise_amplitude(trajectories[t - 1, r]) * sqrt(dt) * ξ
```

###### Milstein

<img width="730" alt="Ekran Resmi 2024-11-20 12 22 33" src="https://github.com/user-attachments/assets/a27bc5b0-fdd6-4e1d-93a2-1c305354dcba">

The Milstein equation starts is the Euler Maruyama, but with added noise. In the Milstein equation on addition to the Euler Maruyama Equation the Derivative of the Noise amplitude is added, a higher order contribution of noise occurs making the relationship non-linear. Note that the derivation is based on Ito calculus.

```julia
 elseif method == "Milstein"

                g = noise_amplitude(trajectories[t - 1, r])
                g_prime = (noise_amplitude(trajectories[t - 1, r] + 1e-6) - g) / 1e-6 
                trajectories[t, r] = trajectories[t - 1, r] +
                                     drift(trajectories[t - 1, r]) * dt +
                                     g * sqrt(dt) * ξ +
                                     0.5 * g * g_prime * (ξ^2 - 1) * dt
```


### Fluctuation Theorems

Fluctuation theorems are theorems which describe forces in NET systems. The fluctuation theorems in the library are Crooks and Jarzynski's equations (so far).

#### Jarzynsiki's Theorem

The Jarzynski's equation is defined as the work values are first multiplied by thermodynamic beta to scale them relative to the system's temperature. These scaled values are then used as exponents in the exponential weight calculation, and the resulting exponential weights are averaged together to form the ensemble average. Finally, the natural logarithm of this average is multiplied by negative one and divided by thermodynamic beta to estimate the equilibrium free energy difference.

In the code the beta calculated is the inverse of thermal energy and used to define a relationship between thermal energy and temperature in the context of Jarzynski's equation. After that is completed to input into the equation the thermal noise with input beta is calculated. After these processes are complete the average work is computed via the Log-Sum-Exp trick(RealSofMax) in order to calculate it's stable average. Later, we calculate the free energy difference, which gives us systems ability to do useful work(negative okay, positive requires energy).

#### Crook's Theorem
Crooks' Theorem relates the probabilities of observing a certain amount of work in a forward process to the reverse process, stating that their ratio is governed by the work difference and the free energy change. It allows estimation of the equilibrium free energy difference (Δ𝐹) from the crossing point of the forward and reverse work distributions, where their probabilities are equal. 

The forward work is from A to B while the reverse work is from B to A- B denoting final and A denoting initial. The formula for thr Crook's Theorem is probability of forward work divided by probability of inverse work is equal to e to the power Beta(reverse thermal energy times work, forward, minus free energy. Below is the code


```julia
function crooks_theorem(work_forward::Vector{Float64}, work_reverse::Vector{Float64}, T::Float64)
    beta = 1 / (k_B * T)  

    log_ratios = Float64[]  
  
    for (w_f, w_r) in zip(work_forward, work_reverse)
        log_numerator = beta * w_f
        log_denominator = beta * w_r

        
        if isfinite(log_numerator) && isfinite(log_denominator) && log_denominator > -700
            push!(log_ratios, log_numerator - log_denominator)
        end
    end

    
    if isempty(log_ratios)
        println("Warning: No valid log-ratios found; returning NaN.")
        return NaN
    end


    mean_log_ratio = mean(log_ratios)
    mean_ratio = exp(mean_log_ratio)
    return mean_ratio
end

```
The final output is the probability ratio of forward to reverse work- likelihood of each type of work happening.


### Langevin Dynamics

Langevin dynamics is a method in statistical mechanics and thermodynamics used to model the motion of particles in a fluid or heat bath. It combines deterministic forces, such as those arising from a potential field, with stochastic thermal forces to account for the random interactions with the surrounding medium. The dynamics are governed by a balance between inertia, friction, and noise, providing a realistic representation of particle behavior in non-equilibrium systems.

In Langevin dynamics, a particle’s motion is influenced by deterministic forces from a potential field, frictional damping from its medium, and random thermal noise. These effects are modeled using the Langevin equation, which balances these forces to describe the system's trajectory. The thermal noise in this equation satisfies the fluctuation-dissipation theorem, ensuring that energy lost due to damping is replenished by thermal fluctuations, maintaining the system at equilibrium. This framework is used to study systems under non-equilibrium conditions, capturing both predictable and random behaviors.

In the code, the BAOAB integrator is applied as a numerical scheme to solve the Langevin equation. The purpose of BAOAB is to accurately and stably simulate Langevin dynamics while preserving equilibrium properties like the Boltzmann distribution. The method splits the simulation into three steps: handling bath interactions (damping and noise), applying deterministic forces, and adjusting particle positions. This ensures each aspect of the dynamics is addressed separately and correctly, making BAOAB an ideal choice for simulations requiring long-term stability and physical consistency.
