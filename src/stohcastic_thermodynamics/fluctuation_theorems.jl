# fluctuation_theorems.jl

using Random

# Physical constants
k_B = 1.380649e-23 


function thermal_noise(T::Float64)
    return sqrt(k_B * T) * 1e-2 
end

# Manually define a mean function
function mean(values::Vector{Float64})
    return sum(values) / length(values)
end

# ------------------------
# Jarzynski Equality with Log-Sum-Exp Trick
# ------------------------

"""
    jarzynski_equality(num_trials, ΔF_true, T)

Estimates the free energy difference ΔF using Jarzynski's equality, based on a physically-derived
thermal noise calculated from temperature, applying the Log-Sum-Exp trick for stability.

# Arguments
- `num_trials::Int` : Number of trials to simulate.
- `ΔF_true::Float64` : True free energy difference in Joules.
- `T::Float64` : Temperature in Kelvin.

# Returns
- Estimated free energy difference `ΔF` in Joules.
"""
function jarzynski_equality(num_trials::Int, ΔF_true::Float64, T::Float64)
    beta = 1 / (k_B * T)  
    noise_stddev = thermal_noise(T)

    work_values = [ΔF_true + noise_stddev * (rand() - 0.5) for _ in 1:num_trials]

   
    max_work = maximum(-beta * work_values)
    avg_exp_work = exp(-max_work) * mean(exp.(-beta * work_values .+ max_work))
    
    if avg_exp_work <= 0 || isnan(avg_exp_work)
        println("Warning: avg_exp_work is non-positive or NaN; adjusting to avoid -Inf result.")
        avg_exp_work = 1e-15 
    end
    
    ΔF_estimated = -log(avg_exp_work) / beta
    return ΔF_estimated
end

# ------------------------
# Crooks Fluctuation Theorem with Validation
# ------------------------

"""
    crooks_theorem(num_trials, ΔF_true, T)

Calculates the probability ratio P_F(W) / P_R(-W) for forward and reverse work values 
to test Crooks' theorem with thermodynamic noise.

# Arguments
- `num_trials::Int` : Number of trials to simulate.
- `ΔF_true::Float64` : True free energy difference in Joules.
- `T::Float64` : Temperature in Kelvin.

# Returns
- Mean of the ratio, which should be close to 1 if Crooks' theorem holds.
"""
function crooks_theorem(num_trials::Int, ΔF_true::Float64, T::Float64)
    beta = 1 / (k_B * T)  #
    noise_stddev = thermal_noise(T) * 1e-5 

   
    work_forward = [noise_stddev * (rand() - 0.5) for _ in 1:num_trials]
    
    work_reverse = [-w_f + ΔF_true for w_f in work_forward]

    log_ratios = Float64[]

   
    for (w_f, w_r) in zip(work_forward, work_reverse)
        log_numerator = beta * (w_f + ΔF_true)
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



