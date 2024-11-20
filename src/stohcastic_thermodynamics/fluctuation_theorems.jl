using Random

# Physical constants
const k_B = 1.380649e-23  # Boltzmann constant (J/K)

# ------------------------
# Utility Functions
# ------------------------


function thermal_noise(T::Float64)
    return sqrt(k_B * T) * 1e-2 
end

# Custom mean function
function mean(values::Vector{Float64})
    return sum(values) / length(values)
end

# ------------------------
# Jarzynski Equality (Using External Work Values)
# ------------------------

"""
    jarzynski_equality(work_values, T)

Estimates the free energy difference ΔF using Jarzynski's equality with externally supplied
work values, applying the Log-Sum-Exp trick for stability.

# Arguments
- `work_values::Vector{Float64}` : Array of externally provided work values in Joules.
- `T::Float64` : Temperature in Kelvin.

# Returns
- Estimated free energy difference `ΔF` in Joules.
"""
function jarzynski_equality(work_values::Vector{Float64}, T::Float64)
    beta = 1 / (k_B * T)


    max_work = maximum(-beta * work_values)
    avg_exp_work = exp(-max_work) * mean(exp.(-beta * work_values .+ max_work))
    
   
    if avg_exp_work <= 0 || isnan(avg_exp_work)
        println("Warning: avg_exp_work is non-positive or NaN; adjusting to avoid -Inf result.")
        avg_exp_work = 1e-15  # Small positive value
    end
    
   
    ΔF_estimated = -log(avg_exp_work) / beta
    return ΔF_estimated
end

# ------------------------
# Crooks Fluctuation Theorem (Using External Work Values)
# ------------------------

"""
    crooks_theorem(work_forward, work_reverse, T)

Calculates the probability ratio P_F(W) / P_R(-W) for forward and reverse work values 
to test Crooks' theorem with externally provided data.

# Arguments
- `work_forward::Vector{Float64}` : Array of work values in the forward direction.
- `work_reverse::Vector{Float64}` : Array of work values in the reverse direction.
- `T::Float64` : Temperature in Kelvin.

# Returns
- Mean of the ratio, which should be close to 1 if Crooks' theorem holds.
"""
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
