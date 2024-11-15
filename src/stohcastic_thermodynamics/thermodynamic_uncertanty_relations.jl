using Statistics

"""
    calculate_tur(entropy_production_rate, variance)

Calculates the thermodynamic uncertainty relation (TUR) bound for a system.

# Arguments
- entropy_production_rate::Float64: Steady-state entropy production rate.
- variance::Float64: Variance of the current or observable.

# Returns
- tur_bound::Float64: TUR inequality value.
"""
function calculate_tur(entropy_production_rate, variance)
    tur_bound = variance / entropy_production_rate
    return tur_bound
end

"""
    analyze_tur(currents, entropy_rates)

Analyzes TUR for multiple observables in a system.

# Arguments
- currents::Vector{Float64}: List of steady-state currents.
- entropy_rates::Vector{Float64}: Corresponding entropy production rates.

# Returns
- results::Vector{Float64}: TUR bounds for each observable.
"""
function analyze_tur(currents, entropy_rates)
    variances = var(currents, dims=1)
    results = [calculate_tur(entropy_rates[i], variances[i]) for i in 1:length(currents)]
    return results
end
