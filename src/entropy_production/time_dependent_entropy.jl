using DifferentialEquations

"""
    analyze_time_dependent_entropy_production(J_func::Function, X_func::Function, t_span::Tuple{Float64, Float64}) -> ODESolution

Analyzes time-dependent entropy production based on time-varying fluxes and forces.
- J_func, X_func: Functions that return J and X matrices as a function of time.
- t_span: Tuple specifying the start and end time for the analysis.
Returns the solution of the ODE problem representing entropy production over time.
"""
function analyze_time_dependent_entropy_production(J_func::Function, X_func::Function, t_span::Tuple{Float64, Float64})
    function dσ_dt(σ, p, t)
        J = J_func(t)
        X = X_func(t)
        return sum(J .* X)  # Rate of entropy production
    end

    σ_0 = 0.0  # Initial entropy production
    prob = ODEProblem(dσ_dt, σ_0, t_span)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    return sol
end
