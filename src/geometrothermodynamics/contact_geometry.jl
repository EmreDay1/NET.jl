using LinearAlgebra

"""
    contact_form(E, X, I)

Compute the contact form for a thermodynamic phase space.

# Parameters
- `E`: The thermodynamic potential, a function of extensive variables.
- `X`: A vector of extensive variables (e.g., volume, entropy).
- `I`: A vector of intensive variables (e.g., pressure, temperature).

# Returns
- A symbolic representation of the contact form: `Θ = dE - Σ Iᵃ dXᵃ`.
"""
function contact_form(E::Function, X::Vector{Symbol}, I::Vector{Symbol})
    n = length(X)
    dE = Symbol("∂$E/∂$(join(X, ','))") 
    dX = ["d$(X[i])" for i in 1:n]
    Θ = string(dE, " - ", join([string(I[i], " * ", dX[i]) for i in 1:n], " - "))
    return Θ
end

"""
    non_degeneracy_condition(contact_form, n)

Check if the manifold satisfies the non-degeneracy condition for thermodynamic systems.

# Parameters
- `contact_form`: The contact form (output from `contact_form` function).
- `n`: Dimensionality of the thermodynamic system.

# Returns
- A condition ensuring `Θ ∧ (dΘ)ⁿ ≠ 0`.
"""
function non_degeneracy_condition(contact_form::String, n::Int)
    return "$(contact_form) ∧ (d$(contact_form))^$n ≠ 0"
end

"""
    thermodynamic_metric(g_ab, variables)

Define a thermodynamic metric for a system's manifold.

# Parameters
- `g_ab`: A symmetric tensor (matrix) representing the metric components.
- `variables`: A vector of thermodynamic variables.

# Returns
- The metric: `g = gᵃᵇ dxᵃ dxᵇ`.
"""
function thermodynamic_metric(g_ab::Matrix, variables::Vector{Symbol})
    dx = ["d$var" for var in variables]
    metric = sum([g_ab[i, j] * "($dx[$i])($dx[$j])" for i in 1:length(dx), j in 1:length(dx)])
    return metric
end

"""
    entropy_based_metric(S, U, X)

Calculate the metric based on the entropy representation.

# Parameters
- `S`: The entropy function.
- `U`: Internal energy.
- `X`: A vector of extensive variables.

# Returns
- The metric incorporating entropy derivatives.
"""
function entropy_based_metric(S::Function, U::Symbol, X::Vector{Symbol})
    d2S_dU2 = Symbol("∂²$S/∂$U²")  #
    d2S_dXiXj = [Symbol("∂²$S/∂$(X[i])∂$(X[j])") for i in 1:length(X), j in 1:length(X)]
    metric_U = "-($d2S_dU2 * d$U^2)"
    metric_X = join([string(d2S_dXiXj[i], " d", X[i], " d", X[j]) for i in 1:length(X)], " + ")
    return metric_U * " + " * metric_X
end
