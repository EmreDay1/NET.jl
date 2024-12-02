using LinearAlgebra
using ForwardDiff

"""
    curvature_based_stability(g_ab)

Determine the stability of a thermodynamic system using the Ricci curvature.

# Parameters
- `g_ab`: The metric tensor of the system.

# Returns
- `"Stable"` if `det(g_ab) > 0`.
- `"Unstable"` if `det(g_ab) ≤ 0`.
"""
function curvature_based_stability(g_ab::Matrix{Float64})
    determinant = det(g_ab)  
    return determinant > 0 ? "Stable" : "Unstable"
end

"""
    flux_force_legendre(S, X, dS_dX)

Perform a Legendre transformation for flux-force duality in non-equilibrium systems.

# Parameters
- `S`: The entropy function value.
- `X`: A vector of extensive variables (forces).
- `dS_dX`: A vector of partial derivatives of entropy with respect to each variable in `X`.

# Returns
- The Legendre-transformed entropy: `Φ = S - Σ (Xᵢ * dS/dXᵢ)`.
"""
function flux_force_legendre(S::Float64, X::Vector{Float64}, dS_dX::Vector{Float64})
    Φ = S - dot(X, dS_dX) 
    return Φ
end

"""
    grand_potential(U, T, S, μ, N)

Compute the grand potential for open thermodynamic systems.

# Parameters
- `U`: Internal energy (scalar).
- `T`: Temperature (scalar).
- `S`: Entropy (scalar).
- `μ`: Chemical potential (scalar).
- `N`: Number of particles (scalar).

# Returns
- The grand potential: `Ω = U - TS - μN`.
"""
function grand_potential(U::Float64, T::Float64, S::Float64, μ::Float64, N::Float64)
    Ω = U - T * S - μ * N 
    return Ω
end

"""
    maxwell_relation(T, V, P, S, tolerance=1e-6)

Verify and compute Maxwell relations derived from thermodynamic potentials.

# Parameters
- `T`: Function for temperature as a function of volume, `T(V)`.
- `V`: Independent variable for volume.
- `P`: Function for pressure as a function of entropy, `P(S)`.
- `S`: Independent variable for entropy.
- `tolerance`: Numerical tolerance for validation (default: `1e-6`).

# Returns
- A Boolean indicating if `-∂P/∂S ≈ ∂T/∂V` holds true within the tolerance.
- The computed numerical values of `-∂P/∂S` and `∂T/∂V` for verification.
"""
function maxwell_relation(T::Function, V::Float64, P::Function, S::Float64, tolerance=1e-6)
    dT_dV = ForwardDiff.derivative(T, V)
    dP_dS = ForwardDiff.derivative(P, S) 
    maxwell_check = abs(-dP_dS - dT_dV) <= tolerance
    return maxwell_check ? "Maxwell relation holds" : "Maxwell relation violated: -∂P/∂S = $(-dP_dS), ∂T/∂V = $dT_dV"
end

"""
    extended_thermodynamic_potential(E, S, X, F, σ, J)

Compute the extended thermodynamic potential including entropy production and heat flux.

# Parameters
- `E`: Internal energy (scalar).
- `S`: Entropy (scalar).
- `X`: Extensive variables (vector).
- `F`: Intensive variables (vector).
- `σ`: Entropy production rate (scalar).
- `J`: Heat flux (scalar).

# Returns
- The extended potential `Φ_ext = E - Σ(Xᵢ * Fᵢ) - S + σ + J`.
"""
function extended_thermodynamic_potential(E::Float64, S::Float64, X::Vector{Float64}, F::Vector{Float64}, σ::Float64, J::Float64)
    Φ_ext = E - dot(X, F) - S + σ + J
    return Φ_ext
end

"""
    legendre_transform(potential, variables, derivatives)

Perform a Legendre transformation for the given potential.

# Parameters
- `potential`: The initial potential function (scalar).
- `variables`: Vector of extensive variables (vector).
- `derivatives`: Vector of partial derivatives of the potential with respect to variables.

# Returns
- The transformed potential after the Legendre transformation.
"""
function legendre_transform(potential::Float64, variables::Vector{Float64}, derivatives::Vector{Float64})
    Φ_trans = potential - dot(variables, derivatives)
    return Φ_trans
end

"""
    contact_form_extended(E, X, I, σ, J)

Compute the extended contact form with dissipation effects.

# Parameters
- `E`: Energy function.
- `X`: Extensive variables.
- `I`: Intensive variables.
- `σ`: Entropy production.
- `J`: Heat flux.

# Returns
- Extended contact form: Θ_ext = dE - Σ(Iᵃ dXᵃ) + σ + J.
"""
function contact_form_extended(E::Float64, X::Vector{Float64}, I::Vector{Float64}, σ::Float64, J::Float64)
    dE = "dE"
    dX = ["dX$i" for i in 1:length(X)]
    terms = [string(I[i], "*", dX[i]) for i in 1:length(X)]
    Θ_ext = string(dE, " - ", join(terms, " - "), " + ", σ, " + ", J)
    return Θ_ext
end

"""
    ricci_curvature(metric_tensor)

Compute the Ricci curvature for the given thermodynamic metric.

# Parameters
- `metric_tensor`: The metric tensor (matrix).

# Returns
- Ricci curvature tensor.
"""
function ricci_curvature(metric_tensor::Matrix)
    return det(metric_tensor)
end
