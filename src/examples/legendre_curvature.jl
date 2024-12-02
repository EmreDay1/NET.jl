using LinearAlgebra
using ForwardDiff
using Plots


time = 0:0.1:10
E = 100.0
S = 20.0
σ = 0.5
J = 1.0
X = [10.0, 15.0]
F = [2.0, 3.0]
g_ab = [1.0 0.5; 0.5 2.0]

Φ_ext_values = []
Φ_trans_values = []
ricci_values = []


for t in time
    E_t = E + 5 * sin(0.2 * t)
    σ_t = σ + 0.1 * cos(0.3 * t)
    J_t = J + 0.2 * sin(0.5 * t)

    Φ_ext = extended_thermodynamic_potential(E_t, S, X, F, σ_t, J_t)
    push!(Φ_ext_values, Φ_ext)


    dΦ_ext_dX = ForwardDiff.gradient(x -> extended_thermodynamic_potential(E_t, S, x, F, σ_t, J_t), X)

    
    Φ_trans = legendre_transform(Φ_ext, X, dΦ_ext_dX)
    push!(Φ_trans_values, Φ_trans)


    ricci = ricci_curvature(g_ab)
    push!(ricci_values, ricci)

   
    if t == 5.0  
        println("Time = $t")
        println("Extended Thermodynamic Potential: Φ_ext = $Φ_ext")
        println("Partial Derivatives of Φ_ext with respect to X: $dΦ_ext_dX")
        println("Legendre Transformed Potential: Φ_trans = $Φ_trans")
        println("Ricci Curvature: $ricci")
    end
end

plot(time, Φ_ext_values, label="Φ_ext (Extended Potential)", xlabel="Time", ylabel="Values", lw=2)
plot!(time, Φ_trans_values, label="Φ_trans (Legendre Transform)", lw=2)
plot!(time, ricci_values, label="Ricci Curvature", lw=2, legend=:topright)
