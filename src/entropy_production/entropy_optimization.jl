"""
    compute_entropy(J::Array{Float64, 2}, X::Array{Float64, 2}) -> Float64

Computes the entropy production for given fluxes J and forces X.
"""
function compute_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})
    return sum(J .* X)
end

"""
    entropy_gradient(J::Array{Float64, 2}, X::Array{Float64, 2}) -> Tuple{Array{Float64, 2}, Array{Float64, 2}}

Calculates the gradient of the entropy function with respect to J and X.
"""
function entropy_gradient(J::Array{Float64, 2}, X::Array{Float64, 2})
    ∇J = X  # Derivative of sum(J .* X) with respect to J is X
    ∇X = J  # Derivative of sum(J .* X) with respect to X is J
    return (∇J, ∇X)
end

"""
    gradient_descent(J::Array{Float64, 2}, X::Array{Float64, 2}, learning_rate::Float64, max_iterations::Int) -> Tuple

Performs gradient descent to minimize the entropy production.
- learning_rate: The step size at each iteration.
- max_iterations: Maximum number of iterations to perform.
Returns the optimized fluxes and forces.
"""
function gradient_descent(J::Array{Float64, 2}, X::Array{Float64, 2}, learning_rate::Float64, max_iterations::Int)
    for i in 1:max_iterations
        ∇J, ∇X = entropy_gradient(J, X)
        J -= learning_rate * ∇J
        X -= learning_rate * ∇X
        
        # Optionally print the entropy to monitor convergence
        @show compute_entropy(J, X)
    end
    return (J, X)
end
