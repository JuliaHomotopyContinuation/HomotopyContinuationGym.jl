struct Cyclooctane <: TestSystem
    F::Vector{DP.Polynomial{true, Float64}}
    parameters::Vector{DP.PolyVar{true}}
    start_solutions::Vector{Vector{ComplexF64}}
    start_parameters::Vector{ComplexF64}
end

function Base.show(io::IO, F::Cyclooctane)
    println(io, "Cyclooctane:")
    println(io, " ◯ solutions → 1408")
    println(io, " ◯ polynomials → 17")
    println(io, " ◯ variables → 17")
    println(io, " ◯ parameters → 36")
end

system(F::Cyclooctane) = F.F
parameters(F::Cyclooctane) = F.parameters
start_solutions(F::Cyclooctane) = F.start_solutions
start_parameters(F::Cyclooctane) = F.start_parameters
nvariables(::Cyclooctane) = 17
npolynomials(::Cyclooctane) = 17
nparameters(::Cyclooctane) = 36

function Cyclooctane()
    c² = 2
    DP.@polyvar z[1:3, 1:6]
    z_vec = vec(z)[1:17] # the 17 variables in a vector
    Z = [zeros(3) z[:,1:5] [z[1,6]; z[2,6]; 0] [√c²; 0; 0]] # the eight points in a matrix

    # define the functions for cyclooctane:
    F1 = [(Z[:, i] - Z[:, i+1]) ⋅ (Z[:, i] - Z[:, i+1]) - c² for i in 1:7]
    F2 = [(Z[:, i] - Z[:, i+2]) ⋅ (Z[:, i] - Z[:, i+2]) - 8c²/3 for i in 1:6]
    F3 = (Z[:, 7] - Z[:, 1]) ⋅ (Z[:, 7] - Z[:, 1]) - 8c²/3
    F4 = (Z[:, 8] - Z[:, 2]) ⋅ (Z[:, 8] - Z[:, 2]) - 8c²/3
    f = [F1; F2; F3; F4]

    n = 2 # dimension of the cyclooctane variety
    N = 17 # ambient dimension
    DP.@polyvar Aᵥ[1:n, 1:N] bᵥ[1:n] # variables for the linear equations
    p = [vec(Aᵥ); bᵥ] # parameters
    F = [f; Aᵥ * z_vec - bᵥ] # the polynomial system we have to solve

    S₀, p₀ = load_start_pairs(@__DIR__)

    Cyclooctane(F, p, S₀, p₀)
end
