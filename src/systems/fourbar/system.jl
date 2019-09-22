export FourBar

struct FourBar <: TestSystem
    F::Vector{DP.Polynomial{true,Int}}
    parameters::Vector{DP.PolyVar{true}}
    start_solutions::Vector{Vector{ComplexF64}}
    start_parameters::Vector{ComplexF64}
end

function Base.show(io::IO, F::FourBar)
    println(io, "FourBar:")
    println(io, " ◯ solutions → 1442")
    println(io, " ◯ polynomials → 24")
    println(io, " ◯ variables → 24")
    println(io, " ◯ parameters → 16")
end
system(F::FourBar) = F.F
parameters(F::FourBar) = F.parameters
start_solutions(F::FourBar) = F.start_solutions
start_parameters(F::FourBar) = F.start_parameters
nvariables(::FourBar) = 24
npolynomials(::FourBar) = 24
nparameters(::FourBar) = 16

function FourBar()
    DP.@polyvar x a y b x̂ â ŷ b̂
    DP.@polyvar γ[1:8] γ̂[1:8] δ[1:8] δ̂[1:8]
    #system of polynomials
    D1 = [(â * x - δ̂[i] * x) * γ[i] + (a * x̂ - δ[i] * x̂) * γ̂[i] + (â - x̂) * δ[i] +
          (a - x) * δ̂[i] - δ[i] * δ̂[i] for i = 1:8]
    D2 = [(b̂ * y - δ̂[i] * y) * γ[i] + (b * ŷ - δ[i] * ŷ) * γ̂[i] + (b̂ - ŷ) * δ[i] +
          (b - y) * δ̂[i] - δ[i] * δ̂[i] for i = 1:8]
    D3 = [γ[i] * γ̂[i] + γ[i] + γ̂[i] for i = 1:8]
    F = [D1; D2; D3]
    params = [δ; δ̂]

    S, p = load_start_pairs(@__DIR__)

    FourBar(F, params, S, p)
end
