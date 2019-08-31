struct Steiner <: TestSystem
    F::Vector{DP.Polynomial{true, Int}}
    parameters::Vector{DP.PolyVar{true}}
    start_solutions::Vector{Vector{ComplexF64}}
    start_parameters::Vector{ComplexF64}
end

function Base.show(io::IO, F::Steiner)
    println(io, "Steiner:")
    println(io, " ◯ solutions → 3264")
    println(io, " ◯ polynomials → 15")
    println(io, " ◯ variables → 15")
    println(io, " ◯ parameters → 30")
end

start_solutions(F::Steiner) = F.start_solutions
start_parameters(F::Steiner) = F.start_parameters
nvariables(::Steiner) = 15
npolynomials(::Steiner) = 15
nparameters(::Steiner) = 30

function Steiner()
    DP.@polyvar x[1:2] a[1:5] c[1:6] y[1:2, 1:5]

    #tangential conics
    f = a[1]*x[1]^2 + a[2]*x[1]*x[2] + a[3]*x[2]^2 + a[4]*x[1]  + a[5]*x[2] + 1;
    ∇ = DP.differentiate(f, x)
    #5 conics
    g = c[1]*x[1]^2 + c[2]*x[1]*x[2] + c[3]*x[2]^2 + c[4]*x[1]  + c[5]*x[2] + c[6];
    ∇_2 = DP.differentiate(g, x)
    #the general system
    #f_a_0 is tangent to g_b₀ at x₀
    function Incidence(f,a₀,g,b₀,x₀)
        fᵢ = f(x=>x₀, a=>a₀)
        ∇ᵢ = [∇ᵢ(x=>x₀, a=>a₀) for ∇ᵢ in ∇]
        Cᵢ = g(x=>x₀, c=>b₀)
        ∇_Cᵢ = [∇ⱼ(x=>x₀, c=>b₀) for ∇ⱼ in ∇_2]

        [fᵢ; Cᵢ;  det([∇ᵢ ∇_Cᵢ])]
    end
    DP.@polyvar v[1:6, 1:5]
    F = vcat(map(i -> Incidence(f,a,g, v[:,i], y[:,i]), 1:5)...)

    start_solutions_real = readdlm(joinpath(@__DIR__, "solutions_real.txt"), '\t', Float64, '\n')
    start_solutions_imag = readdlm(joinpath(@__DIR__, "solutions_imag.txt"), '\t', Float64, '\n')
    start_solutions = [[complex(start_solutions_real[i, j], start_solutions_imag[i,j]) for j=1:15] for i=1:3264]
    start_parameters = readdlm(joinpath(@__DIR__, "parameters.txt"), '\t', ComplexF64)[:,1]

    Steiner(F, vec(v), start_solutions, start_parameters)
end
