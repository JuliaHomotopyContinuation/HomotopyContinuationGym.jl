struct Bacillus <: TestSystem
    F::Vector{DP.Polynomial{true, Float64}}
    parameters::Vector{DP.PolyVar{true}}
    start_solutions::Vector{Vector{ComplexF64}}
    start_parameters::Vector{ComplexF64}
end

function Base.show(io::IO, F::Bacillus)
    println(io, "Bacillus:")
    println(io, " ◯ solutions → 44")
    println(io, " ◯ polynomials → 10")
    println(io, " ◯ variables → 10")
    println(io, " ◯ parameters → 23")
end
system(F::Bacillus) = F.F
parameters(F::Bacillus) = F.parameters
start_solutions(F::Bacillus) = F.start_solutions
start_parameters(F::Bacillus) = F.start_parameters
nvariables(::Bacillus) = 10
npolynomials(::Bacillus) = 10
nparameters(::Bacillus) = 23

function Bacillus()
    DP.@polyvar x[1:10] p[1:23]

    f = [
        (-1 * p[17] * x[1] + -2 * p[1] * (x[1] ^ 2 / 2) + 2 * p[2] * x[2])*(p[20] + x[7]) + p[21] * p[18] * ((1 + p[19] * x[7])),
        -1 * p[17] * x[2] + p[1] * (x[1] ^ 2 / 2) + -1 * p[2] * x[2] + -1 * p[4] * x[2] * x[4] + p[9] * x[3] + p[14] * x[3] + -1 * p[6] * x[2] * x[7] + p[11] * x[8],
        -1 * p[17] * x[3] + p[4] * x[2] * x[4] + -1 * p[9] * x[3] + -1 * p[5] * x[3] * x[4] + p[10] * x[5] + -1 * p[14] * x[3] + p[15] * x[5] + p[7] * x[8] * x[4] + -1 * p[12] * x[3] * x[7],
        (-1 * p[17] * x[4] + -1 * p[4] * x[2] * x[4] + p[9] * x[3] + -1 * p[5] * x[3] * x[4] + p[10] * x[5] + -1 * p[7] * x[8] * x[4] + p[12] * x[3] * x[7] + p[16] * x[9])*(p[20] + x[7]) + p[22] * p[18] * (1 + p[19] * x[7]),
        -1 * p[17] * x[5] + p[5] * x[3] * x[4] + -1 * p[10] * x[5] + -1 * p[15] * x[5],
        -1 * p[17] * x[6] + p[14] * x[3] + p[15] * x[5] + -1 * p[8] * x[6] * x[10] + p[13] * x[9],
        (-1 * p[17] * x[7] + -1 * p[6] * x[2] * x[7] + p[11] * x[8] + p[7] * x[8] * x[4] + -1 * p[12] * x[3] * x[7])*(p[20] + x[7]) + p[18] * (1 + p[19] * x[7]),
        -1 * p[17] * x[8] + p[6] * x[2] * x[7] + -1 * p[11] * x[8] + -1 * p[7] * x[8] * x[4] + p[12] * x[3] * x[7],
        -1 * p[17] * x[9] + p[8] * x[6] * x[10] + -1 * p[13] * x[9] + -1 * p[16] * x[9],
        x[10] + x[9] - p[23]]

    S, p₀ = load_start_pairs(@__DIR__)

    Bacillus(f, p, S, p₀)
end
