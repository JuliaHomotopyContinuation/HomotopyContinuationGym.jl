module HomotopyContinuationGym

using DelimitedFiles, LinearAlgebra, PrettyTables

import DynamicPolynomials
const DP = DynamicPolynomials

include("systems.jl")

struct PathInfo
    # per step info
    s::Vector{Float64}
    Δs::Vector{Float64}
    ω::Vector{Float64}
    accepted_rejected::Vector{Bool}
    norm_x::Vector{Float64}
    cond::Vector{Float64}
    accuracy::Vector{Float64}
    residual::Vector{Float64}
    eval_err::Vector{Float64}
    # total info
    n_factorizations::Int
    n_ldivs::Int
end

function path_info(tracker::CoreTracker, x₀, t₁, t₀)
    state = tracker.state

    s = Float64[]
    Δs = Float64[]
    ω = Float64[]
    accepted_rejected = Bool[]
    norm_x = Float64[]
    cond = Float64[]
    accuracy = Float64[]
    residual = Float64[]
    eval_err = Float64[]

    HC.init!(tracker, x₀, t₁, t₀)
    for _ in tracker
        push!(s, state.s)
        push!(Δs, state.Δs)
        push!(ω, state.ω)
        push!(accepted_rejected, !state.last_step_failed)
        push!(norm_x, maximum(abs, state.x))
        push!(cond, state.jacobian.cond[])
        push!(accuracy, state.accuracy)
        push!(residual, maximum(state.residual))
        push!(eval_err, maximum(state.eval_error))
    end

    PathInfo(
        s,
        Δs,
        ω,
        accepted_rejected,
        norm_x,
        cond,
        accuracy,
        residual,
        eval_err,
        state.jacobian.factorizations[],
        state.jacobian.ldivs[],
    )
end

path_table(info::PathInfo) = path_table(stdout, info)
function path_table(io::IO, info::PathInfo)
    header = ["", "s", "Δs", "ω", "acc", "κ", "ψ", "||x||", "||r||"]
    h1 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✗, crayon = crayon"red")
    h2 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✓, crayon = crayon"green")
    ✓✗ = map(v -> v ? :✓ : :✗, info.accepted_rejected)
    data = hcat(
        ✓✗,
        round.(info.s, sigdigits = 3),
        round.(info.Δs, sigdigits = 3),
        round.(info.ω, sigdigits = 3),
        round.(info.accuracy, sigdigits = 3),
        round.(info.cond, sigdigits = 3),
        round.(info.eval_err, sigdigits = 3),
        round.(info.norm_x, sigdigits = 3),
        round.(info.residual, sigdigits = 3),
    )
    pretty_table(io, data, header, crop = :none, highlighters = (h1, h2))
end

function Base.show(io::IO, info::PathInfo)
    println(io, "PathInfo:")
    println(
        io,
        " • # steps (✓/✗) → ",
        length(info.s),
        " ( ",
        count(info.accepted_rejected),
        " / ",
        count(!, info.accepted_rejected),
        " )",
    )
    println(io, " • # factorizations → ", info.n_factorizations)
    println(io, " • # ldivs → ", info.n_ldivs)
    path_table(io, info)
end


end
