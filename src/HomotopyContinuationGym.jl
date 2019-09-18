module HomotopyContinuationGym

export PathInfo, path_info

using DelimitedFiles, LinearAlgebra, PrettyTables, Random

import DynamicPolynomials
const DP = DynamicPolynomials
import HomotopyContinuation
const HC = HomotopyContinuation

include("systems.jl")

struct PathInfo
    # per step info
    s::Vector{Float64}
    Δs::Vector{Float64}
    ω::Vector{Float64}
    Δx₀::Vector{Float64}
    accepted_rejected::Vector{Bool}
    norm_x::Vector{Float64}
    cond::Vector{Float64}
    accuracy::Vector{Float64}
    residual::Vector{Float64}
    limit_accuracy::Vector{Float64}
    eval_err::Vector{Float64}
    # total info
    return_code::HC.CoreTrackerStatus
    n_factorizations::Int
    n_ldivs::Int
end

function path_info(tracker::HC.CoreTracker, x₀, t₁, t₀)
    state = tracker.state

    s = Float64[]
    Δs = Float64[]
    ω = Float64[]
    Δx₀ = Float64[]
    accepted_rejected = Bool[]
    norm_x = Float64[]
    cond = Float64[]
    accuracy = Float64[]
    limit_accuracy = Float64[]
    residual = Float64[]
    eval_err = Float64[]

    HC.init!(tracker, x₀, t₁, t₀)
    push!(s, state.s)
    push!(Δs, state.Δs)
    push!(ω, state.ω)
    push!(Δx₀, state.norm_Δx₀)
    push!(norm_x, maximum(abs, state.x))
    push!(cond, state.jacobian.cond[])
    push!(accuracy, state.accuracy)
    push!(limit_accuracy, state.limit_accuracy)
    push!(residual, maximum(state.residual))
    push!(eval_err, maximum(state.eval_error))

    first = true
    for _ in tracker
        push!(accepted_rejected, !state.last_step_failed)
        push!(s, state.s)
        push!(Δs, state.Δs)
        push!(ω, state.ω)
        push!(Δx₀, state.norm_Δx₀)
        push!(norm_x, maximum(abs, state.x))
        push!(cond, state.jacobian.cond[])
        push!(accuracy, state.accuracy)
        push!(limit_accuracy, state.limit_accuracy)
        push!(residual, maximum(state.residual))
        push!(eval_err, maximum(state.eval_error))
        first = false
    end
    push!(accepted_rejected, !state.last_step_failed)

    PathInfo(
        s,
        Δs,
        ω,
        Δx₀,
        accepted_rejected,
        norm_x,
        cond,
        accuracy,
        residual,
        limit_accuracy,
        eval_err,
        HC.status(tracker),
        state.jacobian.factorizations[],
        state.jacobian.ldivs[],
    )
end

path_table(info::PathInfo) = path_table(stdout, info)
function path_table(io::IO, info::PathInfo)
    header = ["", "s", "Δs", "ω", "|Δx₀|", "acc", "κ", "ψ", "limit_acc", "|x|", "|r|"]
    h1 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✗, crayon = crayon"red")
    h2 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✓, crayon = crayon"green")
    ✓✗ = map(v -> v ? :✓ : :✗, info.accepted_rejected)
    data = hcat(
        ✓✗,
        round.(info.s, sigdigits = 3),
        round.(info.Δs, sigdigits = 3),
        round.(info.ω, sigdigits = 3),
        round.(info.Δx₀, sigdigits = 2),
        round.(info.accuracy, sigdigits = 2),
        round.(info.cond, sigdigits = 3),
        round.(info.eval_err, sigdigits = 2),
        round.(info.limit_accuracy, sigdigits = 2),
        round.(info.norm_x, sigdigits = 3),
        round.(info.residual, sigdigits = 2),
    )
    pretty_table(io, data, header, crop = :none, highlighters = (h1, h2))
end

function Base.show(io::IO, info::PathInfo)
    println(io, "PathInfo:")
    println(io, " • # return code → ", info.return_code)
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
