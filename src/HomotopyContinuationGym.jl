module HomotopyContinuationGym

export CTPathInfo, path_info

using DelimitedFiles, LinearAlgebra, PrettyTables, Random, Printf

import DynamicPolynomials
const DP = DynamicPolynomials
import HomotopyContinuation
const HC = HomotopyContinuation

include("systems.jl")

################
## CTPathInfo ##
################

struct CTPathInfo
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

    CTPathInfo(
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

path_table(info::CTPathInfo) = path_table(stdout, info)
function path_table(io::IO, info::CTPathInfo)
    header = ["", "s", "Δs", "ω", "|Δx₀|", "acc", "κ", "ψ", "limit_acc", "|x|", "|r|"]
    h1 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✗, crayon = crayon"red")
    h2 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✓, crayon = crayon"green")
    ✓✗ = map(v -> v ? :✓ : :✗, info.accepted_rejected)
    sigdigits3(x) = @sprintf("%.3g", x)
    sigdigits2(x) = @sprintf("%.2g", x)
    data = hcat(
        ✓✗,
        sigdigits3.(info.s),
        sigdigits3.(info.Δs),
        sigdigits3.(info.ω),
        sigdigits2.(info.Δx₀),
        sigdigits2.(info.accuracy),
        sigdigits3.(info.cond),
        sigdigits2.(info.eval_err),
        sigdigits2.(info.limit_accuracy),
        sigdigits3.(info.norm_x),
        sigdigits2.(info.residual),
    )
    pretty_table(io, data, header, crop = :none, highlighters = (h1, h2))
end

function Base.show(io::IO, info::CTPathInfo)
    println(io, "CTPathInfo:")
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


struct PTPathInfo
    # per step info
    s::Vector{Float64}
    Δs::Vector{Float64}
    ω::Vector{Float64}
    Δx₀::Vector{Float64}
    accepted_rejected::Vector{Bool}
    norm_x::Vector{Float64}
    cond::Vector{Float64}
    accuracy::Vector{Float64}
    limit_accuracy::Vector{Float64}
    norm_ν̇::Vector{Float64}
    min_ν::Vector{Float64}
    max_ν::Vector{Float64}
    # total info
    return_code::HC.PathTrackerStatus
    n_factorizations::Int
    n_ldivs::Int
end

function path_info(tracker::HC.PathTracker, x₀, t₀=nothing)
    ct_state = tracker.core_tracker.state

    s = Float64[]
    Δs = Float64[]
    ω = Float64[]
    Δx₀ = Float64[]
    accepted_rejected = Bool[]
    norm_x = Float64[]
    cond = Float64[]
    accuracy = Float64[]
    limit_accuracy = Float64[]
    norm_ν̇ = Float64[]
    min_ν = Float64[]
    max_ν = Float64[]

    if t₀ !== nothing
        HC.init!(tracker, x₀, 0.0, t₀)
    else
        HC.init!(tracker, x₀)
    end

    push!(s, real(HC.current_t(ct_state)))
    push!(Δs, real(HC.current_Δt(ct_state)))
    push!(ω, ct_state.ω)
    push!(Δx₀, ct_state.norm_Δx₀)
    push!(norm_x, maximum(abs, ct_state.x))
    push!(cond, ct_state.jacobian.cond[])
    push!(accuracy, ct_state.accuracy)
    push!(limit_accuracy, ct_state.limit_accuracy)
    push!(norm_ν̇, maximum(abs, tracker.state.valuation.ν̇))
    push!(max_ν, maximum(tracker.state.valuation.ν))
    push!(min_ν, minimum(tracker.state.valuation.ν))

    first = true
    for _ in tracker
        push!(accepted_rejected, !ct_state.last_step_failed)
        push!(s, real(HC.current_t(ct_state)))
        push!(Δs, real(HC.current_Δt(ct_state)))
        push!(ω, ct_state.ω)
        push!(Δx₀, ct_state.norm_Δx₀)
        push!(norm_x, maximum(abs, ct_state.x))
        push!(cond, ct_state.jacobian.cond[])
        push!(accuracy, ct_state.accuracy)
        push!(limit_accuracy, ct_state.limit_accuracy)
        push!(norm_ν̇, maximum(abs, tracker.state.valuation.ν̇))
        push!(max_ν, maximum(tracker.state.valuation.ν))
        push!(min_ν, minimum(tracker.state.valuation.ν))
        first = false
    end
    push!(accepted_rejected, !ct_state.last_step_failed)

    PTPathInfo(
        s,
        Δs,
        ω,
        Δx₀,
        accepted_rejected,
        norm_x,
        cond,
        accuracy,
        limit_accuracy,
        norm_ν̇,
        min_ν,
        max_ν,
        HC.status(tracker),
        ct_state.jacobian.factorizations[],
        ct_state.jacobian.ldivs[],
    )
end

path_table(info::PTPathInfo) = path_table(stdout, info)
function path_table(io::IO, info::PTPathInfo)
    header = [
        "",
        "s",
        "Δs",
        "|ν̇|",
        "min_ν",
        "max_ν",
        "κ",
        "limit_acc",
        "acc",
        "ω",
        "|Δx₀|",
        "|x|",
    ]
    h1 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✗, crayon = crayon"red")
    h2 = Highlighter(f = (data, i, j) -> j == 1 && data[i, 1] == :✓, crayon = crayon"green")
    ✓✗ = map(v -> v ? :✓ : :✗, info.accepted_rejected)
    sigdigits3(x) = @sprintf("%.3g", x)
    sigdigits2(x) = @sprintf("%.2g", x)
    data = hcat(
        ✓✗,
        sigdigits3.(info.s),
        sigdigits3.(info.Δs),
        sigdigits2.(info.norm_ν̇),
        sigdigits2.(info.min_ν),
        sigdigits2.(info.max_ν),
        sigdigits3.(info.cond),
        sigdigits2.(info.limit_accuracy),
        sigdigits2.(info.accuracy),
        sigdigits3.(info.ω),
        sigdigits2.(info.Δx₀),
        sigdigits3.(info.norm_x),
    )
    pretty_table(io, data, header, crop = :none, highlighters = (h1, h2))
end

function Base.show(io::IO, info::PTPathInfo)
    println(io, "PTPathInfo:")
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
