using HomotopyContinuation, PolynomialTestSystems, BenchmarkTools
using Statistics
const HC = HomotopyContinuation

const SEED = 523432

struct BenchResult
    system::String
    time::Float64
    npaths::Int
    nsol::Int
    median_accepted::Float64
    mean_accepted::Float64
    median_rejected::Float64
    mean_rejected::Float64
    mean_finite_accepted::Float64
    mean_finite_rejected::Float64
end

Base.show(io::IO, x::BenchResult) = HC.print_fieldnames(io, x)
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::BenchResult) = x

function default_solve(f)
    solve(f; save_all_paths=true, show_progress=false, seed=SEED)
end
function bench_system(f, name, solve_func=default_solve)
    result = solve_func(f)
    t = @belapsed $solve_func($f)
    npaths = length(result)
    nsols = HC.nnonsingular(result)
    mean_finite_accepted = mean(map(r -> r.accepted_steps,
                                    results(result; onlynonsingular=true)))
    mean_finite_rejected = mean(map(r -> r.rejected_steps,
                                    results(result; onlynonsingular=true)))
    median_accepted = median(map(r -> r.accepted_steps, result))
    mean_accepted = mean(map(r -> r.accepted_steps, result))
    median_rejected = median(map(r -> r.rejected_steps, result))
    mean_rejected = mean(map(r -> r.rejected_steps, result))
    BenchResult(name, t, npaths, nsols, median_accepted, mean_accepted,
            median_rejected, mean_rejected, mean_finite_accepted, mean_finite_rejected)
end


r1 = bench_system(equations(cyclic(5)), "cyclic5")
r2 = bench_system(equations(cyclic(6)), "cyclic6")


function affine_solve(f)
    solve(f; save_all_paths=true, show_progress=false, seed=SEED, affine_tracking=true, accuracy=1e-5)
end
r3 = bench_system(equations(cyclic(6)), "cyclic6", affine_solve)
