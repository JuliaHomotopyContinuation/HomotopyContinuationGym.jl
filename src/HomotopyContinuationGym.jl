module HomotopyContinuationGym

using DelimitedFiles, LinearAlgebra, Random

import DynamicPolynomials
const DP = DynamicPolynomials
import HomotopyContinuation
const HC = HomotopyContinuation

import HomotopyContinuation: parameters, nvariables, npolynomials, nparameters, solve

include("systems.jl")

function solve(sys::TestSystem; kwargs...)
    solve(
        system(sys),
        start_solutions(sys);
        parameters = parameters(sys),
        start_parameters = start_parameters(sys),
        kwargs...,
    )
end

end
