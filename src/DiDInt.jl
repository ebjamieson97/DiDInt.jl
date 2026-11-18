module DiDInt

using CategoricalArrays
using Dates
using DataFrames
using Distributions
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using StatsAPI
using Tables

import StatsAPI: coef, coefnames, stderror

include("data_checks.jl")
include("date_helpers.jl")
include("common_computation.jl")
include("estimate_didint.jl")
include("batched_engine.jl")
include("results.jl")
include("plots.jl")
include("didint_plot.jl")

export didint, didint_plot, DiDIntResult, DiDIntPlot, coef, coefnames, stderror

end
