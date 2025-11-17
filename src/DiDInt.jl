module DiDInt

using CategoricalArrays
using Dates
using DataFrames
using Distributions
using FixedEffectModels
using LinearAlgebra
using Random
using Statistics

include("data_checks.jl")
include("date_helpers.jl")
include("common_computation.jl")
include("estimate_didint.jl")
include("didint_plot.jl")

export didint, didint_plot

end