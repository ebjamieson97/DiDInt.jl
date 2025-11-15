module DiDInt

using CategoricalArrays
using Dates
using DataFrames
using Distributions
using FixedEffectModels
using GLM
using LinearAlgebra
using Random
using Statistics

include("data_checks.jl")
include("helpers.jl")
include("estimate_didint.jl")
include("didint_plot.jl")

export didint, didint_plot

end