module DiDInt

using CategoricalArrays
using Dates
using DataFrames
using Distributions
using FixedEffectModels
using LinearAlgebra
using Random
using Statistics

include("helpers.jl")
include("estimate_didint.jl")

export didint

end