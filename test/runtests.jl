using DiDInt
using Test
using DataFrames
using CSV  

# Load merit data for testing and build treated_states and treated_times
test_data = CSV.read(joinpath(@__DIR__, "data", "merit.csv"), DataFrame)
test_data.state = string.(test_data.state)
const TEST_DATA = test_data
const TREATED_STATES = unique(TEST_DATA[TEST_DATA.merit .== 1, "state"])
treated_times = Vector{Number}(undef, length(TREATED_STATES))
for i in eachindex(TREATED_STATES)
    year_values = TEST_DATA[(string.(TEST_DATA.state) .== TREATED_STATES[i]) .& (TEST_DATA.merit .== 1), "year"]
    treated_times[i] = minimum(year_values)
end
const TREATED_TIMES = treated_times

@testset "DiDInt.jl" begin

    @testset "didint()" begin
        include("test_didint.jl")
    end

    @testset "didint_plot()" begin
        include("test_didint_plot.jl")
    end

end