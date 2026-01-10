using DiDInt
using Test
using DataFrames
using CSV  
using Random

Random.seed!(1234)

# Load merit data for testing and build treated_states and treated_times
test_data = CSV.read(joinpath(@__DIR__, "data", "merit.csv"), DataFrame)
test_data.state = string.(test_data.state)
treated_states = ["34", "57", "58", "59", "61", "64", "71", "72", "85", "88"]
control_states = ["11", "12", "13", "14", "15", "16", "21", "22", "23", "31", 
                  "32", "33", "35", "41", "42", "43", "44", "45", "46", "47",
                  "51", "52", "53", "54", "55", "56", "62", "63", "73", "74",
                  "81", "82", "83", "84", "86", "87", "91", "92", "93", "94", "95"]
selected_controls = shuffle(control_states)[1:10]
selected_states = vcat(selected_controls, treated_states)
test_data = filter(row -> row.state ∈ selected_states, test_data)
const TEST_DATA = test_data
const TREATED_STATES = unique(TEST_DATA[TEST_DATA.merit .== 1, "state"])
treated_times = Vector{Number}(undef, length(TREATED_STATES))
for i in eachindex(TREATED_STATES)
    year_values = TEST_DATA[(string.(TEST_DATA.state) .== TREATED_STATES[i]) .& (TEST_DATA.merit .== 1), "year"]
    treated_times[i] = minimum(year_values)
end
const TREATED_TIMES = treated_times

# Also build test data that has a bunch of mising entries for robustness checks
test_data_missing = allowmissing(copy(TEST_DATA))
cols_to_randomize = ["asian", "black", "male", "merit", "coll", "year", "state"]
for col in cols_to_randomize
    test_data_missing[:, col] = ifelse.(rand(nrow(test_data_missing)) .< 0.15, missing, test_data_missing[:, col])
end
const TEST_DATA_MISSING = test_data_missing

# Make a copy of TEST_DATA that uses a gvar column 
const TEST_DATA_GVAR = transform(groupby(TEST_DATA, :state), [:year, :merit] => ((y, m) -> isempty(y[m .== 1]) ? missing : minimum(y[m .== 1])) => :gvar)

# Make a copy of TEST_DATA that has string columns for time 
const TEST_DATA_STR = transform(TEST_DATA, :year => ByRow(string) => :year)

# Run tests
@testset "DiDInt.jl" begin

    @testset "didint()" begin
        include("test_didint.jl")
    end

    @testset "didint_plot()" begin
        include("test_didint_plot.jl")
    end

end