using DiDInt
using Test
using DataFrames
using CSV  
using Random
using Logging
Random.seed!(1234)

# Load merit data for testing and build treated_states and treated_times
test_data = CSV.read(joinpath(@__DIR__, "data", "merit.csv"), DataFrame)
test_data.state = string.(test_data.state)
treated_states = ["34", "57", "58", "59", "61", "64", "71", "72", "85", "88"]
control_states = ["11", "12", "13", "14", "15", "16", "21", "22", "23", "31", 
                  "32", "33", "35", "41", "42", "43", "44", "45", "46", "47",
                  "51", "52", "53", "54", "55", "56", "62", "63", "73", "74",
                  "81", "82", "83", "84", "86", "87", "91", "92", "93", "94", "95"]
const CONTROL_STATES = control_states
selected_controls = shuffle(control_states)[1:10]
selected_states = vcat(selected_controls, treated_states)
const TEST_DATA_FULL = test_data
TEST_DATA_FULL[(TEST_DATA_FULL.state .== "95") .& (TEST_DATA_FULL.year .== 2000), :asian] .= 1
TEST_DATA_FULL[(TEST_DATA_FULL.state .== "95") .& (TEST_DATA_FULL.year .== 2000), :male] .= 1
TEST_DATA_FULL[(TEST_DATA_FULL.state .== "95") .& (TEST_DATA_FULL.year .== 2000), :black] .= 1
test_data_selected = filter(row -> row.state ∈ selected_states, test_data)
const TEST_DATA = test_data_selected
const TREATED_STATES = unique(TEST_DATA[TEST_DATA.merit .== 1, "state"])
treated_times = Vector{Number}(undef, length(TREATED_STATES))
for i in eachindex(TREATED_STATES)
    year_values = TEST_DATA[(string.(TEST_DATA.state) .== TREATED_STATES[i]) .& (TEST_DATA.merit .== 1), "year"]
    treated_times[i] = minimum(year_values)
end
treated_times_time_collinear = treated_times
treated_times_time_collinear[1] = 1990
const TREATED_TIMES = treated_times
const TREATED_TIMES_T_COLLINEAR = treated_times_time_collinear
# Also build test data that has a bunch of missing entries for robustness checks
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

# Create test datasets that trigger collinearity for each ccc case
# Use first treated state and its times — guaranteed to be in the analysis
target_state = TREATED_STATES[1]
target_times = sort(unique(TEST_DATA[TEST_DATA.state .== target_state, :year]))
target_time = target_times[1]

# ccc = "state": male constant across ALL time periods for one treated state
# issubset condition will trigger since all times for target_state will be missing
const TEST_DATA_STATE_COLLINEAR = copy(TEST_DATA)
TEST_DATA_STATE_COLLINEAR[(TEST_DATA_STATE_COLLINEAR.state .== target_state), :male] .= 1
TEST_DATA_STATE_COLLINEAR[(TEST_DATA_STATE_COLLINEAR.state .== target_state), :asian] .= 1
TEST_DATA_STATE_COLLINEAR[(TEST_DATA_STATE_COLLINEAR.state .== target_state), :black] .= 1

# ccc = "time": male constant across ALL states for one time period
# issubset condition will trigger since all states for target_time will be missing  
const TEST_DATA_TIME_COLLINEAR = copy(TEST_DATA)
TEST_DATA_TIME_COLLINEAR[TEST_DATA_TIME_COLLINEAR.year .== target_time, :male] .= 1

# ccc = "hom": male constant across entire dataset
const TEST_DATA_HOM_COLLINEAR = copy(TEST_DATA)
TEST_DATA_HOM_COLLINEAR[:, :male] .= 1

# Run tests
@testset "DiDInt.jl" begin
    @testset "didint()" begin
        include("test_didint.jl")
    end
    @testset "didint_plot()" begin
        include("test_didint_plot.jl")
    end
end