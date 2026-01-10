# Tests for didint()

# Try different aggregation methods for common adoption and different ccc
@testset "Agg check common - none (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234,
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check common - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234,
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399)
    @test !isnothing(result)
end

# Try different aggregation methods for staggered adoption
@testset "Agg check staggered - cohort (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "state",
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - simple" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "hom",
                           covariates = [:male, :asian, :black],
                           agg = "simple", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - sgt" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234,  ccc = "add",
                           covariates = [:male, :asian, :black],
                           agg = "sgt", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "time",
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - time" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "time", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - none" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399)
    @test !isnothing(result)
end

# Check that gvar works for staggered adoption
@testset "Gvar check - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_GVAR,
                           gvar = :gvar, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black], nperm = 399)
    @test !isnothing(result)
end

# Test that string time column works, also check start_date, end_date, and freq arguments
@testset "String time column, start_date, end_date, freq - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_STR,
                           treatment_times = string.(TREATED_TIMES),
                           date_format = "yyyy", start_date = "1989", end_date = "2000", freq = "yearly",
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black], nperm = 399)
    @test !isnothing(result)
end

# Test usage of notyet and different weighting option
@testset "notyet and weighting = none - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, notyet = true, weighting = "none",
                           covariates = [:male, :asian, :black], nperm = 399)
    @test !isnothing(result)
end

# Test that common or staggered work with missing data
@testset "Missing data - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399)
    @test !isnothing(result)
end

@testset "Missing data - common" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399)
    @test !isnothing(result)
end

# Test that the truejack option works 
@testset "truejack - common none" begin 
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399, truejack = true)
    @test !isnothing(result)
end 

@testset "truejack - common state" begin 
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399, truejack = true)
    @test !isnothing(result)
end 

@testset "truejack check staggered - cohort (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "state",
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399, truejack = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - simple" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "hom",
                           covariates = [:male, :asian, :black],
                           agg = "simple", nperm = 399, truejack = true,
                           notyet = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - sgt" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234,  ccc = "add",
                           covariates = [:male, :asian, :black],
                           agg = "sgt", nperm = 399, truejack = true,
                           notyet = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "time",
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399, truejack = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - time" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "time", nperm = 399, truejack = true)
    @test !isnothing(result)
end


@testset "truejack and notyet check staggered - time" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "time", nperm = 399, truejack = true, notyet = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - none" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399, truejack = true)
    @test !isnothing(result)
end