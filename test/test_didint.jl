# Tests for didint()

# Try different aggregation methods for common adoption
@testset "Agg check common - none (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none")
    @test !isnothing(result)
end

@testset "Agg check common - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "state")
    @test !isnothing(result)
end

# Try different aggregation methods for staggered adoption
@testset "Agg check staggered - cohort (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "cohort")
    @test !isnothing(result)
end

@testset "Agg check staggered - simple" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "simple")
    @test !isnothing(result)
end

@testset "Agg check staggered - sgt" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "sgt")
    @test !isnothing(result)
end

@testset "Agg check staggered - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "state")
    @test !isnothing(result)
end

@testset "Agg check staggered - time" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "time")
    @test !isnothing(result)
end

@testset "Agg check staggered - none" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none")
    @test !isnothing(result)
end

# Check that gvar works for staggered adoption
@testset "Gvar check - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_GVAR,
                           gvar = :gvar, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black])
    @test !isnothing(result)
end

# Test that string time column works, also check start_date, end_date, and freq arguments
@testset "String time column, start_date, end_date, freq - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_STR,
                           treatment_times = string.(TREATED_TIMES),
                           date_format = "yyyy", start_date = "1989", end_date = "2000", freq = "yearly",
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black])
    @test !isnothing(result)
end

# Test that common or staggered work with missing data
@testset "Missing data - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "cohort")
    @test !isnothing(result)
end

@testset "Missing data - common" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "cohort")
    @test !isnothing(result)
end