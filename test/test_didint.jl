# Tests for didint()

@testset "Basic run - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           date_format = "yyyy", seed = 1234, 
                           covariates = [:male, :asian, :black])
    @test !isnothing(result)
end