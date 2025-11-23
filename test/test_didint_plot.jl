# Tests for didint_plot()

# Ensure that didint_plot() works with good data
@testset "Basic run - parallel trends" begin
    test_plot_parallel = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                            treatment_times = TREATED_TIMES, date_format = "yyyy",
                                            treated_states = TREATED_STATES,
                                            event = false, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_parallel)
end

@testset "Basic run - event study" begin
    test_plot_event = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                            treatment_times = TREATED_TIMES, date_format = "yyyy",
                                            treated_states = TREATED_STATES,
                                            event = true, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_event)
end

# Robustness check for data with a bunch of missing values
@testset "Missing run - parallel trends" begin
    test_plot_parallel = DiDInt.didint_plot("coll", "state", "year", TEST_DATA_MISSING,
                                            treatment_times = TREATED_TIMES, date_format = "yyyy",
                                            treated_states = TREATED_STATES,
                                            event = false, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_parallel)
end

@testset "Missing run - event study" begin
    test_plot_event = DiDInt.didint_plot("coll", "state", "year", TEST_DATA_MISSING,
                                            treatment_times = TREATED_TIMES, date_format = "yyyy",
                                            treated_states = TREATED_STATES,
                                            event = true, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_event)
end