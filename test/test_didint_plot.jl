# Tests for didint_plot()

# Ensure that didint_plot() works with good data
@testset "Basic run - parallel trends" begin
    # Integrated parallel trends (the default) and conditional parallel trends
    test_plot_integrated = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                              treatment_times = TREATED_TIMES,
                                              treated_states = TREATED_STATES,
                                              event = false, covariates = ["asian", "black", "male"])
    @test "delta" in names(test_plot_integrated) && nrow(test_plot_integrated) > 0
    test_plot_parallel = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                            treatment_times = TREATED_TIMES, integrated = false,
                                            event = false, covariates = ["asian", "black", "male"])
    @test "lambda" in names(test_plot_parallel) && nrow(test_plot_parallel) > 0
end

@testset "Basic run - event study" begin
    test_plot_event = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                            treatment_times = TREATED_TIMES, 
                                            treated_states = TREATED_STATES,
                                            event = true, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_event)
end

# Robustness check for data with a bunch of missing values
@testset "Missing run - parallel trends" begin
    test_plot_integrated = DiDInt.didint_plot("coll", "state", "year", TEST_DATA_MISSING,
                                              treatment_times = TREATED_TIMES,
                                              treated_states = TREATED_STATES,
                                              event = false, covariates = ["asian", "black", "male"])
    @test "delta" in names(test_plot_integrated) && nrow(test_plot_integrated) > 0
    test_plot_parallel = DiDInt.didint_plot("coll", "state", "year", TEST_DATA_MISSING,
                                            treatment_times = TREATED_TIMES, integrated = false,
                                            event = false, covariates = ["asian", "black", "male"])
    @test "lambda" in names(test_plot_parallel) && nrow(test_plot_parallel) > 0
end

@testset "Missing run - event study" begin
    test_plot_event = DiDInt.didint_plot("coll", "state", "year", TEST_DATA_MISSING,
                                            treatment_times = TREATED_TIMES, 
                                            treated_states = TREATED_STATES,
                                            event = true, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_event)
end
# Event study rows are sorted by ccc option and by period relative to treatment
@testset "Event study rows are sorted" begin
    plot_event = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    event = true, covariates = ["asian", "black", "male"], ccc = ["int", "hom"])
    @test issorted(collect(zip(plot_event.ccc, plot_event.time_since_treatment)))
    @test unique(plot_event.ccc) == ["hom", "int"]
end
