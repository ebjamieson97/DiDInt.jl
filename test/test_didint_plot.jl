# Tests for didint_plot()

# Ensure that didint_plot() works with good data
@testset "Basic run - parallel trends" begin
    # With the treated states: conditional rows (lambda) and the integrated rows of ccc = "int" (delta)
    test_plot_parallel = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                            treatment_times = TREATED_TIMES,
                                            treated_states = TREATED_STATES,
                                            event = false, covariates = ["asian", "black", "male"])
    @test any(!ismissing, test_plot_parallel.lambda) && any(!ismissing, test_plot_parallel.delta)
    # Without them (only the conditional data are asked for)
    test_plot_conditional = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                               treatment_times = TREATED_TIMES, ccc = ["hom", "state"],
                                               event = false, covariates = ["asian", "black", "male"])
    @test any(!ismissing, test_plot_conditional.lambda) && all(ismissing, test_plot_conditional.delta)
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
    test_plot_parallel = DiDInt.didint_plot("coll", "state", "year", TEST_DATA_MISSING,
                                            treatment_times = TREATED_TIMES,
                                            treated_states = TREATED_STATES,
                                            event = false, covariates = ["asian", "black", "male"])
    @test any(!ismissing, test_plot_parallel.lambda) && any(!ismissing, test_plot_parallel.delta)
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
