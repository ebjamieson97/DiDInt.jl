# Tests for didint_plot() on the merit data

plot_merit(data; kwargs...) = DiDInt.didint_plot("coll", "state", "year", data; treatment_times = TREATED_TIMES,
                                                 covariates = ["asian", "black", "male"], kwargs...)

# Ensure that didint_plot() works with good data
@testset "Basic run - parallel trends" begin
    # With the treated states: every kind of plot data
    x = plot_merit(TEST_DATA; treated_states = TREATED_STATES)
    @test DiDInt.available_plot_types(x) == ["event", "conditional", "conditional2", "integrated", "integrated2"]
    @test !any(ismissing, x.conditional.lambda) && nrow(x.integrated) > 0
    # Without them (only the conditional data can be prepared)
    conditional = plot_merit(TEST_DATA; ccc = ["hom", "state"], plot_type = "conditional")
    @test DiDInt.available_plot_types(conditional) == ["conditional"] && sort(unique(conditional.conditional.ccc)) == ["hom", "state"]
end

@testset "Basic run - event study" begin
    x = plot_merit(TEST_DATA; treated_states = TREATED_STATES, plot_type = "event")
    @test nrow(x.event) > 0 && DiDInt.available_plot_types(x) == ["event"]
end

# Robustness check for data with a bunch of missing values
@testset "Missing run - parallel trends" begin
    x = plot_merit(TEST_DATA_MISSING; treated_states = TREATED_STATES)
    @test nrow(x.conditional) > 0 && nrow(x.integrated) > 0 && nrow(x.integrated2) > 0 && nrow(x.conditional2) > 0
end

@testset "Missing run - event study" begin
    x = plot_merit(TEST_DATA_MISSING; treated_states = TREATED_STATES, plot_type = "event")
    @test nrow(x.event) > 0
end

# Event study rows are sorted by ccc option and by period relative to treatment
@testset "Event study rows are sorted" begin
    x = plot_merit(TEST_DATA; treated_states = TREATED_STATES, plot_type = "event", ccc = ["int", "hom"])
    @test issorted(collect(zip(x.event.ccc, x.event.time_since_treatment)))
    @test unique(x.event.ccc) == ["hom", "int"]
end
