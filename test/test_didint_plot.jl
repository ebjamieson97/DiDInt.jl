# Tests for didint_plot()

# Ensure that didint_plot() works with good data
@testset "Basic run - parallel trends" begin
    test_plot_parallel = DiDInt.didint_plot("coll", "state", "year", TEST_DATA,
                                            treatment_times = TREATED_TIMES, 
                                            event = false, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_parallel)
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
                                            event = false, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_parallel)
end

@testset "Missing run - event study" begin
    test_plot_event = DiDInt.didint_plot("coll", "state", "year", TEST_DATA_MISSING,
                                            treatment_times = TREATED_TIMES, 
                                            treated_states = TREATED_STATES,
                                            event = true, covariates = ["asian", "black", "male"])
    @test !isnothing(test_plot_event)
end
# With covariates, every state-period is evaluated at a common covariate profile x*, the treated
# states' mean covariates, so the plot does not depend on how the covariates are centred
@testset "Plot values are fitted outcomes at a common covariate profile" begin
    panel = CSV.read(joinpath(@__DIR__, "data", "didintpy_panel.csv"), DataFrame)
    treated = ["s01", "s02", "s03", "s04", "s05"]
    times = [2004, 2004, 2006, 2006, 2007]
    plot_values(data; event = false) = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        DiDInt.didint_plot("y", "state", "year", data; treated_states = treated, treatment_times = times,
                           covariates = ["x1", "x2"], ccc = ["int", "hom", "time", "state", "add"], event = event)
    end

    # Shifting a covariate by a constant leaves every plotted value unchanged, for every ccc
    shifted = copy(panel)
    shifted.x1 .+= 10.0
    for event in [false, true]
        base = plot_values(panel; event = event)
        moved = plot_values(shifted; event = event)
        column = event ? :y : :lambda
        keep = .!ismissing.(base[!, column])
        @test Float64.(moved[keep, column]) ≈ Float64.(base[keep, column]) atol = 1e-8
    end

    # Under "int", a state-period's value is its own regression's prediction at x*
    values = plot_values(panel)
    treated_rows = in.(panel.state, Ref(treated))
    xstar = [mean(panel.x1[treated_rows]), mean(panel.x2[treated_rows])]
    for (s, yr) in [("s01", 2005), ("s07", 2003)]
        cell = panel[(panel.state .== s) .& (panel.year .== yr), :]
        beta = hcat(ones(nrow(cell)), cell.x1, cell.x2) \ cell.y
        value = values.lambda[isequal.(values.ccc, "int") .& isequal.(values.state, s) .& isequal.(values.time, string(yr))]
        @test only(value) ≈ beta[1] + beta[2] * xstar[1] + beta[3] * xstar[2] atol = 1e-8
    end
end
