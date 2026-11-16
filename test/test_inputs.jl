# Tests of input handling and of degenerate data: covariates without within-cell variation or
# collinear within a cell, treatment times outside the sample, period detection, duplicated
# treated states, gvar columns, float-valued arguments, date formats and the argument shapes
# that didint_plot() accepts.
# The panel is test/data/didintpy_panel.csv (see test_integration.jl): 12 states observed
# 2001-2008, s01 and s02 treated in 2004, s03 and s04 in 2006, s05 in 2007.

using Dates
using Statistics

const INPUT_TIMES = [INTEGRATION_COHORTS[s] for s in INTEGRATION_TREATED]

quiet_call(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

input_didint(data = INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES, kwargs...) = quiet_call() do
    didint_table("y", "state", "year", data; treated_states = treated_states, treatment_times = treatment_times,
                 nperm = 9, seed = 1234, kwargs...)
end

@testset "A covariate that is constant within every cell changes nothing" begin
    # A covariate measured at the state-period level is absorbed by the cell dummies. Its
    # values here are ones whose within-cell demeaning is not exact in floating point, which
    # used to leave rounding noise for the slope regressions to fit.
    panel = copy(INTEGRATION_PANEL)
    cells = collect(enumerate(unique(zip(panel.state, panel.year))))
    unemp = Dict(cell => round(3 + 6 * abs(sin(12.9898 * i)), digits = 1) for (i, cell) in cells)
    income = Dict(cell => round(30000 + 30000 * abs(sin(78.233 * i)), digits = 2) for (i, cell) in cells)
    panel.unemp = [unemp[cell] for cell in zip(panel.state, panel.year)]
    panel.income = [income[cell] for cell in zip(panel.state, panel.year)]
    inexact(col) = any(any(!iszero, cell[!, col] .- mean(cell[!, col])) for cell in groupby(panel, [:state, :year]))
    @test inexact(:unemp) || inexact(:income)

    for ccc in ["hom", "state", "time", "add", "int"]
        without = input_didint(panel; ccc = ccc, agg = "sgt")
        @test input_didint(panel; ccc = ccc, agg = "sgt", covariates = ["unemp"]).att_sgt ≈ without.att_sgt atol = 1e-10
        @test input_didint(panel; ccc = ccc, agg = "sgt", covariates = ["unemp", "income"]).att_sgt ≈ without.att_sgt atol = 1e-10
        base = input_didint(panel; ccc = ccc, agg = "sgt", covariates = ["x1", "x2"])
        both = input_didint(panel; ccc = ccc, agg = "sgt", covariates = ["x1", "x2", "unemp", "income"])
        @test both.att_sgt ≈ base.att_sgt atol = 1e-10
        @test both.agg_att[1] ≈ base.agg_att[1] atol = 1e-10
    end

    # The estimates do not depend on edgecase either (one treated state, one control)
    pair = filter(r -> r.state in ["s01", "s06"], panel)
    for ccc in ["hom", "state", "time"]
        run(edgecase) = input_didint(pair; treated_states = ["s01"], treatment_times = [2004], ccc = ccc, agg = "sgt",
                                     covariates = ["x1", "unemp"], edgecase = edgecase)
        @test run(true).att_sgt ≈ run(false).att_sgt atol = 1e-10
    end

    # State-period level data (one observation per cell): no covariate has within-cell
    # variation, so every ccc gives the estimates without covariates
    aggregated = combine(groupby(INTEGRATION_PANEL, [:state, :year]), :y => mean => :y, :x1 => mean => :x1, :x2 => mean => :x2)
    without = input_didint(aggregated; agg = "cohort")
    for ccc in ["hom", "state", "time", "add", "int"]
        @test input_didint(aggregated; ccc = ccc, agg = "cohort", covariates = ["x1", "x2"]).att_cohort ≈ without.att_cohort atol = 1e-10
    end
end

@testset "ccc = int: covariates that are collinear with the intercept within a cell" begin
    # A factor covariate whose reference category does not occur in one cell: that cell's
    # dummies add up to one. One of them must be dropped there, so that adding a constant to
    # every outcome changes no estimate.
    panel = copy(INTEGRATION_PANEL)
    panel.cat = [["u", "v", "w"][1 + mod(i, 3)] for i in 1:nrow(panel)]
    in_cell = (panel.state .== "s06") .& (panel.year .== 2005)
    panel.cat[in_cell] = [["v", "w"][1 + mod(i, 2)] for i in 1:count(in_cell)]
    cell_dummies = Float64.(hcat(panel.cat[in_cell] .== "v", panel.cat[in_cell] .== "w"))
    @test length(DiDInt.prune_covariate_columns(cell_dummies)) == 1

    # Collinearity is recognized through the rounding error that centering leaves behind, with a
    # tolerance wide enough not to depend on the processor: three dummies that add up to one,
    # and a covariate that is another one plus a large constant. Covariates that are not
    # collinear are all kept.
    n_cell = count(in_cell)
    three_dummies = Float64.(hcat([mod(i, 3) == 0 for i in 1:n_cell], [mod(i, 3) == 1 for i in 1:n_cell], [mod(i, 3) == 2 for i in 1:n_cell]))
    @test length(DiDInt.prune_covariate_columns(three_dummies)) == 2
    x1_cell = panel.x1[in_cell]
    @test length(DiDInt.prune_covariate_columns(hcat(x1_cell, 1000 .+ 3 .* x1_cell))) == 1
    @test length(DiDInt.prune_covariate_columns(hcat(x1_cell, 1e6 .+ x1_cell ./ 7))) == 1
    @test length(DiDInt.prune_covariate_columns(hcat(x1_cell, Float64.(panel.x2[in_cell])))) == 2
    @test length(DiDInt.prune_covariate_columns(hcat(x1_cell, x1_cell .^ 2))) == 2

    run(data) = input_didint(data; ccc = "int", agg = "sgt", covariates = ["cat"], ref = Dict("cat" => "u"))
    shifted = copy(panel)
    shifted.y = shifted.y .+ 100
    @test run(shifted).att_sgt ≈ run(panel).att_sgt atol = 1e-8

    # Numeric covariates with x2 an exact linear function of x1 plus a constant in one cell
    numeric = copy(INTEGRATION_PANEL)
    numeric.x2 = Float64.(numeric.x2)
    numeric.x2[in_cell] = 2 .+ 3 .* numeric.x1[in_cell]
    shifted = copy(numeric)
    shifted.y = shifted.y .+ 100
    run_numeric(data) = input_didint(data; ccc = "int", agg = "sgt", covariates = ["x1", "x2"])
    @test run_numeric(shifted).att_sgt ≈ run_numeric(numeric).att_sgt atol = 1e-8

    # A cell with as many observations as coefficients before the collinear dummy is dropped
    small = panel[.!in_cell .| (cumsum(in_cell) .<= 3), :]
    @test all(isfinite, run(small).att_sgt)
end

@testset "hc4 with weighted regressions and few controls" begin
    # The single treated row of an ATT(s,t) regression has leverage 1, which rounding can put
    # just above 1. Its HC4 standard error is NaN, as for HC2 and HC3, and nothing errors.
    controls = ["s06", "s07", "s08", "s09", "s10", "s11", "s12"]
    for n_controls in 2:7, weighting in ["both", "none"]
        data = filter(r -> r.state in vcat(INTEGRATION_TREATED, controls[1:n_controls]), INTEGRATION_PANEL)
        sgt = input_didint(data; hc = "hc4", weighting = weighting, agg = "sgt")
        @test all(isnan, sgt.se_att_sgt)
        @test isfinite(sgt.se_agg_att[1])
        @test isfinite(input_didint(data; hc = "hc4", weighting = weighting, agg = "cohort").se_agg_att[1])
    end
end

@testset "Treatment times outside the sample" begin
    # A state treated at or before the first period is treated throughout: it is neither a
    # treated state of the analysis nor a valid control, and its observations are dropped
    times = [2004, 2004, 2006, 2006, 2001]
    always_treated = input_didint(; treatment_times = times, agg = "sgt")
    without_s05 = input_didint(filter(r -> r.state != "s05", INTEGRATION_PANEL); treated_states = INTEGRATION_TREATED[1:4],
                               treatment_times = times[1:4], agg = "sgt")
    @test isequal(always_treated, without_s05)
    as_control = input_didint(; treated_states = INTEGRATION_TREATED[1:4], treatment_times = times[1:4], agg = "sgt")
    @test maximum(abs.(always_treated.att_sgt .- as_control.att_sgt)) > 1e-6
    @test_logs (:warn, r"treated at or before the first date") match_mode = :any didint("y", "state", "year", INTEGRATION_PANEL;
        treated_states = INTEGRATION_TREATED, treatment_times = times, nperm = 9, seed = 1234)
    # The same through start_date
    windowed = input_didint(; start_date = 2005, agg = "sgt")
    @test Set(windowed.state) == Set(["s03", "s04", "s05"])
    @test isequal(windowed, input_didint(filter(r -> !(r.state in ["s01", "s02"]) && r.year >= 2005, INTEGRATION_PANEL);
                                         treated_states = ["s03", "s04", "s05"], treatment_times = [2006, 2006, 2007], agg = "sgt"))

    # A state treated after the last date is untreated throughout and serves as a control
    dated = copy(INTEGRATION_PANEL)
    dated.date = Date.(dated.year)
    run_dated(states, dates; kwargs...) = quiet_call() do
        didint_table("y", "state", "date", dated; treated_states = states, treatment_times = dates, nperm = 9, seed = 1234,
                     agg = "sgt", kwargs...)
    end
    dates = Date.(INPUT_TIMES)
    later = run_dated(INTEGRATION_TREATED, vcat(dates[1:4], Date(2008, 6, 1)))
    @test isequal(later, run_dated(INTEGRATION_TREATED[1:4], dates[1:4]))
    @test !("s05" in later.state)
    @test_logs (:warn, r"treated after the last date") match_mode = :any didint("y", "state", "date", dated;
        treated_states = INTEGRATION_TREATED, treatment_times = vcat(dates[1:4], Date(2008, 6, 1)), nperm = 9, seed = 1234)

    # No treated state within the sample
    @test_throws r"None of the treated states" input_didint(; treated_states = ["s01"], treatment_times = [2001])
end

@testset "Automatic period length" begin
    years(values) = Date.(values)
    @test DiDInt.get_max_period(years(2001:2008)) == Year(1)
    @test DiDInt.get_max_period(years([2001, 2003, 2005, 2007])) == Year(2)
    @test DiDInt.get_max_period(years([1995, 2000, 2005, 2010])) == Year(5)      # gaps of 1826, 1827, 1826 days
    @test DiDInt.get_max_period(years([1990, 2000, 2010, 2020])) == Year(10)
    @test DiDInt.get_max_period(years([2001, 2007, 2013])) == Year(6)
    @test DiDInt.get_max_period(years([1999, 2006, 2013])) == Year(7)
    @test DiDInt.get_max_period(years([1896, 1900, 1904])) == Year(4)            # 1900 is not a leap year
    @test DiDInt.get_max_period([Date(2001, 1, 1), Date(2002, 2, 5)]) == Year(2)  # 400 days: rounded up
    @test DiDInt.get_max_period([Date(2001, 1, 1), Date(2001, 2, 1), Date(2001, 3, 1)]) == Month(1)
    @test DiDInt.get_max_period([Date(2001, 1, 1), Date(2001, 1, 8)]) == Week(1)

    # Five-yearly data gives the estimates of the yearly data it was relabelled from
    relabel(year) = 2001 + 5 * (year - 2001)
    spaced = transform(INTEGRATION_PANEL, :year => ByRow(relabel) => :year)
    yearly = input_didint(; agg = "sgt")
    five_yearly = quiet_call() do
        didint("y", "state", "year", spaced; treated_states = INTEGRATION_TREATED, treatment_times = relabel.(INPUT_TIMES),
               agg = "sgt", nperm = 9, seed = 1234)
    end
    @test five_yearly.specs.period == "5 years"
    @test five_yearly.sub.att ≈ yearly.att_sgt atol = 1e-10
end

@testset "A treated state can only be listed once" begin
    @test_throws r"listed more than once" input_didint(; treated_states = vcat(INTEGRATION_TREATED, "s01"), treatment_times = vcat(INPUT_TIMES, 2004))
    @test_throws r"listed more than once" input_didint(; treated_states = vcat(INTEGRATION_TREATED, "s01"), treatment_times = vcat(INPUT_TIMES, 2005))
    @test_throws r"same length" input_didint(; treated_states = INTEGRATION_TREATED, treatment_times = [2004, 2006])

    # A gvar column that is not constant within a state
    panel = copy(INTEGRATION_PANEL)
    panel.gvar = [get(INTEGRATION_COHORTS, s, 0) for s in panel.state]
    panel.gvar[findfirst(panel.state .== "s01")] = 2005
    @test_throws r"more than one value in the 'gvar' column" input_didint(panel; gvar = "gvar", treated_states = nothing, treatment_times = nothing)
end

@testset "gvar columns: dates, and rows with a missing state" begin
    reference = input_didint(; agg = "sgt")

    # A Date time column with a Date gvar column that is missing for never-treated states
    panel = copy(INTEGRATION_PANEL)
    panel.date = Date.(panel.year)
    panel.gdate = [haskey(INTEGRATION_COHORTS, s) ? Date(INTEGRATION_COHORTS[s]) : missing for s in panel.state]
    dated = quiet_call() do
        didint_table("y", "state", "date", panel; gvar = "gdate", agg = "sgt", nperm = 9, seed = 1234)
    end
    @test dated.att_sgt ≈ reference.att_sgt atol = 1e-10

    # Rows with a missing state are dropped, with gvar as with treated_states
    panel = allowmissing(copy(INTEGRATION_PANEL))
    panel.gvar = [get(INTEGRATION_COHORTS, s, 0) for s in panel.state]
    kept = trues(nrow(panel))
    kept[1:17:end] .= false
    panel.state[.!kept] .= missing
    with_gvar = input_didint(panel; gvar = "gvar", treated_states = nothing, treatment_times = nothing, agg = "sgt")
    @test with_gvar.att_sgt ≈ input_didint(INTEGRATION_PANEL[kept, :]; agg = "sgt").att_sgt atol = 1e-10
end

@testset "Numbers given as floats, and date formats" begin
    reference = input_didint(; agg = "cohort", hc = "hc1")
    @test isequal(input_didint(; agg = "cohort", hc = 1.0), reference)
    @test isequal(input_didint(; agg = "cohort", hc = 1), reference)
    @test isequal(input_didint(; agg = "cohort", freq = "year", freq_multiplier = 1.0), reference)
    @test isequal(input_didint(; agg = "cohort", start_date = 2001.0, end_date = 2008.0), reference)
    @test_throws r"freq_multiplier" input_didint(; agg = "cohort", freq = "year", freq_multiplier = 0)
    @test_throws r"'hc' must be one of" input_didint(; agg = "cohort", hc = 7)

    # Date formats are not case sensitive, and an unrecognized one is an error
    panel = copy(INTEGRATION_PANEL)
    panel.date = Date.(panel.year)
    run_dated(date_format) = quiet_call() do
        didint("y", "state", "date", panel; treated_states = INTEGRATION_TREATED, treatment_times = Date.(INPUT_TIMES),
               agg = "cohort", nperm = 9, seed = 1234, freq = "year", date_format = date_format)
    end
    lower = run_dated("yyyy-mm-dd")
    @test lower.specs.start_date == "2001-01-01" && lower.specs.end_date == "2008-01-01"
    @test isequal(run_dated("YYYY-MM-DD").table, lower.table)
    @test_throws r"date_format" run_dated("garbage")

    # String dates in the data and in treatment_times
    panel.text = string.(panel.year) .* "/06/15"
    text = quiet_call() do
        didint_table("y", "state", "text", panel; treated_states = INTEGRATION_TREATED, treatment_times = string.(INPUT_TIMES) .* "/01/01",
                     date_format = "yyyy/mm/dd", freq = "year", agg = "cohort", nperm = 9, seed = 1234)
    end
    @test text.att_cohort ≈ reference.att_cohort atol = 1e-10
end

@testset "didint_plot() accepts the argument shapes that didint() accepts" begin
    plot(; kwargs...) = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; kwargs...)
    end
    covariates = ["x1", "x2"]
    for event in [false, true]
        vectors = plot(treated_states = ["s01", "s02"], treatment_times = [2004, 2004], event = event, covariates = covariates)

        # A single treatment time for several treated states, and a single treated state
        @test isequal(plot(treated_states = ["s01", "s02"], treatment_times = 2004, event = event, covariates = covariates), vectors)
        single = plot(treated_states = "s01", treatment_times = 2004, event = event, covariates = covariates)
        @test isequal(single, plot(treated_states = ["s01"], treatment_times = [2004], event = event, covariates = covariates))

        # Covariates as symbols, or as a single name
        @test isequal(plot(treated_states = ["s01", "s02"], treatment_times = 2004, event = event, covariates = [:x1, :x2]), vectors)
        @test isequal(plot(treated_states = "s01", treatment_times = 2004, event = event, covariates = :x1),
                      plot(treated_states = "s01", treatment_times = 2004, event = event, covariates = "x1"))

        @test_throws r"same length" plot(treated_states = ["s01", "s02", "s03"], treatment_times = [2004, 2006], event = event)
        @test_throws r"listed more than once" plot(treated_states = ["s01", "s01"], treatment_times = [2004, 2004], event = event)
    end
    @test_throws r"not a valid 'ccc' option" plot(treated_states = ["s01"], treatment_times = 2004, ccc = "homm")

    # The caller's data is left as it was
    data = copy(INTEGRATION_PANEL)
    plot_data = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", data; treated_states = ["s01", "s02"], treatment_times = 2004, covariates = covariates, event = true)
    end
    result = input_didint(data; covariates = covariates)
    @test isequal(data, INTEGRATION_PANEL) && names(data) == names(INTEGRATION_PANEL)
end

@testset "Randomization inference needs every state's long difference; jackknife warnings" begin
    # s1 is the only state treated in 2003 and is not observed in 2006. That period only
    # matters for the 2005 cohort, where s1 could not stand in as a treated state: the
    # randomization inference is not performed.
    cohorts = Dict("s1" => 2003, "s3" => 2005)
    data = small_panel(cohorts = cohorts, effect = 1.0, noise = 0.3)
    gap = filter(r -> !(r.state == "s1" && r.year == 2006), data)
    @test !ismissing(quiet_didint(data, cohorts; agg = "cohort").ri_pval_agg_att[1])
    with_gap = quiet_didint(gap, cohorts; agg = "cohort")
    @test all(ismissing, with_gap.ri_pval_att_cohort) && ismissing(with_gap.ri_pval_agg_att[1]) && ismissing(with_gap.nperm[1])

    # The kernel p-value falls back to the ordinary one when the placebo estimates do not vary,
    # also when their standard deviation is only zero up to rounding
    @test DiDInt.ri_pvalue_kernel([0.2, 0.2, 0.2], 0.5) == (1 + 0) / (3 + 1)
    @test DiDInt.ri_pvalue_kernel(fill(0.1 + 0.2, 10), 0.3) == 1.0
    # A p-value needs a finite estimate
    @test ismissing(DiDInt.ri_pvalue([0.1, -0.2], NaN)) && ismissing(DiDInt.ri_pvalue_kernel([0.1, -0.2, 0.3], NaN))

    # Deleting the only treated state of a cohort leaves the aggregate defined over the other
    # cohorts: the standard error is reported, with a warning
    cohorts = Dict("s1" => 2003, "s2" => 2003, "s3" => 2005)
    data = small_panel(cohorts = cohorts, effect = 1.0, noise = 0.3)
    jackknife_warning = r"remove every treated state of a sub-aggregate"
    @test_logs (:warn, jackknife_warning) match_mode = :any didint("y", "state", "year", data; treated_states = ["s1", "s2", "s3"],
        treatment_times = [2003, 2003, 2005], agg = "cohort", nperm = 399, seed = 1234)
    @test !ismissing(quiet_didint(data, cohorts; agg = "cohort").jknifese_agg_att[1])
end
