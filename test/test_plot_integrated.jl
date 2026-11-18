# Tests of the integrated parallel trends data of didint_plot().

using Dates
using LinearAlgebra
using Statistics

# The integrated parallel trends data (Karim & Webb 2025, eqs. 36-38), against a direct computation.
# The panel is INTEGRATION_PANEL (see test_integration.jl): 12 states observed 2001-2008, with
# s01, s02 treated in 2004, s03, s04 in 2006 and s05 in 2007.
@testset "Integrated parallel trends data" begin
    plot_data(; kwargs...) = DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED,
                                                treatment_times = INPUT_TIMES, covariates = INTEGRATION_COVARIATES, kwargs...)
    never_treated = setdiff(sort(unique(INTEGRATION_PANEL.state)), INTEGRATION_TREATED)
    treatment_year = Dict(zip(INTEGRATION_TREATED, INPUT_TIMES))
    # Regression of the outcome on the covariates in one cell, evaluated at covariate values x
    predicted(state, yr, x) = (fit = cell_ols(INTEGRATION_PANEL, state, yr).beta; fit[1] + dot(fit[2:end], x))
    delta(state, target, yr) = (x = cell_xbar(INTEGRATION_PANEL, target, yr); predicted(state, yr, x) - predicted(state, yr - 1, x))

    integrated = plot_data()
    @test names(integrated) == ["treated_state", "state", "time", "delta", "ccc", "period", "treat_period", "start_date", "period_length"]
    @test all(integrated.ccc .== "int")
    # One row per treated state, state (itself and the never-treated states) and period from the second on
    @test nrow(integrated) == length(INTEGRATION_TREATED) * (1 + length(never_treated)) * 7
    @test sort(unique(integrated.treated_state)) == INTEGRATION_TREATED
    for target in INTEGRATION_TREATED
        rows = integrated[integrated.treated_state .== target, :]
        @test sort(unique(rows.state)) == sort(vcat(target, never_treated))
        @test sort(unique(rows.time)) == string.(2002:2008)
        # periods are numbered from 0, and the treated state's treatment period is on every row
        @test all(rows.period .== parse.(Int, rows.time) .- 2001)
        @test all(rows.treat_period .== treatment_year[target] - 2001)
    end
    # Every value is the change in the cell regressions' predictions at the treated state's covariate means
    @test all(isapprox(row.delta, delta(row.state, row.treated_state, parse(Int, row.time)); atol = 1e-8) for row in eachrow(integrated))

    # Relation to didint(): summing the treated state's own changes and its controls' from the
    # period after the base period does not give the long difference (each change is evaluated
    # at a different period's covariate means), but in the first treated period the change is
    # the long difference itself, so ATT(s, g) is the treated state's change minus the weighted
    # mean of its controls' changes
    sgt = integration_didint(ccc = "int", agg = "sgt", weighting = "none")
    for target in INTEGRATION_TREATED
        g = treatment_year[target]
        rows = integrated[(integrated.treated_state .== target) .& (integrated.time .== string(g)), :]
        own = only(rows.delta[rows.state .== target])
        controls = rows.delta[rows.state .!= target]
        att = only(sgt.att_sgt[(sgt.state .== target) .& (year.(sgt.t) .== g)])
        @test own - mean(controls) ≈ att atol = 1e-8
    end

    # notyet = true adds the states not yet treated in the period: for s01 (treated 2004), s03 and
    # s04 until 2005, s05 until 2006, and the other state of its own cohort (s02) until 2003
    with_notyet = plot_data(notyet = true)
    rows = with_notyet[with_notyet.treated_state .== "s01", :]
    @test sort(unique(rows.time[rows.state .== "s03"])) == string.(2002:2005)
    @test sort(unique(rows.time[rows.state .== "s05"])) == string.(2002:2006)
    @test sort(unique(rows.time[rows.state .== "s02"])) == string.(2002:2003)
    @test all(isapprox(row.delta, delta(row.state, row.treated_state, parse(Int, row.time)); atol = 1e-8) for row in eachrow(with_notyet))
    # the rows without notyet are unchanged
    in_default = in.(with_notyet.state, Ref(never_treated)) .| (with_notyet.state .== with_notyet.treated_state)
    @test with_notyet[in_default, :] == integrated

    # Without covariates the changes are those of the cell means, whatever the treated state
    no_covariates = plot_data(covariates = nothing)
    cell_mean(state, yr) = mean(cell_rows(INTEGRATION_PANEL, state, yr).y)
    @test all(isapprox(row.delta, cell_mean(row.state, parse(Int, row.time)) - cell_mean(row.state, parse(Int, row.time) - 1); atol = 1e-10)
              for row in eachrow(no_covariates))

    # integrated = false gives the conditional parallel trends data: the lambdas of every ccc option
    conditional = plot_data(integrated = false)
    @test "lambda" in names(conditional) && !("delta" in names(conditional))
    @test sort(unique(skipmissing(conditional.ccc))) == ["add", "hom", "int", "none", "state", "time"]
    # ... and does not need the treated states
    @test isequal(DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES,
                                     covariates = INTEGRATION_COVARIATES, integrated = false), conditional)
    # The event study ignores `integrated`
    @test isequal(plot_data(event = true), plot_data(event = true, integrated = false))

    # The integrated data need the treated states, and ignore `ccc` with a warning
    @test_throws r"need the treated states" DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES)
    @test_logs (:warn, r"'ccc' is ignored") match_mode = :any plot_data(ccc = "hom")
    @test quiet_call(() -> plot_data(ccc = "hom")) == integrated
    # ... except "int", which is what the integrated data use: no warning
    @test_logs plot_data(ccc = "int")
    @test plot_data(ccc = "int") == integrated
    # gvar gives the treated states as well
    with_gvar = copy(INTEGRATION_PANEL)
    with_gvar.gvar = [get(treatment_year, s, 0) for s in with_gvar.state]
    @test DiDInt.didint_plot("y", "state", "year", with_gvar; gvar = "gvar", covariates = INTEGRATION_COVARIATES) == integrated

    # A state that is not observed in a period has no change into or out of that period, and a
    # treated state that is not observed in a period has no rows for it
    gap = filter(r -> !(r.state == "s07" && r.year == 2004) && !(r.state == "s01" && r.year == 2006), INTEGRATION_PANEL)
    with_gap = DiDInt.didint_plot("y", "state", "year", gap; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                                  covariates = INTEGRATION_COVARIATES)
    s07 = with_gap[(with_gap.treated_state .== "s02") .& (with_gap.state .== "s07"), :]
    @test sort(s07.time) == string.([2002, 2003, 2006, 2007, 2008])
    s01 = with_gap[with_gap.treated_state .== "s01", :]
    @test !("2006" in s01.time)
    @test sort(unique(s01.time[s01.state .== "s01"])) == string.([2002, 2003, 2004, 2005, 2008])
    @test "2007" in s01.time[s01.state .== "s06"]
end
