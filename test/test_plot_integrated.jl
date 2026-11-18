# Tests of the parallel trends data of didint_plot(): the conditional parallel trends data
# (every ccc option except "int") and the integrated parallel trends data (ccc = "int"), as
# they are passed in one table to the Stata and R wrappers (plot_table, see runtests.jl).

using Dates
using LinearAlgebra
using Statistics

# The integrated parallel trends data (Karim & Webb 2025, eqs. 36-38), against a direct computation.
# The panel is INTEGRATION_PANEL (see test_integration.jl): 12 states observed 2001-2008, with
# s01, s02 treated in 2004, s03, s04 in 2006 and s05 in 2007.
@testset "Integrated parallel trends data (ccc = int)" begin
    plot_data(; kwargs...) = plot_table("y", "state", "year", INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED,
                                                treatment_times = INPUT_TIMES, covariates = INTEGRATION_COVARIATES,
                                                ccc = "int", kwargs...)
    never_treated = setdiff(sort(unique(INTEGRATION_PANEL.state)), INTEGRATION_TREATED)
    treatment_year = Dict(zip(INTEGRATION_TREATED, INPUT_TIMES))
    # Regression of the outcome on the covariates in one cell, evaluated at covariate values x
    predicted(state, yr, x) = (fit = cell_ols(INTEGRATION_PANEL, state, yr).beta; fit[1] + dot(fit[2:end], x))
    delta(state, target, yr) = (x = cell_xbar(INTEGRATION_PANEL, target, yr); predicted(state, yr, x) - predicted(state, yr - 1, x))
    cell_n(state, yr) = nrow(cell_rows(INTEGRATION_PANEL, state, yr))

    integrated = plot_data()
    @test names(integrated) == ["state", "time", "lambda", "delta", "n", "treated_state_covars", "ccc", "period",
                                "start_date", "treat_period", "period_length"]
    @test all(integrated.ccc .== "int") && all(ismissing, integrated.lambda)
    # One row per treated state, state (itself and the never-treated states) and period from the
    # second on; no extra rows for the treatment periods, as every row has the treatment period
    @test nrow(integrated) == length(INTEGRATION_TREATED) * (1 + length(never_treated)) * 7
    @test sort(unique(integrated.treated_state_covars)) == INTEGRATION_TREATED
    for target in INTEGRATION_TREATED
        rows = integrated[integrated.treated_state_covars .== target, :]
        @test sort(unique(rows.state)) == sort(vcat(target, never_treated))
        @test sort(unique(rows.time)) == string.(2002:2008)
        # periods are numbered from 0, and the treated state's treatment period is on every row
        @test all(rows.period .== parse.(Int, rows.time) .- 2001)
        @test all(rows.treat_period .== treatment_year[target] - 2001)
    end
    # Every value is the change in the cell regressions' predictions at the covariate means of
    # the treated state, and n counts the state's observations in the two periods
    @test all(isapprox(row.delta, delta(row.state, row.treated_state_covars, parse(Int, row.time)); atol = 1e-8) for row in eachrow(integrated))
    @test all(row.n == cell_n(row.state, parse(Int, row.time)) + cell_n(row.state, parse(Int, row.time) - 1) for row in eachrow(integrated))

    # Relation to didint(): in the first treated period the change is the integrated long
    # difference itself, so ATT(s, g) is the treated state's change minus the mean of its
    # controls' changes: equally weighted under weighting = "none", weighted by n under "both"
    for (weighting, weight) in [("none", rows -> ones(nrow(rows))), ("both", rows -> Float64.(rows.n))]
        sgt = integration_didint(ccc = "int", agg = "sgt", weighting = weighting)
        for target in INTEGRATION_TREATED
            g = treatment_year[target]
            rows = integrated[(integrated.treated_state_covars .== target) .& (integrated.time .== string(g)), :]
            own = only(rows.delta[rows.state .== target])
            controls = rows[rows.state .!= target, :]
            w = weight(controls)
            att = only(sgt.att_sgt[(sgt.state .== target) .& (year.(sgt.t) .== g)])
            @test own - sum(w .* controls.delta) / sum(w) ≈ att atol = 1e-8
        end
    end

    # notyet = true adds the states of later cohorts while they are not yet treated, as in
    # didint(): for s01 (treated 2004), s03 and s04 until 2005 and s05 until 2006, but not the
    # other state of its own cohort (s02)
    with_notyet = plot_data(notyet = true)
    rows = with_notyet[with_notyet.treated_state_covars .== "s01", :]
    @test sort(unique(rows.time[rows.state .== "s03"])) == string.(2002:2005)
    @test sort(unique(rows.time[rows.state .== "s05"])) == string.(2002:2006)
    @test !("s02" in rows.state[rows.state .!= "s01"])
    # ... and for s03 (treated 2006) only s05, not the states treated earlier
    rows = with_notyet[with_notyet.treated_state_covars .== "s03", :]
    @test sort(unique(rows.state)) == sort(vcat("s03", "s05", never_treated))
    @test all(isapprox(row.delta, delta(row.state, row.treated_state_covars, parse(Int, row.time)); atol = 1e-8) for row in eachrow(with_notyet))
    # the rows without notyet are unchanged
    in_default = in.(with_notyet.state, Ref(never_treated)) .| (with_notyet.state .== with_notyet.treated_state_covars)
    @test isequal(with_notyet[in_default, :], integrated)

    # Without covariates the changes are those of the cell means, whatever the treated state
    no_covariates = plot_data(covariates = nothing)
    cell_mean(state, yr) = mean(cell_rows(INTEGRATION_PANEL, state, yr).y)
    @test all(isapprox(row.delta, cell_mean(row.state, parse(Int, row.time)) - cell_mean(row.state, parse(Int, row.time) - 1); atol = 1e-10)
              for row in eachrow(no_covariates))

    # gvar gives the treated states as well
    treated_gvar = copy(INTEGRATION_PANEL)
    treated_gvar.gvar = [get(treatment_year, s, 0) for s in treated_gvar.state]
    @test isequal(plot_table("y", "state", "year", treated_gvar; gvar = "gvar", covariates = INTEGRATION_COVARIATES, ccc = "int"),
                  integrated)

    # Without the treated states: only the conditional data are prepared, with a warning; the
    # integrated data asked for by name are an error
    @test_throws r"need the treated states" DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES,
                                                              ccc = ["hom", "int"], plot_type = "integrated")
    @test_logs (:warn, r"Only the conditional parallel trends data are prepared") match_mode = :any plot_table("y", "state", "year", INTEGRATION_PANEL;
                                                                                                             treatment_times = INPUT_TIMES)
    without_treated = quiet_call(() -> plot_table("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES))
    @test sort(unique(skipmissing(without_treated.ccc))) == ["add", "hom", "none", "state", "time"]

    # A state that is not observed in a period has no change into or out of that period, and a
    # treated state that is not observed in a period has no rows for it
    gap = filter(r -> !(r.state == "s07" && r.year == 2004) && !(r.state == "s01" && r.year == 2006), INTEGRATION_PANEL)
    with_gap = plot_table("y", "state", "year", gap; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                                  covariates = INTEGRATION_COVARIATES, ccc = "int")
    s07 = with_gap[(with_gap.treated_state_covars .== "s02") .& (with_gap.state .== "s07"), :]
    @test sort(s07.time) == string.([2002, 2003, 2006, 2007, 2008])
    s01 = with_gap[with_gap.treated_state_covars .== "s01", :]
    @test !("2006" in s01.time)
    @test sort(unique(s01.time[s01.state .== "s01"])) == string.([2002, 2003, 2004, 2005, 2008])
    @test "2007" in s01.time[s01.state .== "s06"]

    # A state treated at or before the first date is neither a treated state nor a comparison state
    early = quiet_call() do
        plot_table("y", "state", "year", INTEGRATION_PANEL; treated_states = vcat(INTEGRATION_TREATED, "s06"),
                           treatment_times = vcat(INPUT_TIMES, 2001), covariates = INTEGRATION_COVARIATES, ccc = "int")
    end
    @test !("s06" in early.state) && !("s06" in early.treated_state_covars)
    @test isequal(early, integrated[integrated.state .!= "s06", :])
end

@testset "Conditional and integrated parallel trends data together" begin
    plot_data(; kwargs...) = plot_table("y", "state", "year", INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED,
                                                treatment_times = INPUT_TIMES, covariates = INTEGRATION_COVARIATES, kwargs...)
    all_data = plot_data()
    integrated = plot_data(ccc = "int")
    conditional = plot_data(ccc = ["hom", "time", "state", "add", "none"])
    @test names(all_data) == names(integrated) == names(conditional)

    # The default has the conditional rows of every other ccc option, then the integrated rows
    # of "int", then the rows that mark the treatment periods
    @test sort(unique(skipmissing(all_data.ccc))) == ["add", "hom", "int", "none", "state", "time"]
    @test isequal(all_data[isequal.(all_data.ccc, "int"), :], integrated)
    @test isequal(all_data[.!isequal.(all_data.ccc, "int"), :], conditional)

    # Conditional rows: one per state, period and ccc option, with lambda and the cell's number
    # of observations; delta and treated_state_covars are missing
    rows = conditional[.!ismissing.(conditional.ccc), :]
    @test nrow(rows) == 5 * 12 * 8
    @test all(ismissing, rows.delta) && all(ismissing, rows.treated_state_covars) && all(ismissing, rows.treat_period)
    @test all(row.n == nrow(cell_rows(INTEGRATION_PANEL, row.state, parse(Int, row.time))) for row in eachrow(rows))
    none = rows[rows.ccc .== "none", :]
    @test all(isapprox(row.lambda, mean(cell_rows(INTEGRATION_PANEL, row.state, parse(Int, row.time)).y); atol = 1e-10) for row in eachrow(none))
    # followed by one row per treatment period, with only treat_period filled
    marks = conditional[ismissing.(conditional.ccc), :]
    @test sort(marks.treat_period) == [3, 5, 6]
    @test all(ismissing, marks.state) && all(ismissing, marks.lambda)
    # ... which the integrated data alone do not need
    @test !any(ismissing, integrated.ccc)

    # The conditional data do not need the treated states
    without_treated = quiet_call() do
        plot_table("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES, covariates = INTEGRATION_COVARIATES,
                   ccc = ["hom", "time", "state", "add", "none"])
    end
    @test isequal(without_treated, conditional)

    # The event study still uses the lambdas of "int"
    event = plot_data(event = true, ccc = ["hom", "int"])
    @test unique(event.ccc) == ["hom", "int"] && "y" in names(event)
end
