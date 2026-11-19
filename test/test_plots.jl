# Tests of DiDIntPlot, the object that didint_plot() returns: its tables of plot data, what each
# kind of plot shows (the layouts of src/plots.jl) and that Makie draws them.
# The panel is INTEGRATION_PANEL (see test_integration.jl): 12 states observed 2001-2008 (periods 0
# to 7), with s01, s02 treated in 2004, s03, s04 in 2006 and s05 in 2007; s06 to s12 are never treated.

using Dates
using LinearAlgebra
using Statistics
import CairoMakie

plot_object(; kwargs...) = quiet_call() do
    DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                       covariates = INTEGRATION_COVARIATES, kwargs...)
end
const PLOT_NEVER_TREATED = setdiff(sort(unique(INTEGRATION_PANEL.state)), INTEGRATION_TREATED)
const PLOT_TREATMENT_PERIODS = INPUT_TIMES .- 2001
const PLOT_TREATMENT_PERIOD = Dict(zip(INTEGRATION_TREATED, PLOT_TREATMENT_PERIODS))
# Titles of the panels of the treated states
const PLOT_STATE_TITLES = [string(s, " (treated ", g, ")") for (s, g) in zip(INTEGRATION_TREATED, INPUT_TIMES)]
const PLOT_TYPES = ["event", "conditional", "conditional2", "integrated", "integrated2"]

# The control states of the treated state s in a period: the never-treated states and, with
# notyet, the states of later cohorts that are not yet treated, as in didint()
function plot_controls(s, period, notyet)
    notyet || return PLOT_NEVER_TREATED
    later = [j for j in INTEGRATION_TREATED if PLOT_TREATMENT_PERIOD[j] > max(period, PLOT_TREATMENT_PERIOD[s])]
    return vcat(PLOT_NEVER_TREATED, later)
end

# The points of one panel of a layout, optionally of one line or of one kind of line
function panel_series(layout, panel; group = nothing, role = nothing)
    keep = layout.series.panel .== panel
    isnothing(group) || (keep = keep .& (layout.series.group .== group))
    isnothing(role) || (keep = keep .& (layout.series.role .== role))
    return layout.series[keep, :]
end

# The lambdas and numbers of observations of one ccc option in the conditional table, by (state, period)
function conditional_values(plot, ccc)
    rows = plot.conditional[plot.conditional.ccc .== ccc, :]
    return Dict(zip(zip(rows.state, rows.period), rows.lambda)), Dict(zip(zip(rows.state, rows.period), rows.n))
end

# In a layout of changes since the pre-treatment period, the treated state's change minus the
# average change of its controls is ATT(s,t) of didint() in every treated period, and both are
# zero in the pre-treatment period
function test_gaps_are_atts(layout, sgt)
    for (panel, s) in enumerate(INTEGRATION_TREATED)
        own = panel_series(layout, panel; role = :treated)
        average = panel_series(layout, panel; role = :control_mean)
        @test own.x == 0:7 && average.x == 0:7 && all(own.group .== s)
        base = PLOT_TREATMENT_PERIOD[s] - 1
        @test only(own.y[own.x .== base]) == 0 && only(average.y[average.x .== base]) == 0
        gap(period) = only(own.y[own.x .== period]) - only(average.y[average.x .== period])
        estimates = sgt[sgt.state .== s, :]
        @test nrow(estimates) == 8 - PLOT_TREATMENT_PERIOD[s]
        @test all(isapprox(gap(year(row.t) - 2001), row.att_sgt; atol = 1e-8) for row in eachrow(estimates))
    end
end

@testset "didint_plot() returns a DiDIntPlot with a table for each kind of plot" begin
    x = plot_object()
    @test x isa DiDIntPlot
    @test DiDInt.available_plot_types(x) == PLOT_TYPES && x.specs.plot_types == PLOT_TYPES
    @test names(x.conditional) == ["state", "time", "period", "ccc", "lambda", "n"]
    @test names(x.conditional2) == ["state", "time", "period", "ccc", "delta", "n", "treated_state", "base_period", "treat_period"]
    @test names(x.integrated) == ["state", "time", "period", "delta", "n", "treated_state_covars", "treat_period"]
    @test names(x.integrated2) == ["state", "time", "period", "delta", "n", "treated_state_covars", "base_period", "treat_period"]
    @test names(x.event) == ["ccc", "time_since_treatment", "y", "se", "ci_lower", "ci_upper", "ngroup"]
    @test names(x.event_lambdas) == ["ccc", "state", "time_since_treatment", "lambda", "n"]
    # no table has missing values, except the event study where an interval cannot be computed
    @test all(!any(ismissing, column) for table in (x.conditional, x.conditional2, x.integrated, x.integrated2) for column in eachcol(table))
    # dates and period numbers go together
    @test all(all(table.time .== string.(2001 .+ table.period)) for table in (x.conditional, x.conditional2, x.integrated, x.integrated2))

    # The specification that the plots use
    specs = x.specs
    @test specs.outcome == "y" && !specs.notyet && specs.covariates == INTEGRATION_COVARIATES
    @test specs.treated_states == INTEGRATION_TREATED && specs.treatment_periods == PLOT_TREATMENT_PERIODS
    @test specs.treatment_times == string.(INPUT_TIMES) && isempty(specs.always_treated)
    @test specs.period == "1 year" && specs.start_date == "2001" && specs.end_date == "2008"
    @test specs.period_labels == string.(2001:2008) && sort(specs.ccc) == ["add", "hom", "int", "none", "state", "time"]
    @test plot_object(notyet = true).specs.notyet

    # plot_type selects the kinds of plot data: only their tables are filled, with the same data
    some = plot_object(plot_type = ["integrated2", "Event"])
    @test DiDInt.available_plot_types(some) == ["event", "integrated2"] && some.specs.plot_types == ["event", "integrated2"]
    @test isequal(some.event, x.event) && isequal(some.event_lambdas, x.event_lambdas) && isequal(some.integrated2, x.integrated2)
    @test all(nrow(table) == 0 for table in (some.conditional, some.conditional2, some.integrated))
    one = plot_object(plot_type = "conditional2")
    @test DiDInt.available_plot_types(one) == ["conditional2"] && isequal(one.conditional2, x.conditional2)
    @test sort(one.specs.ccc) == ["add", "hom", "none", "state", "time"]
    @test isequal(plot_object(plot_type = "integrated").integrated, x.integrated)
    @test isequal(plot_object(plot_type = "conditional").conditional, x.conditional)

    # With "all", the kinds that the ccc options do not allow are left out; asked for by name, they are an error
    @test DiDInt.available_plot_types(plot_object(ccc = "hom")) == ["event", "conditional", "conditional2"]
    @test DiDInt.available_plot_types(plot_object(ccc = "int")) == ["event", "integrated", "integrated2"]
    @test unique(plot_object(ccc = "int").event.ccc) == ["int"]
    @test_throws r"needs \"int\" among the 'ccc' options" plot_object(ccc = "hom", plot_type = "integrated2")
    @test_throws r"needs a 'ccc' option other than \"int\"" plot_object(ccc = "int", plot_type = ["event", "conditional"])
    @test_throws r"not a valid 'plot_type'" plot_object(plot_type = "histogram")
    @test_throws r"not a valid 'ccc' option" plot_object(ccc = "homm")

    # Without the treated states: only the conditional data, with a warning under "all"; an
    # error if another kind is asked for by name
    without_treated(; kwargs...) = DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES,
                                                      covariates = INTEGRATION_COVARIATES, ccc = ["hom", "none"], kwargs...)
    @test_logs (:warn, r"Only the conditional parallel trends data are prepared") without_treated()
    unknown = @test_logs without_treated(plot_type = "conditional")
    @test DiDInt.available_plot_types(unknown) == ["conditional"]
    @test isempty(unknown.specs.treated_states) && unknown.specs.treatment_periods == [3, 5, 6]
    @test unknown.specs.treatment_times == ["2004", "2006", "2007"]
    @test_throws r"need the treated states" without_treated(plot_type = ["conditional", "event"])
    @test_throws r"need the treated states" without_treated(plot_type = "integrated", ccc = "int")

    # A state treated at or before the first date is remembered as such
    early = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = vcat(INTEGRATION_TREATED, "s06"),
                           treatment_times = vcat(INPUT_TIMES, 2001), ccc = "none")
    end
    @test early.specs.treated_states == INTEGRATION_TREATED && early.specs.always_treated == ["s06"]

    # Dates are strings in every table and in the specifications, as the Stata and R wrappers need
    @test !any(eltype(column) <: DiDInt.Dates.TimeType for field in DiDInt.PLOT_TYPES for column in eachcol(getfield(x, Symbol(field))))
    @test !any(eltype(column) <: DiDInt.Dates.TimeType for column in eachcol(x.event_lambdas))
    @test !any(value isa DiDInt.Dates.TimeType || (value isa AbstractVector && eltype(value) <: DiDInt.Dates.TimeType) for value in values(x.specs))

    # Display
    text = sprint(show, MIME("text/plain"), x)
    @test occursin("DID-INT plot data", text) && occursin("s01 (2004)", text) && occursin("never-treated states", text)
    @test all(occursin("  " * rpad(type, 12) * "   " * string(nrow(DiDInt.plot_table(x, type))) * " rows", text) for type in PLOT_TYPES)
    @test all(length(line) <= 110 for line in split(text, "\n"))
    @test DiDInt.wrapped_list(["aa", "bb", "cc"], 2; width = 6) == "aa, bb,\n  cc"
    @test occursin("not specified (treatment times: 2004, 2006, 2007)", sprint(show, MIME("text/plain"), unknown))
    @test sprint(show, some) == "DiDIntPlot(plot data for: event, integrated2)"
end

@testset "Tables of changes: conditional2, integrated and integrated2" begin
    predicted(state, yr, x) = (fit = cell_ols(INTEGRATION_PANEL, state, yr).beta; fit[1] + dot(fit[2:end], x))
    cell_n(state, period) = nrow(cell_rows(INTEGRATION_PANEL, state, 2001 + period))

    for notyet in [false, true]
        x = plot_object(notyet = notyet)

        # conditional2: the change in lambda since the pre-treatment period of the treated state,
        # for that state and its control states, in every other period
        for option in ["none", "hom", "time", "state", "add"]
            lambda, n = conditional_values(x, option)
            rows = x.conditional2[x.conditional2.ccc .== option, :]
            @test all(row.delta == lambda[(row.state, row.period)] - lambda[(row.state, row.base_period)] for row in eachrow(rows))
            @test all(row.n == n[(row.state, row.period)] + n[(row.state, row.base_period)] for row in eachrow(rows))
            expected = [(s, j, period) for s in INTEGRATION_TREATED for period in 0:7 if period != PLOT_TREATMENT_PERIOD[s] - 1
                        for j in vcat(s, plot_controls(s, period, notyet))]
            @test collect(zip(rows.treated_state, rows.state, rows.period)) == expected
        end
        @test all(x.conditional2.base_period .== x.conditional2.treat_period .- 1)
        @test all(row.treat_period == PLOT_TREATMENT_PERIOD[row.treated_state] for row in eachrow(x.conditional2))

        # integrated: the change in the predicted outcome since the previous period, at the
        # covariate means of the treated state; integrated2: the same since the pre-treatment period
        for (table, has_base) in [(x.integrated, false), (x.integrated2, true)]
            for row in eachrow(table)
                reference = has_base ? row.base_period : row.period - 1
                covariate_means = cell_xbar(INTEGRATION_PANEL, row.treated_state_covars, 2001 + row.period)
                @test row.delta ≈ predicted(row.state, 2001 + row.period, covariate_means) -
                                  predicted(row.state, 2001 + reference, covariate_means) atol = 1e-8
                @test row.n == cell_n(row.state, row.period) + cell_n(row.state, reference)
                @test row.treat_period == PLOT_TREATMENT_PERIOD[row.treated_state_covars]
            end
            periods(s) = has_base ? [period for period in 0:7 if period != PLOT_TREATMENT_PERIOD[s] - 1] : 1:7
            expected = [(s, j, period) for s in INTEGRATION_TREATED for period in periods(s) for j in vcat(s, plot_controls(s, period, notyet))]
            @test collect(zip(table.treated_state_covars, table.state, table.period)) == expected
        end
        @test all(x.integrated2.base_period .== x.integrated2.treat_period .- 1)
    end

    # The control states with notyet = true are those of didint(): the states of later cohorts
    # while they are not yet treated. For s01 (treated in period 3): s03 and s04 until period 4,
    # s05 until period 5, and never s02 (the same cohort). For s03 (treated in period 5): s05
    # in period 5 and before, and never s01, s02 (treated earlier) or s04 (the same cohort).
    x = plot_object(notyet = true, ccc = "int", plot_type = ["integrated", "integrated2"])
    rows(table, s) = table[(table.treated_state_covars .== s) .& (table.state .!= s), :]
    s01 = rows(x.integrated, "s01")
    @test s01.period[s01.state .== "s03"] == 1:4 && s01.period[s01.state .== "s05"] == 1:5 && !("s02" in s01.state)
    s03 = rows(x.integrated2, "s03")
    @test s03.period[s03.state .== "s05"] == [0, 1, 2, 3, 5] && isempty(intersect(s03.state, ["s01", "s02", "s04"]))

    # In the first treated period the change since the pre-treatment period is the change since
    # the previous period
    x = plot_object()
    for s in INTEGRATION_TREATED
        g = PLOT_TREATMENT_PERIOD[s]
        @test x.integrated.delta[(x.integrated.treated_state_covars .== s) .& (x.integrated.period .== g)] ≈
              x.integrated2.delta[(x.integrated2.treated_state_covars .== s) .& (x.integrated2.period .== g)] atol = 1e-10
    end

    # A treated state that is not observed in its pre-treatment period has no changes since that period
    gap = filter(r -> !(r.state == "s01" && r.year == 2003), INTEGRATION_PANEL)
    with_gap = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", gap; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                           covariates = INTEGRATION_COVARIATES, plot_type = ["conditional2", "integrated2"])
    end
    @test !("s01" in with_gap.conditional2.treated_state) && !("s01" in with_gap.integrated2.treated_state_covars)
    @test "s02" in with_gap.conditional2.treated_state && "s02" in with_gap.integrated2.treated_state_covars
end

@testset "Event study plot" begin
    x = plot_object()
    options = ["none", "hom", "time", "state", "add", "int"]

    # The first stage is that of the other plot data and of didint(), on all states: the lambdas
    # of the treated states are those of the conditional table
    for option in ["none", "hom", "time", "state", "add"]
        lambda, n = conditional_values(x, option)
        rows = x.event_lambdas[x.event_lambdas.ccc .== option, :]
        @test nrow(rows) == 5 * 8
        @test all(row.lambda == lambda[(row.state, PLOT_TREATMENT_PERIOD[row.state] + row.time_since_treatment)] for row in eachrow(rows))
        @test all(row.n == n[(row.state, PLOT_TREATMENT_PERIOD[row.state] + row.time_since_treatment)] for row in eachrow(rows))
    end
    @test sort(unique(x.event_lambdas.state)) == INTEGRATION_TREATED

    # The estimates are the means of the treated states' lambdas, weighted by their numbers of observations
    for group in groupby(x.event_lambdas, [:ccc, :time_since_treatment])
        row = only(eachrow(x.event[(x.event.ccc .== first(group.ccc)) .& (x.event.time_since_treatment .== first(group.time_since_treatment)), :]))
        @test row.y ≈ sum(group.n .* group.lambda) / sum(group.n) atol = 1e-12
        @test row.ngroup == nrow(group)
    end
    @test issorted(collect(zip(x.event.ccc, x.event.time_since_treatment)))

    # One panel per ccc option with the estimates of the data; an interval needs 3 treated states
    layout = DiDInt.plot_layout(x; plot_type = "event")
    @test layout.plot_type == "event" && layout.title == "Event study" && layout.subtitle == "" && isnothing(layout.xticks)
    @test layout.panels == [DiDInt.PLOT_CCC_LABELS[option] for option in options]
    @test layout.vlines.panel == 1:6 && all(layout.vlines.x .== 0)
    for (panel, option) in enumerate(options)
        rows = x.event[x.event.ccc .== option, :]
        series = panel_series(layout, panel)
        @test all(series.role .== :estimate) && series.x == rows.time_since_treatment && series.y == rows.y
        @test all(ismissing(lower) ? ngroup < 3 : lower == ci for (lower, ci, ngroup) in zip(series.lower, rows.ci_lower, rows.ngroup))
        @test all(ismissing(upper) ? ngroup < 3 : upper == ci for (upper, ci, ngroup) in zip(series.upper, rows.ci_upper, rows.ngroup))
    end
    @test any(ismissing, layout.series.lower) && any(!ismissing, layout.series.lower)
    @test first.(layout.legend_extras) == [:band, :no_interval] && layout.legend_extras[1][2] == "95% confidence interval"
    with_two = DiDInt.plot_layout(x; plot_type = "event", groupmin = 2)
    @test isequal(with_two.series.lower, x.event.ci_lower[sortperm([findfirst(==(c), options) for c in x.event.ccc])])
    wider = plot_object(plot_type = "event", ci = 0.9)
    @test wider.event.ci_lower != x.event.ci_lower && wider.event.y == x.event.y
    @test DiDInt.plot_layout(wider).legend_extras[1][2] == "90% confidence interval"
    @test DiDInt.plot_layout(wider).plot_type == "event"
    @test_throws r"'ci' must be between 0 and 1" plot_object(ci = 1.5)

    # The outcome is labelled as adjusted for covariates unless only ccc = "none" is shown
    @test layout.ylabel == "y (covariate-adjusted)" && DiDInt.plot_layout(x; plot_type = "event", ccc = "none").ylabel == "y"

    # Selected treated states: the means are those of these states. Under "none" and "int" the
    # lambdas of a state do not depend on the other states, so the result is that of
    # didint_plot() for these states alone.
    alone = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = ["s01", "s03"], treatment_times = [2004, 2006],
                           covariates = INTEGRATION_COVARIATES, plot_type = "event", ccc = ["none", "int"])
    end
    selected = DiDInt.plot_layout(x; plot_type = "event", treated_states = ["s03", "s01"], ccc = ["int", "none"], groupmin = 2)
    @test selected.panels == [DiDInt.PLOT_CCC_LABELS["none"], DiDInt.PLOT_CCC_LABELS["int"]]
    for (panel, option) in enumerate(["none", "int"])
        rows = alone.event[alone.event.ccc .== option, :]
        series = panel_series(selected, panel)
        @test series.x == rows.time_since_treatment
        @test series.y ≈ rows.y atol = 1e-10
        @test collect(skipmissing(series.lower)) ≈ collect(skipmissing(rows.ci_lower)) atol = 1e-10
        @test collect(skipmissing(series.upper)) ≈ collect(skipmissing(rows.ci_upper)) atol = 1e-10
    end

    # weights = false: the treated states count equally, whatever was chosen in didint_plot()
    none = x.event_lambdas[x.event_lambdas.ccc .== "none", :]
    unweighted = DiDInt.plot_layout(x; plot_type = "event", weights = false, ccc = "none")
    expected = combine(groupby(none, :time_since_treatment), :lambda => mean => :y)
    sort!(expected, :time_since_treatment)
    @test unweighted.series.x == expected.time_since_treatment
    @test unweighted.series.y ≈ expected.y atol = 1e-12
    @test unweighted.series.y ≈ plot_object(plot_type = "event", weights = false, ccc = "none").event.y atol = 1e-12

    # A single treated state (no intervals) and a window of periods; the legend only explains
    # what is in the window
    one = DiDInt.plot_layout(x; plot_type = "event", treated_states = "s05", window = (-2, 1))
    @test sort(unique(one.series.x)) == -2:1 && all(ismissing, one.series.lower)
    @test first.(DiDInt.plot_layout(x; plot_type = "event", window = (-1, 1)).legend_extras) == [:band]

    # What cannot be drawn
    @test_throws r"'plot_type' must be" DiDInt.plot_layout(x; plot_type = "trends")
    @test_throws r"plot_type = \"event\" were not prepared" DiDInt.plot_layout(plot_object(plot_type = "conditional"); plot_type = "event")
    @test_throws r"Not among the treated states" DiDInt.plot_layout(x; plot_type = "event", treated_states = "s09")
    @test_throws r"No data to plot for ccc = \"time\"" DiDInt.plot_layout(plot_object(ccc = "hom"); plot_type = "event", ccc = "time")
    @test_throws r"'window' must be" DiDInt.plot_layout(x; plot_type = "event", window = (3, 1))
    @test_throws r"No data within 'window'" DiDInt.plot_layout(x; plot_type = "event", window = (50, 60))
end

@testset "Conditional parallel trends plot, one panel per ccc option" begin
    x = plot_object()
    options = ["none", "hom", "time", "state", "add"]

    # The data: under "none" the lambdas are the cell means
    none = x.conditional[x.conditional.ccc .== "none", :]
    @test nrow(x.conditional) == 5 * 12 * 8
    @test all(isapprox(row.lambda, mean(cell_rows(INTEGRATION_PANEL, row.state, 2001 + row.period).y); atol = 1e-10) for row in eachrow(none))
    @test all(row.n == nrow(cell_rows(INTEGRATION_PANEL, row.state, 2001 + row.period)) for row in eachrow(x.conditional))

    # The default plot
    layout = DiDInt.plot_layout(x)
    @test layout.plot_type == "conditional" && layout.legend_title == "Treated state"
    @test layout.title == "Conditional parallel trends: y means (covariate-adjusted) by period" && layout.subtitle == ""
    @test layout.ylabel == "y (covariate-adjusted)"
    @test layout.panels == [DiDInt.PLOT_CCC_LABELS[option] for option in options]
    @test layout.xticks == (collect(0:7), string.(2001:2008))
    unadjusted = DiDInt.plot_layout(x; ccc = "none")
    @test unadjusted.title == "Conditional parallel trends: y means by period" && unadjusted.ylabel == "y"

    # A line for every state: the treated states each in their colour, the others as controls
    lines = unique(layout.series[!, [:panel, :group, :label, :role, :color]])
    @test nrow(lines) == 5 * 12
    for (panel, option) in enumerate(options)
        lambda, _ = conditional_values(x, option)
        for (index, s) in enumerate(INTEGRATION_TREATED)
            series = panel_series(layout, panel; group = s)
            @test all(series.role .== :treated) && all(series.label .== s) && all(series.color .== index)
            @test series.x == 0:7 && series.y == [lambda[(s, period)] for period in 0:7]
        end
        for c in PLOT_NEVER_TREATED
            series = panel_series(layout, panel; group = c)
            @test all(series.role .== :control) && all(series.label .== "Control states") && all(series.color .== 0)
            @test series.y == [lambda[(c, period)] for period in 0:7]
        end
        @test layout.vlines.x[layout.vlines.panel .== panel] == [3, 5, 6]
    end

    # Selected treated states keep their colours; the other treated states are left out
    selected = DiDInt.plot_layout(x; treated_states = ["s03"], ccc = "hom")
    @test selected.panels == [DiDInt.PLOT_CCC_LABELS["hom"]]
    @test sort(unique(selected.series.group)) == sort(vcat("s03", PLOT_NEVER_TREATED))
    @test all(panel_series(selected, 1; group = "s03").color .== 3) && selected.vlines.x == [5]
    @test DiDInt.plot_layout(x; ccc = "all").panels == layout.panels
    @test DiDInt.plot_layout(x; window = (2, 5)).xticks == (collect(2:5), string.(2003:2006))
    @test isempty(DiDInt.plot_layout(x; treated_states = "s05", window = (0, 3)).legend_extras)

    # Every treated state has its own colour and entry in the legend, however many there are,
    # and keeps its colour when only some treated states are selected
    many = quiet_call() do
        DiDInt.didint_plot("coll", "state", "year", TEST_DATA_FULL; treated_states = ["88", "57", "58", "59", "61", "64", "71", "72", "85", "34"],
                           treatment_times = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000], ccc = "none",
                           plot_type = ["conditional", "conditional2"])
    end
    merit_states = ["34", "57", "58", "59", "61", "64", "71", "72", "85", "88"]
    all_states = DiDInt.plot_layout(many)
    @test unique(all_states.series[all_states.series.role .== :treated, [:group, :label, :color]]) ==
          DataFrame(group = merit_states, label = merit_states, color = 1:10)
    @test all_states.legend_title == "Treated state"
    # a state can be selected by its number when the states are numbers written as strings
    two = DiDInt.plot_layout(many; treated_states = [71, 58])
    @test unique(two.series[two.series.role .== :treated, [:group, :color]]) == DataFrame(group = ["58", "71"], color = [3, 7])

    # Panels of treated states run from the earliest treated state to the latest
    @test DiDInt.plot_layout(many; facet = "state").panels ==
          ["71 (treated 1991)", "58 (treated 1993)", "64 (treated 1996)", "59 (treated 1997)", "85 (treated 1997)",
           "57 (treated 1998)", "72 (treated 1998)", "61 (treated 1999)", "34 (treated 2000)", "88 (treated 2000)"]
    @test DiDInt.plot_layout(many; plot_type = "conditional2", treated_states = [34, 71, 64]).panels ==
          ["71 (treated 1991)", "64 (treated 1996)", "34 (treated 2000)"]

    # Without the treated states every state is drawn alike, and all treatment periods are marked
    unknown = DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES, covariates = INTEGRATION_COVARIATES,
                                 ccc = ["hom", "none"], plot_type = "conditional")
    alike = DiDInt.plot_layout(unknown)
    @test all(alike.series.role .== :state) && length(unique(alike.series.group)) == 12 && alike.vlines.x == [3, 5, 6, 3, 5, 6]
    @test_throws r"were not given to didint_plot" DiDInt.plot_layout(unknown; treated_states = "s01")
    @test_throws r"needs the treated states" DiDInt.plot_layout(unknown; facet = "state", ccc = "hom")
    @test_throws r"plot_type = \"conditional2\" were not prepared" DiDInt.plot_layout(unknown; plot_type = "conditional2")

    # A state treated throughout is neither a treated state nor a control state
    early = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = vcat(INTEGRATION_TREATED, "s06"),
                           treatment_times = vcat(INPUT_TIMES, 2001), ccc = "none")
    end
    @test !("s06" in DiDInt.plot_layout(early).series.group)
    @test !("s06" in DiDInt.plot_layout(early; facet = "state", controls = "separate").series.group)
    @test !("s06" in early.conditional2.state) && !("s06" in early.event_lambdas.state)

    @test_throws r"plot_type = \"conditional\" were not prepared" DiDInt.plot_layout(plot_object(ccc = "int"); plot_type = "conditional")
    @test_throws r"'facet' must be" DiDInt.plot_layout(x; facet = "cohort")
    @test_throws r"'controls' must be" DiDInt.plot_layout(x; controls = "mean")
end

@testset "Conditional parallel trends plot, one panel per treated state" begin
    x = plot_object(ccc = ["hom", "int"])
    lambda, n = conditional_values(x, "hom")

    # The treated state next to the average control state. The average starts at the controls'
    # weighted mean and then adds the weighted mean of their changes, with weights n_t + n_{t-1}
    layout = DiDInt.plot_layout(x; plot_type = "conditional", facet = "state")
    @test layout.panels == PLOT_STATE_TITLES && layout.subtitle == "ccc = \"hom\""
    @test layout.title == "Conditional parallel trends: y means (covariate-adjusted) by period"
    @test layout.vlines.panel == 1:5 && layout.vlines.x == PLOT_TREATMENT_PERIODS
    for (panel, s) in enumerate(INTEGRATION_TREATED)
        own = panel_series(layout, panel; role = :treated)
        @test all(own.group .== s) && own.x == 0:7 && own.y == [lambda[(s, period)] for period in 0:7]
        average = panel_series(layout, panel; role = :control_mean)
        @test average.x == 0:7 && all(average.label .== "Control states (mean, weighted by observations)")
        @test average.y[1] ≈ sum(n[(c, 0)] * lambda[(c, 0)] for c in PLOT_NEVER_TREATED) / sum(n[(c, 0)] for c in PLOT_NEVER_TREATED) atol = 1e-12
        for period in 1:7
            weight(c) = n[(c, period)] + n[(c, period - 1)]
            change = sum(weight(c) * (lambda[(c, period)] - lambda[(c, period - 1)]) for c in PLOT_NEVER_TREATED) / sum(weight, PLOT_NEVER_TREATED)
            @test average.y[period + 1] - average.y[period] ≈ change atol = 1e-12
        end
    end

    # With equal weights and the same controls in every period, it is the mean of the controls' lambdas
    equal = DiDInt.plot_layout(x; plot_type = "conditional", facet = "state", weights = false, treated_states = "s04")
    @test equal.panels == ["s04 (treated 2006)"]
    average = panel_series(equal, 1; role = :control_mean)
    @test all(average.label .== "Control states (mean)")
    @test average.y ≈ [mean(lambda[(c, period)] for c in PLOT_NEVER_TREATED) for period in 0:7] atol = 1e-12

    # One line per control state instead
    separate = DiDInt.plot_layout(x; plot_type = "conditional", facet = "state", controls = "separate")
    for panel in 1:5
        controls = panel_series(separate, panel; role = :control)
        @test sort(unique(controls.group)) == PLOT_NEVER_TREATED
        @test all(row.y == lambda[(row.group, row.x)] for row in eachrow(controls))
        @test nrow(panel_series(separate, panel; role = :treated)) == 8
    end

    # notyet = true: a state of a later cohort is a control up to the period before its
    # treatment, as in didint(). For s01 (treated in period 3): s03 and s04 until period 4, s05
    # until period 5, and never s02 (the same cohort).
    with_notyet = plot_object(ccc = "hom", notyet = true)
    lambda, n = conditional_values(with_notyet, "hom")
    separate = DiDInt.plot_layout(with_notyet; facet = "state", controls = "separate", treated_states = "s01")
    @test panel_series(separate, 1; group = "s03", role = :control).x == 0:4
    @test panel_series(separate, 1; group = "s05", role = :control).x == 0:5
    @test panel_series(separate, 1; group = "s07", role = :control).x == 0:7
    @test nrow(panel_series(separate, 1; group = "s02", role = :control)) == 0
    controls_at(period) = plot_controls("s01", period, true)
    average = panel_series(DiDInt.plot_layout(with_notyet; facet = "state", weights = false, treated_states = "s01"), 1; role = :control_mean)
    @test average.y[1] ≈ mean(lambda[(c, 0)] for c in controls_at(0)) atol = 1e-12
    @test diff(average.y) ≈ [mean(lambda[(c, period)] - lambda[(c, period - 1)] for c in controls_at(period)) for period in 1:7] atol = 1e-12

    # The average only depends on the controls' changes: shifting the lambdas of each state by a
    # constant moves the whole line by a constant, also when the set of controls changes
    shift = Dict(s => 10.0 * i for (i, s) in enumerate(sort(unique(INTEGRATION_PANEL.state))))
    shifted = Dict(key => value + shift[key[1]] for (key, value) in lambda)
    for weights in [true, false]
        periods, values = DiDInt.average_of_controls(lambda, n, controls_at, 0:7, weights)
        _, shifted_values = DiDInt.average_of_controls(shifted, n, controls_at, 0:7, weights)
        @test periods == 0:7
        @test diff(shifted_values) ≈ diff(values) atol = 1e-10
    end

    # Where no control has a change, the line starts again from the controls' mean; periods
    # without any control have no value
    gaps = Dict(("a", 0) => 1.0, ("a", 1) => 2.0, ("b", 0) => 3.0, ("b", 1) => 3.0, ("a", 3) => 7.0, ("b", 3) => 9.0, ("a", 4) => 8.0)
    periods, values = DiDInt.average_of_controls(gaps, Dict(key => 1 for key in keys(gaps)), period -> ["a", "b"], 0:4, true)
    @test periods == [0, 1, 3, 4] && values ≈ [2.0, 2.5, 8.0, 9.0]

    # Several ccc options: a row of panels per treated state, a column per option, each panel
    # as the plot of its option alone
    y = plot_object()
    @test DiDInt.plot_layout(y; facet = "state", ccc = "time").subtitle == "ccc = \"time\""
    options = DiDInt.selected_ccc(unique(y.conditional.ccc), nothing)
    @test length(options) >= 2
    grid = DiDInt.plot_layout(y; facet = "state")
    @test grid.subtitle == "" && grid.ncols == length(options) && grid.title == layout.title
    @test grid.panels == [state * ", ccc = \"" * option * "\"" for state in PLOT_STATE_TITLES for option in options]
    for (i, option) in enumerate(options), controls in ["average", "separate"]
        alone = DiDInt.plot_layout(y; facet = "state", ccc = option, controls = controls)
        grid = DiDInt.plot_layout(y; facet = "state", controls = controls)
        @test isnothing(alone.ncols)
        for state in eachindex(PLOT_STATE_TITLES)
            grid_panel = grid.series[grid.series.panel .== (state - 1) * length(options) + i, Not(:panel)]
            @test isequal(grid_panel, alone.series[alone.series.panel .== state, Not(:panel)])
        end
    end
    two = DiDInt.plot_layout(y; facet = "state", ccc = ["state", "hom"], treated_states = "s05")
    @test two.panels == ["s05 (treated 2007), ccc = \"hom\"", "s05 (treated 2007), ccc = \"state\""] && two.ncols == 2
    @test two.vlines == DataFrame(panel = [1, 2], x = [6, 6])
end

@testset "Plots of changes since the pre-treatment period (conditional2 and integrated2)" begin
    # The treated state's change minus the average of its controls' changes is ATT(s,t) of
    # didint() in every treated period: with the controls weighted by their numbers of
    # observations under weighting = "both" and equally under "none", with and without the
    # not-yet-treated states as controls
    for notyet in [false, true], (weights, weighting) in [(true, "both"), (false, "none")]
        x = plot_object(notyet = notyet, ccc = ["hom", "add", "int"], plot_type = ["conditional2", "integrated2"])
        for option in ["hom", "add"]
            layout = DiDInt.plot_layout(x; plot_type = "conditional2", ccc = option, weights = weights)
            @test layout.plot_type == "conditional2" && layout.panels == PLOT_STATE_TITLES
            @test layout.vlines.x == PLOT_TREATMENT_PERIODS && layout.xticks == (collect(0:7), string.(2001:2008))
            test_gaps_are_atts(layout, integration_didint(ccc = option, agg = "sgt", weighting = weighting, notyet = notyet))
        end
        layout = DiDInt.plot_layout(x; plot_type = "integrated2", weights = weights)
        @test layout.plot_type == "integrated2" && layout.panels == PLOT_STATE_TITLES
        @test layout.vlines.x == PLOT_TREATMENT_PERIODS && layout.xticks == (collect(0:7), string.(2001:2008))
        test_gaps_are_atts(layout, integration_didint(ccc = "int", agg = "sgt", weighting = weighting, notyet = notyet))
    end

    # Titles and labels
    x = plot_object()
    conditional2 = DiDInt.plot_layout(x; plot_type = "conditional2", ccc = "state")
    @test conditional2.title == "Conditional parallel trends: changes since the pre-treatment period"
    @test conditional2.subtitle == "ccc = \"state\"" && conditional2.ylabel == "Change in y (covariate-adjusted)"
    @test DiDInt.plot_layout(x; plot_type = "conditional2", ccc = "none").ylabel == "Change in y"
    integrated2 = DiDInt.plot_layout(x; plot_type = "integrated2")
    @test integrated2.title == "Integrated parallel trends: changes since the pre-treatment period"
    @test integrated2.subtitle == "" && integrated2.ylabel == "Change in predicted y"
    # a title of one's own replaces both lines
    custom = DiDInt.plot_layout(x; plot_type = "conditional2", ccc = "state", title = "My title")
    @test custom.title == "My title" && custom.subtitle == ""
    @test DiDInt.plot_layout(x; plot_type = "integrated2", title = "").title == ""

    # One line per control state, each through zero in the pre-treatment period
    for (plot_type, table, column, arguments) in [("conditional2", x.conditional2[x.conditional2.ccc .== "hom", :], :treated_state, (ccc = "hom",)),
                                                  ("integrated2", x.integrated2, :treated_state_covars, (;))]
        separate = DiDInt.plot_layout(x; plot_type = plot_type, controls = "separate", treated_states = "s05", arguments...)
        @test separate.panels == ["s05 (treated 2007)"]
        @test sort(unique(separate.series.group)) == sort(vcat("s05", PLOT_NEVER_TREATED))
        @test all(separate.series.y[separate.series.x .== 5] .== 0) && count(separate.series.x .== 5) == 1 + length(PLOT_NEVER_TREATED)
        delta = Dict((row.state, row.period) => row.delta for row in eachrow(table[table[!, column] .== "s05", :]))
        @test all(row.y == delta[(row.group, row.x)] for row in eachrow(separate.series) if row.x != 5)
        @test all((row.role == :treated) == (row.group == "s05") for row in eachrow(separate.series))
    end

    # Several ccc options: a row of panels per treated state, a column per option, each panel
    # as the plot of its option alone
    options = DiDInt.selected_ccc(unique(x.conditional2.ccc), nothing)
    @test length(options) >= 2
    grid = DiDInt.plot_layout(x; plot_type = "conditional2")
    @test grid.subtitle == "" && grid.ncols == length(options)
    @test grid.panels == [state * ", ccc = \"" * option * "\"" for state in PLOT_STATE_TITLES for option in options]
    for (i, option) in enumerate(options)
        alone = DiDInt.plot_layout(x; plot_type = "conditional2", ccc = option)
        @test isnothing(alone.ncols)
        for state in eachindex(PLOT_STATE_TITLES)
            grid_panel = grid.series[grid.series.panel .== (state - 1) * length(options) + i, Not(:panel)]
            @test isequal(grid_panel, alone.series[alone.series.panel .== state, Not(:panel)])
        end
    end
    two = DiDInt.plot_layout(x; plot_type = "conditional2", ccc = ["state", "hom"], treated_states = "s05")
    @test two.panels == ["s05 (treated 2007), ccc = \"hom\"", "s05 (treated 2007), ccc = \"state\""] && two.ncols == 2

    # What cannot be drawn
    @test_throws r"cannot select another option" DiDInt.plot_layout(x; plot_type = "integrated2", ccc = "hom")
    @test_throws r"cannot select another option" DiDInt.plot_layout(x; plot_type = "integrated", ccc = ["int", "hom"])
    for ccc in ["int", "all", "ALL"], plot_type in ["integrated", "integrated2"]
        @test isequal(DiDInt.plot_layout(x; plot_type = plot_type, ccc = ccc).series, DiDInt.plot_layout(x; plot_type = plot_type).series)
    end
    @test_throws r"plot_type = \"integrated2\" were not prepared" DiDInt.plot_layout(plot_object(ccc = "hom"); plot_type = "integrated2")
end

@testset "Integrated parallel trends plot, changes since the previous period" begin
    x = plot_object()

    # In its first treated period, a treated state's change minus the average of its controls'
    # changes is its ATT(s,g) under ccc = "int": with the controls weighted by their numbers
    # of observations under weighting = "both", and equally under weighting = "none"
    for (weights, weighting) in [(true, "both"), (false, "none")]
        layout = DiDInt.plot_layout(x; plot_type = "integrated", weights = weights)
        @test layout.plot_type == "integrated" && layout.panels == PLOT_STATE_TITLES
        @test layout.title == "Integrated parallel trends: changes since the previous period" && layout.ylabel == "Change in predicted y"
        @test layout.vlines.x == PLOT_TREATMENT_PERIODS && layout.xticks == (collect(1:7), string.(2002:2008))
        sgt = integration_didint(ccc = "int", agg = "sgt", weighting = weighting)
        for (panel, s) in enumerate(INTEGRATION_TREATED)
            g = PLOT_TREATMENT_PERIODS[panel]
            own = panel_series(layout, panel; role = :treated)
            average = panel_series(layout, panel; role = :control_mean)
            @test own.x == 1:7 && average.x == 1:7 && all(own.group .== s)
            att = only(sgt.att_sgt[(sgt.state .== s) .& (year.(sgt.t) .== INPUT_TIMES[panel])])
            @test only(own.y[own.x .== g]) - only(average.y[average.x .== g]) ≈ att atol = 1e-8
        end
    end

    # One line per control state instead, with the deltas of the data
    separate = DiDInt.plot_layout(x; plot_type = "integrated", controls = "separate", treated_states = ["s05", "s02"])
    @test separate.panels == ["s02 (treated 2004)", "s05 (treated 2007)"]
    for (panel, s) in enumerate(["s02", "s05"])
        delta = Dict((row.state, row.period) => row.delta for row in eachrow(x.integrated[x.integrated.treated_state_covars .== s, :]))
        series = panel_series(separate, panel)
        @test sort(unique(series.group)) == sort(vcat(s, PLOT_NEVER_TREATED))
        @test all(row.y == delta[(row.group, row.x)] for row in eachrow(series))
        @test all((row.role == :treated) == (row.group == s) for row in eachrow(series))
    end

    # notyet = true adds the not-yet-treated states of later cohorts to the controls of the average
    with_notyet = plot_object(ccc = "int", notyet = true, plot_type = "integrated")
    layout = DiDInt.plot_layout(with_notyet; weights = false, treated_states = "s01")
    @test layout.plot_type == "integrated"
    others = with_notyet.integrated[(with_notyet.integrated.treated_state_covars .== "s01") .& (with_notyet.integrated.state .!= "s01"), :]
    @test panel_series(layout, 1; role = :control_mean).y ≈ [mean(others.delta[others.period .== period]) for period in 1:7] atol = 1e-12
    @test count(others.period .== 1) == length(PLOT_NEVER_TREATED) + 3 && count(others.period .== 7) == length(PLOT_NEVER_TREATED)

    @test DiDInt.plot_layout(x; plot_type = "integrated", ccc = "int").panels == PLOT_STATE_TITLES
    @test_throws r"cannot select another option" DiDInt.plot_layout(x; plot_type = "integrated", ccc = "hom")
    @test_throws r"plot_type = \"integrated\" were not prepared" DiDInt.plot_layout(plot_object(ccc = "hom"); plot_type = "integrated")
end

@testset "plot() draws every kind of plot with Makie" begin
    @test !isnothing(Base.get_extension(DiDInt, :DiDIntMakieExt))
    x = plot_object(notyet = true)
    blocks(figure, type) = [block for block in figure.content if block isa type]

    cases = [((;), 5),
             ((plot_type = "event",), 6),
             ((plot_type = "event", ccc = ["none", "int"], treated_states = ["s01", "s03", "s05"], window = (-3, 2), groupmin = 2), 2),
             # a confidence interval in a single period (drawn as a bar)
             ((plot_type = "event", ccc = "hom", window = (1, 4), groupmin = 5), 1),
             ((plot_type = "conditional", treated_states = "s03", ccc = "hom"), 1),
             ((facet = "state", ccc = "add"), 5),
             ((facet = "state", ccc = ["hom", "add"], treated_states = ["s01", "s04", "s05"]), 6),
             ((facet = "state", ccc = "state", controls = "separate", weights = false), 5),
             ((plot_type = "conditional2", ccc = "hom"), 5),
             ((plot_type = "conditional2", ccc = "time", controls = "separate", treated_states = "s04"), 1),
             ((plot_type = "conditional2", treated_states = ["s01", "s04"]), 10),
             ((plot_type = "integrated",), 5),
             ((plot_type = "integrated", controls = "separate", treated_states = ["s01", "s05"]), 2),
             ((plot_type = "integrated2",), 5),
             ((plot_type = "integrated2", controls = "separate", weights = false, treated_states = ["s03", "s02"]), 2)]
    for (arguments, n_panels) in cases
        figure = CairoMakie.plot(x; arguments...)
        @test figure isa CairoMakie.Figure
        @test length(blocks(figure, CairoMakie.Axis)) == n_panels && length(blocks(figure, CairoMakie.Legend)) == 1
        file = tempname() * ".png"
        CairoMakie.save(file, figure)
        @test filesize(file) > 5_000
        rm(file; force = true)
    end

    # The title is drawn with the ccc option on a line under it; a title of one's own, or none, replaces both
    texts(figure) = [block.text[] for block in blocks(figure, CairoMakie.Label)]
    drawn = texts(CairoMakie.plot(x; plot_type = "conditional2", ccc = "hom"))
    @test "Conditional parallel trends: changes since the pre-treatment period" in drawn && "ccc = \"hom\"" in drawn
    drawn = texts(CairoMakie.plot(x; plot_type = "conditional2", ccc = "hom", title = "Mine"))
    @test "Mine" in drawn && !("ccc = \"hom\"" in drawn)
    @test !any(startswith("Integrated parallel trends"), texts(CairoMakie.plot(x; plot_type = "integrated2", title = "")))

    # Appearance: number of columns, size and theme
    figure = CairoMakie.plot(x; plot_type = "integrated", ncols = 1, size = (500, 1100), title = "", theme = CairoMakie.Theme())
    @test figure isa CairoMakie.Figure && length(blocks(figure, CairoMakie.Axis)) == 5
    @test Tuple(figure.scene.viewport[].widths) == (500, 1100)
    @test_throws r"'ncols' must be" CairoMakie.plot(x; ncols = 0)
    @test_throws r"were not prepared" CairoMakie.plot(plot_object(plot_type = "conditional"); plot_type = "event")

    # Every coloured line can be told from the others, also when there are more of them than colours
    extension = Base.get_extension(DiDInt, :DiDIntMakieExt)
    @test allunique((extension.series_color(index), extension.series_linestyle(index)) for index in 1:48)
    @test length(unique(extension.series_color.(1:16))) == 16 && all(extension.series_linestyle.(1:16) .== :solid)
end
