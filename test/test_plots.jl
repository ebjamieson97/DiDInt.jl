# Tests of DiDIntPlot, the object that didint_plot() returns, and of its plots: what each kind of
# plot shows (the layouts of src/plots.jl) and that Makie draws them.
# The panel is INTEGRATION_PANEL (see test_integration.jl): 12 states observed 2001-2008 (periods 0
# to 7), with s01, s02 treated in 2004, s03, s04 in 2006 and s05 in 2007; s06 to s12 are never treated.

using Dates
using Statistics
import CairoMakie

plot_object(; kwargs...) = quiet_call() do
    DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                       covariates = INTEGRATION_COVARIATES, kwargs...)
end
const PLOT_NEVER_TREATED = setdiff(sort(unique(INTEGRATION_PANEL.state)), INTEGRATION_TREATED)
const PLOT_TREATMENT_PERIODS = INPUT_TIMES .- 2001
# Titles of the panels of the treated states
const PLOT_STATE_TITLES = [string(s, " (treated ", g, ")") for (s, g) in zip(INTEGRATION_TREATED, INPUT_TIMES)]

# The points of one panel of a layout, optionally of one line or of one kind of line
function panel_series(layout, panel; group = nothing, role = nothing)
    keep = layout.series.panel .== panel
    isnothing(group) || (keep = keep .& (layout.series.group .== group))
    isnothing(role) || (keep = keep .& (layout.series.role .== role))
    return layout.series[keep, :]
end

# The lambdas and numbers of observations of the conditional rows of one ccc option, by (state, period)
function conditional_values(plot, ccc)
    rows = plot.data[isequal.(plot.data.ccc, ccc) .& .!ismissing.(plot.data.lambda), :]
    return Dict(zip(zip(rows.state, rows.period), rows.lambda)), Dict(zip(zip(rows.state, rows.period), rows.n))
end

@testset "didint_plot() returns a DiDIntPlot" begin
    parallel = plot_object()
    @test parallel isa DiDIntPlot
    @test isequal(DataFrame(parallel), parallel.data)
    @test DiDInt.Tables.istable(parallel) && DiDInt.Tables.schema(parallel).names == Tuple(Symbol.(names(parallel.data)))
    specs = parallel.specs
    @test nrow(parallel.lambdas) == 0 && nrow(parallel.long_differences) > 0
    @test DiDInt.available_plot_types(parallel) == ["conditional", "conditional2", "integrated", "integrated2"]
    @test DiDInt.available_plot_types(plot_object(ccc = "hom")) == ["conditional", "conditional2"]
    @test nrow(plot_object(ccc = "hom").long_differences) == 0
    @test specs.period_labels == string.(2001:2008)

    # The specification that the plots use
    @test specs.outcome == "y" && !specs.event && !specs.notyet && specs.covariates == INTEGRATION_COVARIATES
    @test specs.treated_states == INTEGRATION_TREATED && specs.treatment_periods == PLOT_TREATMENT_PERIODS
    @test specs.treatment_times == string.(INPUT_TIMES) && isempty(specs.always_treated)
    @test specs.period == "1 year" && specs.start_date == "2001" && specs.end_date == "2008"
    @test plot_object(notyet = true).specs.notyet

    # Event study data come with the lambdas of the treated states
    event = plot_object(event = true)
    @test event isa DiDIntPlot && event.specs.event && DiDInt.available_plot_types(event) == ["event"]
    @test nrow(event.long_differences) == 0
    @test names(event.lambdas) == ["ccc", "state", "time_since_treatment", "lambda", "n"]
    @test sort(unique(event.lambdas.state)) == INTEGRATION_TREATED

    # The Stata and R wrappers receive the table itself
    for wrapper in ["r", "stata"]
        @test isequal(plot_object(wrapper = wrapper), parallel.data)
        @test isequal(plot_object(wrapper = wrapper, event = true), event.data)
    end
    @test plot_object(wrapper = "julia") isa DiDIntPlot

    # Without the treated states, the specification has the distinct treatment periods
    unknown = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES, ccc = ["hom", "none"])
    end
    @test isempty(unknown.specs.treated_states) && unknown.specs.treatment_periods == [3, 5, 6]
    @test unknown.specs.treatment_times == ["2004", "2006", "2007"]

    # A state treated at or before the first date is remembered as such
    early = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = vcat(INTEGRATION_TREATED, "s06"),
                           treatment_times = vcat(INPUT_TIMES, 2001), ccc = "none")
    end
    @test early.specs.treated_states == INTEGRATION_TREATED && early.specs.always_treated == ["s06"]

    # Display
    text = sprint(show, MIME("text/plain"), parallel)
    @test occursin("DID-INT plot data: parallel trends", text) && occursin("s01 (2004)", text)
    @test occursin("\"conditional\", \"conditional2\", \"integrated\", \"integrated2\"", text)
    @test occursin("(drawn with plot", text) && occursin("never-treated states", text)
    @test all(length(line) <= 90 for line in split(first(split(text, "\n\nData")), "\n"))
    @test DiDInt.wrapped_list(["aa", "bb", "cc"], 2; width = 6) == "aa, bb,\n  cc"
    @test occursin("DID-INT plot data: event study", sprint(show, MIME("text/plain"), event))
    @test occursin("not specified (treatment times: 2004, 2006, 2007)", sprint(show, MIME("text/plain"), unknown))
    @test sprint(show, parallel) == "DiDIntPlot(parallel trends plot data, $(nrow(parallel.data)) rows)"
end

@testset "Event study plot" begin
    event = plot_object(event = true)
    options = ["none", "hom", "time", "state", "add", "int"]

    # The data are the means of the treated states' lambdas, weighted by their numbers of observations
    for group in groupby(event.lambdas, [:ccc, :time_since_treatment])
        row = only(eachrow(event.data[(event.data.ccc .== first(group.ccc)) .& (event.data.time_since_treatment .== first(group.time_since_treatment)), :]))
        @test row.y ≈ sum(group.n .* group.lambda) / sum(group.n) atol = 1e-12
        @test row.ngroup == nrow(group)
    end
    none = event.lambdas[event.lambdas.ccc .== "none", :]
    @test all(isapprox(row.lambda, mean(cell_rows(INTEGRATION_PANEL, row.state, INTEGRATION_COHORTS[row.state] + row.time_since_treatment).y); atol = 1e-12)
              for row in eachrow(none))

    # One panel per ccc option with the estimates of the data; an interval needs 3 treated states
    layout = DiDInt.plot_layout(event)
    @test layout.plot_type == "event" && layout.title == "Event study" && isnothing(layout.xticks)
    @test layout.panels == [DiDInt.PLOT_CCC_LABELS[option] for option in options]
    @test layout.vlines.panel == 1:6 && all(layout.vlines.x .== 0)
    for (panel, option) in enumerate(options)
        rows = event.data[event.data.ccc .== option, :]
        series = panel_series(layout, panel)
        @test all(series.role .== :estimate) && series.x == rows.time_since_treatment && series.y == rows.y
        @test all(ismissing(lower) ? ngroup < 3 : lower == ci for (lower, ci, ngroup) in zip(series.lower, rows.ci_lower, rows.ngroup))
        @test all(ismissing(upper) ? ngroup < 3 : upper == ci for (upper, ci, ngroup) in zip(series.upper, rows.ci_upper, rows.ngroup))
    end
    @test any(ismissing, layout.series.lower) && any(!ismissing, layout.series.lower)
    @test first.(layout.legend_extras) == [:band, :no_interval] && layout.legend_extras[1][2] == "95% confidence interval"
    with_two = DiDInt.plot_layout(event; groupmin = 2)
    @test isequal(with_two.series.lower, event.data.ci_lower[sortperm([findfirst(==(c), options) for c in event.data.ccc])])
    @test plot_object(event = true, ci = 0.9).data.ci_lower != event.data.ci_lower
    @test DiDInt.plot_layout(plot_object(event = true, ci = 0.9)).legend_extras[1][2] == "90% confidence interval"

    # Selected treated states: the means are those of these states. Under "none" and "int" the
    # lambdas of a state do not depend on the other states, so the result is that of
    # didint_plot() for these states alone.
    alone = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = ["s01", "s03"], treatment_times = [2004, 2006],
                           covariates = INTEGRATION_COVARIATES, event = true, ccc = ["none", "int"])
    end
    selected = DiDInt.plot_layout(event; treated_states = ["s03", "s01"], ccc = ["int", "none"], groupmin = 2)
    @test selected.panels == [DiDInt.PLOT_CCC_LABELS["none"], DiDInt.PLOT_CCC_LABELS["int"]]
    for (panel, option) in enumerate(["none", "int"])
        rows = alone.data[alone.data.ccc .== option, :]
        series = panel_series(selected, panel)
        @test series.x == rows.time_since_treatment
        @test series.y ≈ rows.y atol = 1e-10
        @test collect(skipmissing(series.lower)) ≈ collect(skipmissing(rows.ci_lower)) atol = 1e-10
        @test collect(skipmissing(series.upper)) ≈ collect(skipmissing(rows.ci_upper)) atol = 1e-10
    end

    # weights = false: the treated states count equally, whatever was chosen in didint_plot()
    unweighted = DiDInt.plot_layout(event; weights = false, ccc = "none")
    expected = combine(groupby(none, :time_since_treatment), :lambda => mean => :y)
    sort!(expected, :time_since_treatment)
    @test unweighted.series.x == expected.time_since_treatment
    @test unweighted.series.y ≈ expected.y atol = 1e-12
    @test unweighted.series.y ≈ plot_object(event = true, weights = false, ccc = "none").data.y atol = 1e-12

    # A single treated state (no intervals) and a window of periods
    one = DiDInt.plot_layout(event; treated_states = "s05", window = (-2, 1))
    @test sort(unique(one.series.x)) == -2:1 && all(ismissing, one.series.lower)
    # the legend only explains what is in the window
    @test first.(DiDInt.plot_layout(event; window = (-1, 1)).legend_extras) == [:band]
    @test isempty(DiDInt.plot_layout(plot_object(); treated_states = "s05", window = (0, 3)).legend_extras)

    # The outcome is labelled as adjusted for covariates unless only ccc = "none" is shown
    @test layout.ylabel == "y (adjusted for covariates)" && DiDInt.plot_layout(event; ccc = "none").ylabel == "y"

    # What cannot be drawn
    @test_throws r"event = false" DiDInt.plot_layout(event; plot_type = "integrated")
    @test_throws r"event = true" DiDInt.plot_layout(plot_object(); plot_type = "event")
    @test_throws r"'plot_type' must be" DiDInt.plot_layout(event; plot_type = "trends")
    @test_throws r"Not among the treated states" DiDInt.plot_layout(event; treated_states = "s09")
    @test_throws r"No data to plot for ccc = \"time\"" DiDInt.plot_layout(plot_object(event = true, ccc = "hom"); ccc = "time")
    @test_throws r"'window' must be" DiDInt.plot_layout(event; window = (3, 1))
    @test_throws r"No data within 'window'" DiDInt.plot_layout(event; window = (50, 60))
end

@testset "Conditional parallel trends plot, one panel per ccc option" begin
    parallel = plot_object()
    options = ["none", "hom", "time", "state", "add"]
    layout = DiDInt.plot_layout(parallel)
    @test layout.plot_type == "conditional" && layout.legend_title == "Treated state"
    @test layout.panels == [DiDInt.PLOT_CCC_LABELS[option] for option in options]
    @test layout.xticks == (collect(0:7), string.(2001:2008))

    # A line for every state: the treated states each in their colour, the others as controls
    lines = unique(layout.series[!, [:panel, :group, :label, :role, :color]])
    @test nrow(lines) == 5 * 12
    for (panel, option) in enumerate(options)
        lambda, _ = conditional_values(parallel, option)
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
    selected = DiDInt.plot_layout(parallel; treated_states = ["s03"], ccc = "hom")
    @test selected.panels == [DiDInt.PLOT_CCC_LABELS["hom"]]
    @test sort(unique(selected.series.group)) == sort(vcat("s03", PLOT_NEVER_TREATED))
    @test all(panel_series(selected, 1; group = "s03").color .== 3) && selected.vlines.x == [5]
    @test DiDInt.plot_layout(parallel; ccc = "all").panels == layout.panels
    @test DiDInt.plot_layout(parallel; window = (2, 5)).xticks == (collect(2:5), string.(2003:2006))

    # Every treated state has its own colour and entry in the legend, however many there are,
    # and keeps its colour when only some treated states are selected
    many = quiet_call() do
        DiDInt.didint_plot("coll", "state", "year", TEST_DATA_FULL; treated_states = ["88", "57", "58", "59", "61", "64", "71", "72", "85", "34"],
                           treatment_times = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000], ccc = "none")
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
    unknown = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES, ccc = ["hom", "none"])
    end
    alike = DiDInt.plot_layout(unknown)
    @test all(alike.series.role .== :state) && length(unique(alike.series.group)) == 12 && alike.vlines.x == [3, 5, 6, 3, 5, 6]
    @test_throws r"were not given to didint_plot" DiDInt.plot_layout(unknown; treated_states = "s01")
    @test_throws r"needs the treated states" DiDInt.plot_layout(unknown; facet = "state", ccc = "hom")

    # A state treated throughout is neither a treated state nor a control state
    early = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treated_states = vcat(INTEGRATION_TREATED, "s06"),
                           treatment_times = vcat(INPUT_TIMES, 2001), ccc = "none")
    end
    @test !("s06" in DiDInt.plot_layout(early).series.group)
    @test !("s06" in DiDInt.plot_layout(early; facet = "state", controls = "separate").series.group)

    @test_throws r"other than \"int\"" DiDInt.plot_layout(plot_object(ccc = "int"); plot_type = "conditional")
    @test_throws r"'facet' must be" DiDInt.plot_layout(parallel; facet = "cohort")
    @test_throws r"'controls' must be" DiDInt.plot_layout(parallel; controls = "mean")
end

@testset "Conditional parallel trends plot, one panel per treated state" begin
    parallel = plot_object(ccc = ["hom", "int"])
    lambda, n = conditional_values(parallel, "hom")

    # The treated state next to the average control state. The average starts at the controls'
    # weighted mean and then adds the weighted mean of their changes, with weights n_t + n_{t-1}
    layout = DiDInt.plot_layout(parallel; plot_type = "conditional", facet = "state")
    @test layout.panels == PLOT_STATE_TITLES && occursin("ccc = \"hom\"", layout.title)
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
    equal = DiDInt.plot_layout(parallel; plot_type = "conditional", facet = "state", weights = false, treated_states = "s04")
    @test equal.panels == ["s04 (treated 2006)"]
    average = panel_series(equal, 1; role = :control_mean)
    @test all(average.label .== "Control states (mean)")
    @test average.y ≈ [mean(lambda[(c, period)] for c in PLOT_NEVER_TREATED) for period in 0:7] atol = 1e-12

    # One line per control state instead
    separate = DiDInt.plot_layout(parallel; plot_type = "conditional", facet = "state", controls = "separate")
    for panel in 1:5
        controls = panel_series(separate, panel; role = :control)
        @test sort(unique(controls.group)) == PLOT_NEVER_TREATED
        @test all(row.y == lambda[(row.group, row.x)] for row in eachrow(controls))
        @test nrow(panel_series(separate, panel; role = :treated)) == 8
    end

    # notyet = true: a state that is not yet treated is a control up to the period before its
    # treatment. For s01 (treated in period 3): s02 until period 2, s03 and s04 until 4, s05 until 5.
    with_notyet = plot_object(ccc = "hom", notyet = true)
    lambda, n = conditional_values(with_notyet, "hom")
    separate = DiDInt.plot_layout(with_notyet; facet = "state", controls = "separate", treated_states = "s01")
    @test panel_series(separate, 1; group = "s02", role = :control).x == 0:2
    @test panel_series(separate, 1; group = "s03", role = :control).x == 0:4
    @test panel_series(separate, 1; group = "s05", role = :control).x == 0:5
    @test panel_series(separate, 1; group = "s07", role = :control).x == 0:7
    controls_at(period) = vcat(PLOT_NEVER_TREATED, [j for (j, g) in zip(INTEGRATION_TREATED, PLOT_TREATMENT_PERIODS) if j != "s01" && g > period])
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

    # A single ccc option
    @test_throws r"single 'ccc' option" DiDInt.plot_layout(plot_object(); facet = "state")
    @test DiDInt.plot_layout(plot_object(); facet = "state", ccc = "time").title == "Conditional parallel trends by treated state, ccc = \"time\""
end

@testset "Integrated parallel trends plot" begin
    parallel = plot_object()
    rows = parallel.data[isequal.(parallel.data.ccc, "int"), :]

    # In its first treated period, a treated state's change minus the average of its controls'
    # changes is its ATT(s,g) under ccc = "int": with the controls weighted by their numbers
    # of observations under weighting = "both", and equally under weighting = "none"
    for (weights, weighting) in [(true, "both"), (false, "none")]
        layout = DiDInt.plot_layout(parallel; plot_type = "integrated", weights = weights)
        @test layout.plot_type == "integrated" && layout.panels == PLOT_STATE_TITLES
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
    separate = DiDInt.plot_layout(parallel; plot_type = "integrated", controls = "separate", treated_states = ["s05", "s02"])
    @test separate.panels == ["s02 (treated 2004)", "s05 (treated 2007)"]
    for (panel, s) in enumerate(["s02", "s05"])
        delta = Dict((row.state, row.period) => row.delta for row in eachrow(rows[rows.treated_state_covars .== s, :]))
        series = panel_series(separate, panel)
        @test sort(unique(series.group)) == sort(vcat(s, PLOT_NEVER_TREATED))
        @test all(row.y == delta[(row.group, row.x)] for row in eachrow(series))
        @test all((row.role == :treated) == (row.group == s) for row in eachrow(series))
    end

    # notyet = true adds the not-yet-treated states to the controls of the average
    with_notyet = plot_object(ccc = "int", notyet = true)
    @test DiDInt.available_plot_types(with_notyet) == ["integrated", "integrated2"]
    layout = DiDInt.plot_layout(with_notyet; weights = false, treated_states = "s01")
    @test layout.plot_type == "integrated"
    notyet_rows = with_notyet.data[with_notyet.data.treated_state_covars .== "s01", :]
    others = notyet_rows[notyet_rows.state .!= "s01", :]
    @test panel_series(layout, 1; role = :control_mean).y ≈ [mean(others.delta[others.period .== period]) for period in 1:7] atol = 1e-12
    @test count(others.period .== 1) == length(PLOT_NEVER_TREATED) + 4 && count(others.period .== 7) == length(PLOT_NEVER_TREATED)

    # Titles, and what cannot be drawn
    @test DiDInt.plot_layout(parallel; plot_type = "integrated", title = "My title").title == "My title"
    @test DiDInt.plot_layout(parallel; plot_type = "integrated", title = "").title == ""
    @test DiDInt.plot_layout(parallel; plot_type = "integrated", ccc = "int").panels == PLOT_STATE_TITLES
    @test_throws r"cannot select another option" DiDInt.plot_layout(parallel; plot_type = "integrated", ccc = "hom")
    @test_throws r"no integrated parallel trends" DiDInt.plot_layout(plot_object(ccc = "hom"); plot_type = "integrated")
end

@testset "Conditional parallel trends relative to the period before treatment (conditional2)" begin
    # The average control state is measured from the period before the treated state's
    # treatment with the differences and weights of didint(). The gap between the treated
    # state and that average, less the gap in the period before treatment, is therefore
    # ATT(s,t) in every treated period: with the controls weighted by their numbers of
    # observations under weighting = "both" and equally under "none", with and without the
    # not-yet-treated states as controls.
    for notyet in [false, true], (weights, weighting) in [(true, "both"), (false, "none")], option in ["hom", "add"]
        object = plot_object(ccc = option, notyet = notyet)
        layout = DiDInt.plot_layout(object; plot_type = "conditional2", weights = weights)
        @test layout.plot_type == "conditional2" && layout.panels == PLOT_STATE_TITLES
        @test layout.vlines.x == PLOT_TREATMENT_PERIODS && layout.xticks == (collect(0:7), string.(2001:2008))
        sgt = integration_didint(ccc = option, agg = "sgt", weighting = weighting, notyet = notyet)
        for (panel, s) in enumerate(INTEGRATION_TREATED)
            base = PLOT_TREATMENT_PERIODS[panel] - 1
            own = panel_series(layout, panel; role = :treated)
            average = panel_series(layout, panel; role = :control_mean)
            @test own.x == 0:7 && average.x == 0:7
            gap(period) = only(own.y[own.x .== period]) - only(average.y[average.x .== period])
            estimates = sgt[sgt.state .== s, :]
            @test nrow(estimates) == 8 - PLOT_TREATMENT_PERIODS[panel]
            @test all(isapprox(gap(year(row.t) - 2001) - gap(base), row.att_sgt; atol = 1e-8) for row in eachrow(estimates))
        end
    end

    # In the period before treatment the line is at the controls' average lambda; with equal
    # weights and the same controls in every period it is their average in every period, as
    # is the line of plot_type = "conditional" with facet = "state"
    parallel = plot_object(ccc = "hom")
    lambda, n = conditional_values(parallel, "hom")
    weighted = DiDInt.plot_layout(parallel; plot_type = "conditional2", treated_states = "s03")
    average = panel_series(weighted, 1; role = :control_mean)
    @test only(average.y[average.x .== 4]) ≈ sum(n[(c, 4)] * lambda[(c, 4)] for c in PLOT_NEVER_TREATED) / sum(n[(c, 4)] for c in PLOT_NEVER_TREATED) atol = 1e-12
    equal = DiDInt.plot_layout(parallel; plot_type = "conditional2", weights = false)
    chained = DiDInt.plot_layout(parallel; facet = "state", weights = false)
    @test equal.series.y ≈ chained.series.y atol = 1e-12
    @test panel_series(equal, 1; role = :control_mean).y ≈ [mean(lambda[(c, period)] for c in PLOT_NEVER_TREATED) for period in 0:7] atol = 1e-12

    # One line per control state: the lambdas themselves
    @test isequal(DiDInt.plot_layout(parallel; plot_type = "conditional2", controls = "separate").series,
                  DiDInt.plot_layout(parallel; facet = "state", controls = "separate").series)
    # ... where a not-yet-treated state of another cohort is a control while it is untreated in
    # both periods. For s03 (treated in period 5, compared with period 4): s05 (treated in
    # period 6) until period 5; s01, s02 (treated before) and s04 (same cohort) never.
    separate = DiDInt.plot_layout(plot_object(ccc = "hom", notyet = true); plot_type = "conditional2", controls = "separate", treated_states = "s03")
    @test panel_series(separate, 1; group = "s05", role = :control).x == 0:5
    @test sort(unique(separate.series.group)) == sort(vcat("s03", "s05", PLOT_NEVER_TREATED))

    # Controls that are missing in a period or in the base period are left out; without any
    # control in the base period there is no line
    gaps = Dict(("a", 0) => 1.0, ("a", 1) => 2.0, ("b", 0) => 3.0, ("b", 1) => 3.0, ("a", 3) => 7.0, ("b", 3) => 9.0, ("a", 4) => 8.0)
    counts = Dict(key => 1 for key in keys(gaps))
    periods, values = DiDInt.anchored_average_of_controls(gaps, counts, period -> ["a", "b"], 0:4, 1, true)
    @test periods == [0, 1, 3, 4] && values ≈ [2.0, 2.5, 8.0, 8.5]
    @test DiDInt.anchored_average_of_controls(gaps, counts, period -> ["a", "b"], 0:4, 2, true) == (Int[], Float64[])

    # It needs the treated states and a single ccc option
    unknown = quiet_call() do
        DiDInt.didint_plot("y", "state", "year", INTEGRATION_PANEL; treatment_times = INPUT_TIMES, ccc = "hom")
    end
    @test DiDInt.available_plot_types(unknown) == ["conditional"]
    @test_throws r"\"conditional2\" needs the treated states" DiDInt.plot_layout(unknown; plot_type = "conditional2")
    @test_throws r"\"conditional2\" shows a single 'ccc' option" DiDInt.plot_layout(plot_object(); plot_type = "conditional2")
    @test_throws r"other than \"int\"" DiDInt.plot_layout(plot_object(ccc = "int"); plot_type = "conditional2")
end

@testset "Integrated parallel trends since the period before treatment (integrated2)" begin
    predicted(state, yr, x) = (fit = cell_ols(INTEGRATION_PANEL, state, yr).beta; fit[1] + dot(fit[2:end], x))
    cell_n(state, yr) = nrow(cell_rows(INTEGRATION_PANEL, state, yr))
    treatment_period = Dict(zip(INTEGRATION_TREATED, PLOT_TREATMENT_PERIODS))

    for notyet in [false, true]
        object = plot_object(ccc = "int", notyet = notyet)
        long = object.long_differences
        @test names(long) == ["state", "time", "delta", "n", "treated_state_covars", "period", "base_period", "treat_period"]
        @test all(row.time == string(2001 + row.period) && row.period != row.base_period for row in eachrow(long))
        @test all(row.treat_period == treatment_period[row.treated_state_covars] && row.base_period == row.treat_period - 1 for row in eachrow(long))

        # Every value is the change in the cell regressions' predictions between the period and
        # the period before the treated state's treatment, at the covariate means of the
        # treated state in the period; n counts the state's observations in the two periods
        for row in eachrow(long)
            yr, base_yr = 2001 + row.period, 2001 + row.base_period
            x = cell_xbar(INTEGRATION_PANEL, row.treated_state_covars, yr)
            @test row.delta ≈ predicted(row.state, yr, x) - predicted(row.state, base_yr, x) atol = 1e-8
            @test row.n == cell_n(row.state, yr) + cell_n(row.state, base_yr)
        end

        # The comparison states: the never-treated states and, with notyet, the states of other
        # cohorts that are untreated in both periods
        if !notyet
            @test nrow(long) == 5 * (1 + length(PLOT_NEVER_TREATED)) * 7
        else
            s01 = long[long.treated_state_covars .== "s01", :]
            @test s01.period[s01.state .== "s03"] == [0, 1, 3, 4] && s01.period[s01.state .== "s05"] == [0, 1, 3, 4, 5]
            @test !("s02" in s01.state[s01.state .!= "s01"])
        end

        # The treated state's change minus the average of its controls' changes is ATT(s,t) in
        # every treated period, and both are zero in the period before treatment
        for (weights, weighting) in [(true, "both"), (false, "none")]
            layout = DiDInt.plot_layout(object; plot_type = "integrated2", weights = weights)
            @test layout.plot_type == "integrated2" && layout.panels == PLOT_STATE_TITLES
            @test layout.vlines.x == PLOT_TREATMENT_PERIODS && layout.xticks == (collect(0:7), string.(2001:2008))
            sgt = integration_didint(ccc = "int", agg = "sgt", weighting = weighting, notyet = notyet)
            for (panel, s) in enumerate(INTEGRATION_TREATED)
                own = panel_series(layout, panel; role = :treated)
                average = panel_series(layout, panel; role = :control_mean)
                @test own.x == 0:7 && average.x == 0:7
                base = treatment_period[s] - 1
                @test only(own.y[own.x .== base]) == 0 && only(average.y[average.x .== base]) == 0
                gap(period) = only(own.y[own.x .== period]) - only(average.y[average.x .== period])
                estimates = sgt[sgt.state .== s, :]
                @test nrow(estimates) == 8 - treatment_period[s]
                @test all(isapprox(gap(year(row.t) - 2001), row.att_sgt; atol = 1e-8) for row in eachrow(estimates))
            end
        end
    end

    # In the first treated period the change since the period before treatment is the change
    # since the previous period: the two integrated plots agree there
    parallel = plot_object()
    since_previous = DiDInt.plot_layout(parallel; plot_type = "integrated")
    since_base = DiDInt.plot_layout(parallel; plot_type = "integrated2")
    for (panel, g) in enumerate(PLOT_TREATMENT_PERIODS), role in [:treated, :control_mean]
        @test only(panel_series(since_previous, panel; role = role).y[panel_series(since_previous, panel; role = role).x .== g]) ≈
              only(panel_series(since_base, panel; role = role).y[panel_series(since_base, panel; role = role).x .== g]) atol = 1e-10
    end

    # One line per control state, each through zero in the period before treatment
    separate = DiDInt.plot_layout(parallel; plot_type = "integrated2", controls = "separate", treated_states = "s05")
    @test separate.panels == ["s05 (treated 2007)"]
    @test sort(unique(separate.series.group)) == sort(vcat("s05", PLOT_NEVER_TREATED))
    @test all(separate.series.y[separate.series.x .== 5] .== 0) && count(separate.series.x .== 5) == 1 + length(PLOT_NEVER_TREATED)

    @test_throws r"cannot select another option" DiDInt.plot_layout(parallel; plot_type = "integrated2", ccc = "hom")
    @test_throws r"no integrated parallel trends" DiDInt.plot_layout(plot_object(ccc = "hom"); plot_type = "integrated2")
end

@testset "plot() draws every kind of plot with Makie" begin
    @test !isnothing(Base.get_extension(DiDInt, :DiDIntMakieExt))
    parallel = plot_object(notyet = true)
    event = plot_object(event = true)
    axes_of(figure) = [block for block in figure.content if block isa CairoMakie.Axis]
    legends_of(figure) = [block for block in figure.content if block isa CairoMakie.Legend]

    cases = [(event, (;), 6),
             (event, (ccc = ["none", "int"], treated_states = ["s01", "s03", "s05"], window = (-3, 2), groupmin = 2), 2),
             (parallel, (;), 5),
             (parallel, (plot_type = "conditional", treated_states = "s03", ccc = "hom"), 1),
             (parallel, (facet = "state", ccc = "add"), 5),
             (parallel, (facet = "state", ccc = "state", controls = "separate", weights = false), 5),
             (parallel, (plot_type = "integrated",), 5),
             (parallel, (plot_type = "integrated", controls = "separate", treated_states = ["s01", "s05"]), 2),
             (parallel, (plot_type = "conditional2", ccc = "hom"), 5),
             (parallel, (plot_type = "conditional2", ccc = "time", controls = "separate", treated_states = "s04"), 1),
             (parallel, (plot_type = "integrated2",), 5),
             (parallel, (plot_type = "integrated2", controls = "separate", weights = false, treated_states = ["s03", "s02"]), 2)]
    for (object, arguments, n_panels) in cases
        figure = CairoMakie.plot(object; arguments...)
        @test figure isa CairoMakie.Figure
        @test length(axes_of(figure)) == n_panels && length(legends_of(figure)) == 1
        file = tempname() * ".png"
        CairoMakie.save(file, figure)
        @test filesize(file) > 5_000
        rm(file; force = true)
    end

    # Appearance: number of columns, size, title and theme
    figure = CairoMakie.plot(parallel; plot_type = "integrated", ncols = 1, size = (500, 1100), title = "", theme = CairoMakie.Theme())
    @test figure isa CairoMakie.Figure && length(axes_of(figure)) == 5
    @test Tuple(figure.scene.viewport[].widths) == (500, 1100)
    @test_throws r"'ncols' must be" CairoMakie.plot(parallel; ncols = 0)
    @test_throws r"event = true" CairoMakie.plot(parallel; plot_type = "event")

    # Every coloured line can be told from the others, also when there are more of them than colours
    extension = Base.get_extension(DiDInt, :DiDIntMakieExt)
    @test allunique((extension.series_color(index), extension.series_linestyle(index)) for index in 1:48)
    @test length(unique(extension.series_color.(1:16))) == 16 && all(extension.series_linestyle.(1:16) .== :solid)
end
