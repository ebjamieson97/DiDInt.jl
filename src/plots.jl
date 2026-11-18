###########################################################################################
# plots.jl: the object that didint_plot() returns and the layouts of its plots            #
#                                                                                         #
#   Section 1: DiDIntPlot                                                                 #
#   Section 2: Display and conversion to a table                                          #
#   Section 3: Plot layouts: the arguments of plot() and what each plot shows             #
#   Section 4: Event study layout                                                         #
#   Section 5: Conditional parallel trends layouts                                        #
#   Section 6: Integrated parallel trends layout                                          #
#                                                                                         #
# The layouts are plain tables of lines and reference marks. They are drawn by the Makie  #
# extension (ext/DiDIntMakieExt.jl), which is loaded together with a Makie backend.       #
###########################################################################################


###########################################################################################
# Section 1: DiDIntPlot                                                                   #
###########################################################################################

"""
    DiDIntPlot

The data that [`didint_plot`](@ref) prepares, which can be drawn with `plot`. Displaying it prints a summary.

# Fields
- `data::DataFrame`: the plot data, as described under "Returns" in [`didint_plot`](@ref): the parallel trends
  data (`event = false`) or the event study data (`event = true`).
- `lambdas::DataFrame`: for an event study, the `λ(s,t)` of each treated state behind `data`, by `ccc` option and
  period relative to treatment (`time_since_treatment`), with the state's number of observations (`n`). Empty
  for parallel trends data.
- `specs::NamedTuple`: the specification (`outcome`, `event`, `ccc`, `covariates`, `treated_states` with their
  `treatment_times` and `treatment_periods`, `notyet`, `weights`, `ci`, the period length and date range, and so on).

`DataFrame(x)` returns the plot data (`x.data`). A `DiDIntPlot` is also a Tables.jl table, so it can be passed
directly to functions such as `CSV.write`.

# Plotting
Plots are drawn with [Makie](https://docs.makie.org): load one of its backends and call `plot`, which returns a
Makie `Figure`.

```julia
using CairoMakie
x = didint_plot("coll", "state", "year", data;
                treated_states = states, treatment_times = times,
                covariates = ["asian", "black", "male"])

# Conditional parallel trends, one panel per ccc option
plot(x)

# Integrated parallel trends, one panel per treated state
plot(x; plot_type = "integrated")

# Two treated states next to their control states, under ccc = "hom"
figure = plot(x; facet = "state", ccc = "hom", treated_states = ["71", "58"])
save("trends.png", figure)
```

    plot(x::DiDIntPlot; plot_type, treated_states, ccc, window,
         facet, controls, weights, groupmin, ncols, size, title, theme)

## What is plotted
- `plot_type::Union{AbstractString, Nothing} = nothing`
    `"event"`, `"conditional"` or `"integrated"`. An event study plot needs data made with `event = true`, and
    the other two data made with `event = false`. By default: the event study for event study data; otherwise
    the conditional parallel trends if the data have them and the integrated parallel trends if not.
- `treated_states = nothing`
    The treated state(s) to show, as a single value or a vector. By default all of them.
- `ccc = nothing`
    The `ccc` option(s) to show in an event study or conditional parallel trends plot, as a string or a vector
    of strings. By default all of those in the data.
- `window = nothing`
    A pair `(first, last)` that limits the horizontal axis: periods relative to treatment for an event study,
    period numbers (the `period` column) otherwise.

## Parallel trends plots
- `facet::AbstractString = "ccc"`
    For conditional parallel trends: `"ccc"` draws one panel per `ccc` option with a line for every state;
    `"state"` draws one panel per treated state, comparing it with its control states under a single `ccc`
    option (chosen with `ccc` when the data have several).
- `controls::AbstractString = "average"`
    In the panels of a treated state (integrated parallel trends, and conditional parallel trends with
    `facet = "state"`): `"average"` draws one line for the average of the control states, `"separate"` a line
    for each control state.
- `weights::Union{Bool, Nothing} = nothing`
    Whether the average of the control states weights them by their numbers of observations (`true`, the
    default) or equally (`false`). For an event study: whether the treated states are weighted by their numbers
    of observations; the default is the choice made in `didint_plot()`.

## Event study plots
- `groupmin::Integer = 3`
    The smallest number of treated states behind an estimate for its confidence interval to be drawn. Intervals
    that rest on very few states are too wide to be of use.

## Appearance
- `ncols::Union{Integer, Nothing} = nothing`: number of columns of panels.
- `size::Union{Tuple, Nothing} = nothing`: size of the figure, as `(width, height)`.
- `title::Union{AbstractString, Nothing} = nothing`: title of the figure; `""` for none.
- `theme`: a Makie `Theme` that replaces the default look.

The control states of a treated state are the never-treated states and, if `didint_plot()` was called with
`notyet = true`, the states not yet treated in the period in question. A dotted vertical line marks the first
treated period.
"""
struct DiDIntPlot
    data::DataFrame
    lambdas::DataFrame
    specs::NamedTuple
end


###########################################################################################
# Section 2: Display and conversion to a table                                            #
###########################################################################################

# The ccc options of the conditional parallel trends data, in the order of their panels, and
# the titles of the panels
const PLOT_CCC_ORDER = ["none", "hom", "time", "state", "add", "int"]
const PLOT_CCC_LABELS = Dict("none" => "No covariates (ccc = \"none\")", "hom" => "Homogeneous (ccc = \"hom\")",
                             "time" => "Time-varying (ccc = \"time\")", "state" => "Region-varying (ccc = \"state\")",
                             "add" => "Two one-way (ccc = \"add\")", "int" => "Two-way (ccc = \"int\")")

# Number of rows of the data printed in the summary
const PLOT_ROWS_SHOWN = 5

# The plot types that the data of a DiDIntPlot allow
function available_plot_types(plot::DiDIntPlot)
    plot.specs.event && return ["event"]
    types = String[]
    any(!ismissing, plot.data.lambda) && push!(types, "conditional")
    any(!ismissing, plot.data.delta) && push!(types, "integrated")
    return types
end

function Base.show(io::IO, plot::DiDIntPlot)

    # One-line form, used inside collections
    print(io, "DiDIntPlot(", join(available_plot_types(plot), " and "), " plot data, ", nrow(plot.data), " rows)")
end

function wrapped_list(items, indent; width = 72)

    # The items separated by commas, on lines of at most `width` characters, with the lines
    # after the first one indented
    lines = String[]
    line = ""
    for item in items
        if !isempty(line) && length(line) + 2 + length(item) > width
            push!(lines, line * ",")
            line = item
        else
            line = isempty(line) ? item : line * ", " * item
        end
    end
    push!(lines, line)
    return join(lines, "\n" * " " ^ indent)
end

function Base.show(io::IO, ::MIME"text/plain", plot::DiDIntPlot)

    # Summary: the specification, the plots that can be drawn and the first rows of the data
    specs = plot.specs
    types = available_plot_types(plot)
    println(io, "DID-INT plot data: ", specs.event ? "event study" : "parallel trends")
    println(io, "  Outcome:          ", specs.outcome)
    if !isempty(specs.covariates)
        println(io, "  Covariates:       ", join(specs.covariates, ", "))
    end
    println(io, "  ccc options:      ", join(sort(specs.ccc, by = c -> findfirst(==(c), PLOT_CCC_ORDER)), ", "))
    if isempty(specs.treated_states)
        println(io, "  Treated states:   not specified (treatment times: ", join(specs.treatment_times, ", "), ")")
    else
        treated = [string(s, " (", t, ")") for (s, t) in zip(specs.treated_states, specs.treatment_times)]
        println(io, "  Treated states:   ", wrapped_list(treated, 20))
    end
    if "integrated" in types
        println(io, "  Controls:         ", specs.notyet ? "never-treated and not-yet-treated states" : "never-treated states")
    end
    println(io, "  Periods:          ", specs.period, ", from ", specs.start_date, " to ", specs.end_date)
    println(io, "  Plot types:       ", join(("\"" * type * "\"" for type in types), ", "), " (drawn with plot; see ?DiDIntPlot)")
    println(io)
    println(io, "Data (first ", min(PLOT_ROWS_SHOWN, nrow(plot.data)), " of ", nrow(plot.data), " rows; see the `data` field)")
    show(io, first(plot.data, PLOT_ROWS_SHOWN); summary = false, eltypes = false, show_row_number = false)
end

# A DiDIntPlot is a Tables.jl table: its plot data. This is what makes DataFrame(plot) and,
# for example, CSV.write(file, plot) work.
Tables.istable(::Type{DiDIntPlot}) = true
Tables.columnaccess(::Type{DiDIntPlot}) = true
Tables.columns(plot::DiDIntPlot) = Tables.columns(plot.data)
Tables.schema(plot::DiDIntPlot) = Tables.schema(plot.data)


###########################################################################################
# Section 3: Plot layouts: the arguments of plot() and what each plot shows               #
###########################################################################################

# What a plot shows, independently of the plotting package that draws it.
#   panels:  the titles of the panels, in order
#   series:  one row per point, with the panel (its number), the line the point belongs to
#            (`group`), the entry of the legend that the line belongs to (`label`), what the
#            line is (`role`: :estimate, :treated, :state, :control or :control_mean), the
#            number of its colour among the coloured lines (`color`, 0 for control lines), its
#            coordinates (`x`, `y`) and, for an event study, the confidence bounds
#   vlines:  the vertical reference lines (first treated periods), by panel
#   xticks:  for a period axis, the period numbers and the labels of their dates
#   legend_extras: entries of the legend for marks other than the lines, as kind => label,
#            where the kind is :band (confidence interval), :no_interval (estimate without
#            one) or :vline (vertical reference line)
struct PlotLayout
    plot_type::String
    title::String
    xlabel::String
    ylabel::String
    legend_title::String
    panels::Vector{String}
    series::DataFrame
    vlines::DataFrame
    xticks::Union{Nothing, Tuple{Vector{Int}, Vector{String}}}
    legend_extras::Vector{Pair{Symbol, String}}
end

# Largest number of treated states that are given a colour each; with more, they share one
const MAX_STATE_COLORS = 8

empty_series() = DataFrame(panel = Int[], group = String[], label = String[], role = Symbol[], color = Int[], x = Int[],
                           y = Float64[], lower = Union{Missing, Float64}[], upper = Union{Missing, Float64}[])

function add_line!(series, panel, group, label, role, color, x, y; lower = missing, upper = missing)

    # Add the points of one line to the table of series
    for i in eachindex(x)
        push!(series, (panel = panel, group = group, label = label, role = role, color = color, x = x[i], y = y[i],
                       lower = lower isa AbstractVector ? lower[i] : lower, upper = upper isa AbstractVector ? upper[i] : upper))
    end
    return series
end

function plot_layout(plot::DiDIntPlot; plot_type::Union{AbstractString, Nothing} = nothing, treated_states = nothing,
                     ccc::Union{AbstractString, AbstractVector{<:AbstractString}, Nothing} = nothing, facet::AbstractString = "ccc",
                     controls::AbstractString = "average", weights::Union{Bool, Nothing} = nothing, window = nothing,
                     groupmin::Integer = 3, title::Union{AbstractString, Nothing} = nothing)

    # The layout of the plot that the arguments of plot() ask for
    plot_type = check_plot_type(plot, plot_type)
    treated = selected_treated_states(plot, treated_states)
    facet = lowercase(strip(facet))
    facet in ("ccc", "state") || error("'facet' must be \"ccc\" or \"state\".")
    controls = lowercase(strip(controls))
    controls in ("average", "separate") || error("'controls' must be \"average\" or \"separate\".")
    window = check_window(window)

    if plot_type == "event"
        layout = event_layout(plot, treated, ccc, weights, groupmin)
    elseif plot_type == "integrated"
        layout = integrated_layout(plot, treated, ccc, controls, something(weights, true))
    elseif facet == "ccc"
        layout = conditional_layout(plot, treated, ccc)
    else
        layout = conditional_state_layout(plot, treated, ccc, controls, something(weights, true))
    end
    layout = apply_window(layout, window)
    return isnothing(title) ? layout : with_title(layout, title)
end

function check_plot_type(plot::DiDIntPlot, plot_type)

    # The plot type asked for, or the default one for the data: the event study for event
    # study data; otherwise the conditional parallel trends if the data have them and the
    # integrated parallel trends if not
    available = available_plot_types(plot)
    isempty(available) && error("The plot data are empty.")
    isnothing(plot_type) && return first(available)
    plot_type = lowercase(strip(plot_type))
    plot_type in ("event", "conditional", "integrated") || error("'plot_type' must be \"event\", \"conditional\" or \"integrated\".")
    plot_type in available && return plot_type
    if plot_type == "event"
        error("These are parallel trends data. An event study plot needs data from didint_plot() with event = true.")
    elseif plot.specs.event
        error("These are event study data. A parallel trends plot needs data from didint_plot() with event = false.")
    elseif plot_type == "integrated"
        error("The data have no integrated parallel trends: call didint_plot() with \"int\" among the 'ccc' options " *
              "(and with the treated states).")
    end
    error("The data have no conditional parallel trends: call didint_plot() with a 'ccc' option other than \"int\".")
end

function selected_treated_states(plot::DiDIntPlot, treated_states)

    # The treated states to show: all of them, or those given in `treated_states` (a single
    # value or a vector, of strings or numbers)
    known = plot.specs.treated_states
    isnothing(treated_states) && return known
    isempty(known) && error("'treated_states' cannot be used here: the treated states were not given to didint_plot().")
    selected = unique(string.(treated_states isa AbstractVector ? treated_states : [treated_states]))
    isempty(selected) && error("'treated_states' is empty.")
    unknown = setdiff(selected, known)
    isempty(unknown) || error("Not among the treated states of the plot data: $(join(unknown, ", ")). The treated states are: " *
                              join(known, ", ") * ".")
    # in the order of the plot data
    return filter(in(selected), known)
end

function selected_ccc(available, ccc)

    # The ccc options to show, in the order of the panels: those in the data, or those given in `ccc`
    available = sort(available, by = c -> findfirst(==(c), PLOT_CCC_ORDER))
    isnothing(ccc) && return available
    selected = unique(replace.(lowercase.(ccc isa AbstractVector ? ccc : [ccc]), r"\s" => ""))
    selected == ["all"] && return available
    unknown = setdiff(selected, available)
    isempty(unknown) || error("No data to plot for ccc = $(join(("\"" * c * "\"" for c in unknown), ", ")). The options that can be " *
                              "plotted here are: " * join(available, ", ") * ".")
    return filter(in(selected), available)
end

function check_window(window)

    # `window` as a pair of numbers (first, last), or nothing
    isnothing(window) && return nothing
    if !(window isa Union{Tuple, AbstractVector}) || length(window) != 2 || !all(w -> w isa Real, window) || window[1] > window[2]
        error("'window' must be a pair of numbers (first, last), with first <= last.")
    end
    return (window[1], window[2])
end

function apply_window(layout::PlotLayout, window)

    # Keep what lies within the window of the horizontal axis
    isnothing(window) && return layout
    inside(x) = window[1] <= x <= window[2]
    series = layout.series[inside.(layout.series.x), :]
    nrow(series) == 0 && error("No data within 'window' = $window.")
    xticks = layout.xticks
    if !isnothing(xticks)
        keep = inside.(xticks[1])
        xticks = (xticks[1][keep], xticks[2][keep])
    end
    # the legend only explains marks that are still there
    vlines = layout.vlines[inside.(layout.vlines.x), :]
    legend_extras = filter(layout.legend_extras) do (kind, _)
        kind == :no_interval ? any(ismissing, series.lower) : kind == :vline ? nrow(vlines) > 0 : true
    end
    return PlotLayout(layout.plot_type, layout.title, layout.xlabel, layout.ylabel, layout.legend_title, layout.panels, series,
                      vlines, xticks, legend_extras)
end

with_title(layout::PlotLayout, title) = PlotLayout(layout.plot_type, String(title), layout.xlabel, layout.ylabel, layout.legend_title,
                                                   layout.panels, layout.series, layout.vlines, layout.xticks, layout.legend_extras)

function outcome_label(plot::DiDIntPlot, options)

    # Label of the vertical axis: the outcome, which is adjusted for covariates unless there
    # are none or only ccc = "none" is shown
    adjusted = !isempty(plot.specs.covariates) && any(!=("none"), options)
    return adjusted ? plot.specs.outcome * " (adjusted for covariates)" : plot.specs.outcome
end

sorted_states(states) = sort(String.(states), by = custom_sort_order)

function period_ticks(rows)

    # The period numbers of the rows with the labels of their dates
    periods = sort(unique(rows.period))
    label = Dict(zip(rows.period, rows.time))
    return (Int.(periods), String[label[p] for p in periods])
end

function weighted_mean(values, weights)
    return sum(values .* weights) / sum(weights)
end


###########################################################################################
# Section 4: Event study layout                                                           #
###########################################################################################

function event_layout(plot::DiDIntPlot, treated, ccc, weights, groupmin)

    # One panel per ccc option: the mean of the treated states' lambdas by period relative to
    # treatment, with its confidence interval where at least `groupmin` treated states are
    # behind it (an interval from very few states is too wide to be of use). When only some of
    # the treated states are selected (or another weighting is asked for), the means are
    # computed again from the lambdas of those states.
    specs = plot.specs
    weights = something(weights, specs.weights)
    if treated == specs.treated_states && weights == specs.weights
        table = plot.data
    else
        table = event_study_data(plot.lambdas[in.(plot.lambdas.state, Ref(treated)), :], weights, specs.ci, specs.hc)
    end
    options = selected_ccc(unique(table.ccc), ccc)

    series = empty_series()
    for (panel, option) in enumerate(options)
        rows = table[table.ccc .== option, :]
        enough = rows.ngroup .>= groupmin
        add_line!(series, panel, option, "Mean of the treated states", :estimate, 1, rows.time_since_treatment, rows.y;
                  lower = ifelse.(enough, rows.ci_lower, missing), upper = ifelse.(enough, rows.ci_upper, missing))
    end
    vlines = DataFrame(panel = collect(eachindex(options)), x = zeros(Int, length(options)))

    # The legend explains the band and, if some estimates have no interval, their marker
    level = 100 * specs.ci
    legend_extras = [:band => string(isinteger(level) ? Int(level) : level, "% confidence interval")]
    if any(ismissing, series.lower)
        push!(legend_extras, :no_interval => "No interval (fewer than $(max(groupmin, 2)) treated states)")
    end
    return PlotLayout("event", "Event study", "Periods relative to treatment", outcome_label(plot, options), "",
                      [PLOT_CCC_LABELS[option] for option in options], series, vlines, nothing, legend_extras)
end


###########################################################################################
# Section 5: Conditional parallel trends layouts                                          #
###########################################################################################

function conditional_rows(plot::DiDIntPlot)

    # The conditional parallel trends rows of the data: the lambda of every state and period
    # under every ccc option other than "int"
    data = plot.data
    return data[.!ismissing.(data.lambda), :]
end

function conditional_layout(plot::DiDIntPlot, treated, ccc)

    # One panel per ccc option, with the lambdas of every state as a line (Karim & Webb 2025,
    # eq. 45). The selected treated states are coloured and the states that are never treated
    # are drawn as control states; treated states that are not selected, and states treated
    # throughout, are left out. When the treated states are not known, every state is drawn
    # alike. Vertical lines mark the first treated periods.
    specs = plot.specs
    rows = conditional_rows(plot)
    options = selected_ccc(unique(rows.ccc), ccc)
    known = !isempty(specs.treated_states)

    # A colour for each treated state, kept when only some of them are selected; when there
    # are too many, the selected ones share a colour
    if length(specs.treated_states) <= MAX_STATE_COLORS
        color = Dict(s => i for (i, s) in enumerate(specs.treated_states))
    elseif length(treated) <= MAX_STATE_COLORS
        color = Dict(s => i for (i, s) in enumerate(treated))
    else
        color = Dict(s => 1 for s in treated)
    end
    one_color = length(treated) > MAX_STATE_COLORS

    series = empty_series()
    for (panel, option) in enumerate(options)
        option_rows = rows[rows.ccc .== option, :]
        for state_rows in groupby(option_rows, :state; sort = false)
            state = first(state_rows.state)
            state_rows = sort(state_rows, :period)
            if !known
                add_line!(series, panel, state, "States", :state, 1, state_rows.period, state_rows.lambda)
            elseif state in treated
                add_line!(series, panel, state, one_color ? "Treated states" : state, :treated, color[state], state_rows.period,
                          state_rows.lambda)
            elseif !(state in specs.treated_states) && !(state in specs.always_treated)
                add_line!(series, panel, state, "Control states", :control, 0, state_rows.period, state_rows.lambda)
            end
        end
    end

    # First treated periods of the selected treated states (or all treatment periods)
    treatment_periods = known ? [p for (s, p) in zip(specs.treated_states, specs.treatment_periods) if s in treated] : specs.treatment_periods
    treatment_periods = sort(unique(treatment_periods))
    vlines = DataFrame(panel = repeat(eachindex(options), inner = length(treatment_periods)),
                       x = repeat(treatment_periods, outer = length(options)))
    return PlotLayout("conditional", "Conditional parallel trends", "Time", outcome_label(plot, options),
                      known && !one_color ? "Treated state" : "", [PLOT_CCC_LABELS[option] for option in options], series, vlines,
                      period_ticks(rows), [:vline => (length(treatment_periods) == 1 ? "First treated period" : "First treated periods")])
end

function control_states_of(plot::DiDIntPlot, states)

    # The control states of a treated state s, as a function of the period: among `states`, the
    # never-treated states and, with `notyet`, the treated states other than s that are not
    # yet treated in that period. States treated throughout are never controls.
    specs = plot.specs
    treatment_period = Dict(zip(specs.treated_states, specs.treatment_periods))
    never_treated = sorted_states(setdiff(states, specs.treated_states, specs.always_treated))
    not_yet_candidates = sorted_states(intersect(specs.treated_states, states))
    return function (s, period)
        specs.notyet || return never_treated
        return vcat(never_treated, [j for j in not_yet_candidates if j != s && treatment_period[j] > period])
    end
end

function average_of_controls(lambda, n, controls_at, periods, weights)

    # The line of the average control state, built from changes, as didint() compares the
    # change of a treated state with the average of its controls' changes. It starts at the
    # controls' average lambda in the first period, and each later period adds the average of
    # the controls' changes since the previous period: over the controls of that period
    # (`controls_at(period)`) observed in both periods, weighted by their numbers of
    # observations in the two periods (as in didint()) or equally. With the same controls and
    # weights in every period this is the average of the controls' lambdas. Unlike an average
    # of the lambdas, it does not move when the set of controls or their weights change, nor
    # with shifts of the lambdas that are constant within a state.
    # Where no control has a change (a gap in the data), the line starts again from the
    # controls' average lambda.
    # Returns the periods that have a value, and the values.
    out_periods, values = Int[], Float64[]
    level = nothing
    for period in periods
        states = controls_at(period)
        changing = [c for c in states if haskey(lambda, (c, period)) && haskey(lambda, (c, period - 1))]
        if !isnothing(level) && !isempty(changing)
            w = weights ? Float64[n[(c, period)] + n[(c, period - 1)] for c in changing] : ones(length(changing))
            level += weighted_mean([lambda[(c, period)] - lambda[(c, period - 1)] for c in changing], w)
        else
            observed = [c for c in states if haskey(lambda, (c, period))]
            if isempty(observed)
                level = nothing
                continue
            end
            w = weights ? Float64[n[(c, period)] for c in observed] : ones(length(observed))
            level = weighted_mean([lambda[(c, period)] for c in observed], w)
        end
        push!(out_periods, period)
        push!(values, level)
    end
    return out_periods, values
end

function control_mean_label(weights)
    return weights ? "Control states (mean, weighted by observations)" : "Control states (mean)"
end

function conditional_state_layout(plot::DiDIntPlot, treated, ccc, controls, weights)

    # One panel per treated state, under a single ccc option: the lambdas of the treated state
    # and of its control states, the latter as their average (see average_of_controls) or as
    # one line per control state. A state that is not yet treated is a control up to the
    # period before its treatment.
    specs = plot.specs
    isempty(specs.treated_states) && error("facet = \"state\" needs the treated states: call didint_plot() with 'treated_states' or 'gvar'.")
    rows = conditional_rows(plot)
    options = selected_ccc(unique(rows.ccc), ccc)
    if length(options) != 1
        error("facet = \"state\" shows a single 'ccc' option: choose one of ", join(options, ", "), " with the 'ccc' argument.")
    end
    option = only(options)
    rows = rows[rows.ccc .== option, :]

    lambda = Dict(zip(zip(rows.state, rows.period), rows.lambda))
    n = Dict(zip(zip(rows.state, rows.period), rows.n))
    periods = sort(unique(rows.period))
    controls_at = control_states_of(plot, unique(rows.state))

    series = empty_series()
    vlines = DataFrame(panel = Int[], x = Int[])
    for (panel, s) in enumerate(treated)
        if controls == "average"
            mean_periods, mean_values = average_of_controls(lambda, n, period -> controls_at(s, period), periods, weights)
            add_line!(series, panel, "Control states", control_mean_label(weights), :control_mean, 0, mean_periods, mean_values)
        else
            control_periods = Dict{String, Vector{Int}}()
            for period in periods, c in controls_at(s, period)
                haskey(lambda, (c, period)) && push!(get!(control_periods, c, Int[]), period)
            end
            for c in sorted_states(collect(keys(control_periods)))
                add_line!(series, panel, c, "Control states", :control, 0, control_periods[c],
                          [lambda[(c, period)] for period in control_periods[c]])
            end
        end
        own_periods = [period for period in periods if haskey(lambda, (s, period))]
        add_line!(series, panel, s, "Treated state", :treated, 1, own_periods, [lambda[(s, period)] for period in own_periods])
        push!(vlines, (panel = panel, x = specs.treatment_periods[findfirst(==(s), specs.treated_states)]))
    end
    return PlotLayout("conditional", "Conditional parallel trends by treated state, ccc = \"$option\"", "Time", outcome_label(plot, options),
                      "", copy(treated), series, vlines, period_ticks(rows), [:vline => "First treated period"])
end


###########################################################################################
# Section 6: Integrated parallel trends layout                                            #
###########################################################################################

function integrated_layout(plot::DiDIntPlot, treated, ccc, controls, weights)

    # One panel per treated state s (Karim & Webb 2025, eqs. 36-38): the change since the
    # previous period in the predicted outcome of s and of its control states, all evaluated at
    # the covariate means of s. The controls are drawn as their average in each period,
    # weighted by their numbers of observations in the two periods (the weights of didint()) or
    # equally, or as one line per control state.
    if !isnothing(ccc) && unique(lowercase.(ccc isa AbstractVector ? ccc : [ccc])) != ["int"]
        error("The integrated parallel trends plot is that of ccc = \"int\"; 'ccc' cannot select another option.")
    end
    data = plot.data
    rows = data[.!ismissing.(data.delta), :]
    series = empty_series()
    vlines = DataFrame(panel = Int[], x = Int[])
    panels = String[]
    for s in treated
        state_rows = sort(rows[rows.treated_state_covars .== s, :], :period)
        own = state_rows[state_rows.state .== s, :]
        others = state_rows[state_rows.state .!= s, :]
        nrow(own) == 0 && continue
        push!(panels, s)
        panel = length(panels)
        if controls == "average"
            periods = sort(unique(others.period))
            means = [(in_period = others.period .== period;
                      weighted_mean(others.delta[in_period], weights ? Float64.(others.n[in_period]) : ones(count(in_period))))
                     for period in periods]
            add_line!(series, panel, "Control states", control_mean_label(weights), :control_mean, 0, periods, means)
        else
            for c in sorted_states(unique(others.state))
                control_rows = others[others.state .== c, :]
                add_line!(series, panel, c, "Control states", :control, 0, control_rows.period, control_rows.delta)
            end
        end
        add_line!(series, panel, s, "Treated state", :treated, 1, own.period, own.delta)
        push!(vlines, (panel = panel, x = first(own.treat_period)))
    end
    isempty(panels) && error("No integrated parallel trends data for the selected treated states.")
    return PlotLayout("integrated", "Integrated parallel trends by treated state",
                      "Time", "Change in predicted " * plot.specs.outcome * " since the previous period", "", panels, series, vlines,
                      period_ticks(rows), [:vline => "First treated period"])
end
