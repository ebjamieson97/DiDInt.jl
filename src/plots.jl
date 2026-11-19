###########################################################################################
# plots.jl: the object that didint_plot() returns and the layouts of its plots            #
#                                                                                         #
#   Section 1: DiDIntPlot                                                                 #
#   Section 2: Display                                                                    #
#   Section 3: Plot layouts: the arguments of plot() and what each plot shows             #
#   Section 4: Event study layout                                                         #
#   Section 5: Conditional parallel trends layouts                                        #
#   Section 6: Layouts of changes: conditional2, integrated and integrated2               #
#                                                                                         #
# The layouts are plain tables of lines and reference marks. They are drawn by the Makie  #
# extension (ext/DiDIntMakieExt.jl), which is loaded together with a Makie backend.       #
###########################################################################################


###########################################################################################
# Section 1: DiDIntPlot                                                                   #
###########################################################################################

"""
    DiDIntPlot

The data that [`didint_plot`](@ref) prepares, which are drawn with `plot`. Displaying it prints a summary.

# Fields
One table (a DataFrame) for each kind of plot, described under "Returns" in [`didint_plot`](@ref). A table is empty
if that kind of plot data was not prepared.

- `conditional`: the mean outcome of each state and period, adjusted for covariates (`lambda`), by `ccc` option.
- `conditional2`: for each `ccc` option and treated state, the change in `lambda` of that state and of its control
  states since the period before its treatment (`delta`).
- `integrated`: for each treated state, the change in the predicted outcome of that state and of its control
  states since the previous period, at the treated state's covariate means (`delta`).
- `integrated2`: the same change measured from the period before the treated state's treatment.
- `event`: the event study estimates by `ccc` option and period relative to treatment.
- `event_lambdas`: the `lambda` of each treated state behind `event`, by `ccc` option and period relative to
  treatment (`time_since_treatment`), with the state's number of observations (`n`).
- `specs::NamedTuple`: the specification (`outcome`, `plot_types`, `ccc`, `covariates`, `treated_states` with their
  `treatment_times` and `treatment_periods`, `notyet`, `weights`, `ci`, the period length and date range, and so on).

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
    One of the kinds of plot data that `x` holds:
    - `"event"`: the event study, one panel per `ccc` option.
    - `"conditional"`: the conditional parallel trends, the mean outcome (adjusted for covariates) by period
      (see `facet`).
    - `"conditional2"`: the conditional parallel trends as changes since the pre-treatment period, one panel per
      treated state: the change in its adjusted mean outcome since the period before its treatment, next to the
      average change of its control states.
    - `"integrated"`: the integrated parallel trends of `ccc = "int"`, one panel per treated state, as changes
      since the previous period.
    - `"integrated2"`: the same as changes since the pre-treatment period.

    In `"conditional2"` and `"integrated2"` both lines are zero in the pre-treatment period. The gap between them
    is the treated state's ATT in each treated period, as `didint()` estimates it, and a placebo estimate in each
    period before treatment, which parallel trends requires to be close to zero.

    By default: the conditional parallel trends if `x` has them, otherwise the first kind that it has.
- `treated_states = nothing`
    The treated state(s) to show, as a single value or a vector. By default all of them.
- `ccc = nothing`
    The `ccc` option(s) to show in an event study or a conditional parallel trends plot, as a string or a vector
    of strings. By default all of those in the data.
- `window = nothing`
    A pair `(first, last)` that limits the horizontal axis: periods relative to treatment for an event study,
    period numbers (the `period` column) otherwise.

## Parallel trends plots
- `facet::AbstractString = "ccc"`
    For `plot_type = "conditional"`: `"ccc"` draws one panel per `ccc` option with a line for every state;
    `"state"` draws one panel per treated state and `ccc` option (a row of panels per treated state),
    comparing the treated state with its control states. `"conditional2"` is always drawn that way. `ccc`
    selects the options in both.
- `controls::AbstractString = "average"`
    In the panels of a treated state (every plot type except `"event"` and `"conditional"` with
    `facet = "ccc"`): `"average"` draws one line for the average of the control states, `"separate"` a line
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
- `ncols::Union{Integer, Nothing} = nothing`: number of columns of panels. By default one per `ccc` option in a
  `"conditional2"` plot of several options, otherwise about the square root of the number of panels.
- `size::Union{Tuple, Nothing} = nothing`: size of the figure, as `(width, height)`.
- `title::Union{AbstractString, Nothing} = nothing`: title of the figure; `""` for none.
- `theme`: a Makie `Theme` that replaces the default look.

The control states of a treated state are the never-treated states and, if `didint_plot()` was called with
`notyet = true`, the states of later cohorts while they are not yet treated, as in `didint()`. Panels of treated
states are ordered from the earliest treated state to the latest. A dashed vertical line marks the first treated
period.
"""
struct DiDIntPlot
    conditional::DataFrame
    conditional2::DataFrame
    integrated::DataFrame
    integrated2::DataFrame
    event::DataFrame
    event_lambdas::DataFrame
    specs::NamedTuple
end


###########################################################################################
# Section 2: Display                                                                      #
###########################################################################################

# The kinds of plot data of didint_plot() and of plots of plot(), each with its table in a
# DiDIntPlot
const PLOT_TYPES = ["event", "conditional", "conditional2", "integrated", "integrated2"]

# The ccc options in the order of their panels, and the titles of the panels
const PLOT_CCC_ORDER = ["none", "hom", "time", "state", "add", "int"]
const PLOT_CCC_LABELS = Dict("none" => "No covariates (ccc = \"none\")", "hom" => "Homogeneous (ccc = \"hom\")",
                             "time" => "Time-varying (ccc = \"time\")", "state" => "Region-varying (ccc = \"state\")",
                             "add" => "Two one-way (ccc = \"add\")", "int" => "Two-way (ccc = \"int\")")

# The table of one kind of plot data
plot_table(plot::DiDIntPlot, plot_type::AbstractString) = getfield(plot, Symbol(plot_type))

# The kinds of plots that the data of a DiDIntPlot allow
available_plot_types(plot::DiDIntPlot) = [type for type in PLOT_TYPES if nrow(plot_table(plot, type)) > 0]

function Base.show(io::IO, plot::DiDIntPlot)

    # One-line form, used inside collections
    print(io, "DiDIntPlot(plot data for: ", join(available_plot_types(plot), ", "), ")")
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

    # Summary: the specification and the tables of plot data
    specs = plot.specs
    types = available_plot_types(plot)
    println(io, "DID-INT plot data")
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
        println(io, "  Controls:         ", specs.notyet ? "never-treated and not-yet-treated states" : "never-treated states")
    end
    println(io, "  Periods:          ", specs.period, ", from ", specs.start_date, " to ", specs.end_date)
    println(io)
    println(io, "Tables of plot data, each drawn with plot(x; plot_type = ...) (see ?DiDIntPlot)")
    width = maximum(length, types; init = 0)
    for (i, type) in enumerate(types)
        table = plot_table(plot, type)
        print(io, "  ", rpad(type, width), "   ", nrow(table), " rows: ", join(names(table), ", "))
        i < length(types) && println(io)
    end
end


###########################################################################################
# Section 3: Plot layouts: the arguments of plot() and what each plot shows               #
###########################################################################################

# What a plot shows, independently of the plotting package that draws it.
#   title, subtitle: the title of the figure and a second line under it (may be empty)
#   panels:  the titles of the panels, in order
#   series:  one row per point, with the panel (its number), the line the point belongs to
#            (`group`), the entry of the legend that the line belongs to (`label`), what the
#            line is (`role`: :estimate, :treated, :state, :control or :control_mean), the
#            number of its colour among the coloured lines (`color`, 0 for control lines; a
#            treated state keeps its number whichever treated states are shown), its
#            coordinates (`x`, `y`) and, for an event study, the confidence bounds
#   vlines:  the vertical reference lines (first treated periods), by panel
#   xticks:  for a period axis, the period numbers and the labels of their dates
#   legend_extras: entries of the legend for marks other than the lines, as kind => label,
#            where the kind is :band (confidence interval), :no_interval (estimate without
#            one) or :vline (vertical reference line)
#   ncols:   the number of columns of panels that suits the plot (a grid of treated states by
#            ccc options has one column per option), or nothing to let the drawing decide
struct PlotLayout
    plot_type::String
    title::String
    subtitle::String
    xlabel::String
    ylabel::String
    legend_title::String
    panels::Vector{String}
    series::DataFrame
    vlines::DataFrame
    xticks::Union{Nothing, Tuple{Vector{Int}, Vector{String}}}
    legend_extras::Vector{Pair{Symbol, String}}
    ncols::Union{Nothing, Int}
end

PlotLayout(plot_type, title, subtitle, xlabel, ylabel, legend_title, panels, series, vlines, xticks, legend_extras) =
    PlotLayout(plot_type, title, subtitle, xlabel, ylabel, legend_title, panels, series, vlines, xticks, legend_extras, nothing)

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
    elseif plot_type == "conditional"
        layout = facet == "ccc" ? conditional_layout(plot, treated, ccc) :
                                  conditional_state_layout(plot, treated, ccc, controls, something(weights, true))
    else
        layout = changes_layout(plot, plot_type, treated, ccc, controls, something(weights, true))
    end
    layout = apply_window(layout, window)
    return isnothing(title) ? layout : with_title(layout, title)
end

# What each kind of plot data needs from didint_plot(), for the message when it is missing
const PLOT_TYPE_NEEDS = Dict("event" => "the treated states",
                             "conditional" => "a 'ccc' option other than \"int\"",
                             "conditional2" => "the treated states and a 'ccc' option other than \"int\"",
                             "integrated" => "the treated states and \"int\" among the 'ccc' options",
                             "integrated2" => "the treated states and \"int\" among the 'ccc' options")

function check_plot_type(plot::DiDIntPlot, plot_type)

    # The plot type asked for, or the default one for the data: the conditional parallel trends
    # if the data have them, otherwise the first kind that they have
    available = available_plot_types(plot)
    isempty(available) && error("The plot data are empty.")
    isnothing(plot_type) && return "conditional" in available ? "conditional" : first(available)
    plot_type = lowercase(strip(plot_type))
    if !(plot_type in PLOT_TYPES)
        error("'plot_type' must be one of: ", join(("\"" * type * "\"" for type in PLOT_TYPES), ", "), ".")
    end
    plot_type in available && return plot_type
    error("The data for plot_type = \"$plot_type\" were not prepared. Call didint_plot() with plot_type = \"all\" (the default) or " *
          "including \"$plot_type\"; these data need ", PLOT_TYPE_NEEDS[plot_type], ".")
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

function grid_panel_title(plot::DiDIntPlot, s, option, n_options)

    # Title of the panel of treated state s under a ccc option, in a plot of n_options options
    # (the option is named only when there are several; `option` is nothing for "int")
    state_title = only(treated_state_titles(plot, [s]))
    return n_options == 1 ? state_title : state_title * ", " * ccc_subtitle(option)
end

# A grid of treated states by ccc options has a column per option
grid_columns(options) = length(options) == 1 ? nothing : length(options)

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
    return PlotLayout(layout.plot_type, layout.title, layout.subtitle, layout.xlabel, layout.ylabel, layout.legend_title,
                      layout.panels, series, vlines, xticks, legend_extras, layout.ncols)
end

# A title given by the user replaces the title and the line under it
with_title(layout::PlotLayout, title) = PlotLayout(layout.plot_type, String(title), "", layout.xlabel, layout.ylabel,
                                                   layout.legend_title, layout.panels, layout.series, layout.vlines, layout.xticks,
                                                   layout.legend_extras, layout.ncols)

# Whether the outcome shown is adjusted for covariates: not without covariates, nor when only
# ccc = "none" is shown
is_adjusted(plot::DiDIntPlot, options) = !isempty(plot.specs.covariates) && any(!=("none"), options)

# Label of a vertical axis that shows the outcome (`prefix` precedes its name)
outcome_label(plot::DiDIntPlot, options; prefix = "") = prefix * plot.specs.outcome * (is_adjusted(plot, options) ? " (covariate-adjusted)" : "")

# The line under the title of a plot that shows a single ccc option
ccc_subtitle(option) = "ccc = \"" * option * "\""

sorted_states(states) = sort(String.(states), by = custom_sort_order)

function period_ticks(plot::DiDIntPlot, periods)

    # The period numbers that are plotted, with the labels of their dates
    periods = sort(unique(Int.(periods)))
    return (periods, String[plot.specs.period_labels[period + 1] for period in periods])
end

function by_treatment_time(plot::DiDIntPlot, treated)

    # The treated states from the earliest treated to the latest (those treated at the same
    # time in the order of their names), the order of their panels
    period = Dict(zip(plot.specs.treated_states, plot.specs.treatment_periods))
    return sort(treated, by = s -> (period[s], custom_sort_order(s)))
end

function treated_state_titles(plot::DiDIntPlot, treated)

    # Titles of the panels of treated states: the state and when it is first treated
    time = Dict(zip(plot.specs.treated_states, plot.specs.treatment_times))
    return [string(s, " (treated ", time[s], ")") for s in treated]
end

function weighted_mean(values, weights)
    return sum(values .* weights) / sum(weights)
end

function control_mean_label(weights)
    return weights ? "Control states (mean, weighted by observations)" : "Control states (mean)"
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
        table = plot.event
    else
        table = event_study_data(plot.event_lambdas[in.(plot.event_lambdas.state, Ref(treated)), :], weights, specs.ci, specs.hc)
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
    return PlotLayout("event", "Event study", "", "Periods relative to treatment", outcome_label(plot, options), "",
                      [PLOT_CCC_LABELS[option] for option in options], series, vlines, nothing, legend_extras)
end


###########################################################################################
# Section 5: Conditional parallel trends layouts                                          #
###########################################################################################

function conditional_title(plot::DiDIntPlot, options)

    # Title of the plots of the lambdas themselves
    return "Conditional parallel trends: " * plot.specs.outcome * " means" * (is_adjusted(plot, options) ? " (covariate-adjusted)" : "") *
           " by period"
end

function conditional_layout(plot::DiDIntPlot, treated, ccc)

    # One panel per ccc option, with the lambdas of every state as a line (Karim & Webb 2025,
    # eq. 45). The selected treated states are each in their own colour and the states that
    # are never treated are drawn as control states; treated states that are not selected, and
    # states treated throughout, are left out. When the treated states are not known, every
    # state is drawn alike. Vertical lines mark the first treated periods.
    specs = plot.specs
    rows = plot.conditional
    options = selected_ccc(unique(rows.ccc), ccc)
    known = !isempty(specs.treated_states)

    # A colour for each treated state, which it keeps when only some of them are selected
    color = Dict(s => i for (i, s) in enumerate(sorted_states(specs.treated_states)))

    series = empty_series()
    for (panel, option) in enumerate(options)
        option_rows = rows[rows.ccc .== option, :]
        for state_rows in groupby(option_rows, :state; sort = false)
            state = first(state_rows.state)
            state_rows = sort(state_rows, :period)
            if !known
                add_line!(series, panel, state, "States", :state, 1, state_rows.period, state_rows.lambda)
            elseif state in treated
                add_line!(series, panel, state, state, :treated, color[state], state_rows.period, state_rows.lambda)
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
    return PlotLayout("conditional", conditional_title(plot, options), "", "Time", outcome_label(plot, options),
                      known ? "Treated state" : "", [PLOT_CCC_LABELS[option] for option in options], series, vlines,
                      period_ticks(plot, series.x), [:vline => (length(treatment_periods) == 1 ? "First treated period" : "First treated periods")])
end

function control_states_of(plot::DiDIntPlot, states)

    # The control states of a treated state s, as a function of the period: among `states`, the
    # never-treated states and, with `notyet`, the states of later cohorts that are not yet
    # treated in that period, as in didint() (see control_states_function). States treated
    # throughout are never controls.
    specs = plot.specs
    treatment_period = Dict(zip(specs.treated_states, specs.treatment_periods))
    never_treated = sorted_states(setdiff(states, specs.treated_states, specs.always_treated))
    later = sorted_states(intersect(specs.treated_states, states))
    return function (s, period)
        specs.notyet || return never_treated
        return vcat(never_treated, [j for j in later if treatment_period[j] > max(period, treatment_period[s])])
    end
end

function average_of_controls(lambda, n, controls_at, periods, weights)

    # The line of the average control state, built from changes, as didint() compares the
    # change of a treated state with the average of its controls' changes. It starts at the
    # controls' average lambda in the first period, and each later period adds the average of
    # the controls' changes since the previous period: over the controls of that period
    # (`controls_at(period)`) observed in both periods, weighted by their numbers of
    # observations in the two periods or equally. With the same controls and weights in every
    # period this is the average of the controls' lambdas. Unlike an average of the lambdas,
    # it does not move when the set of controls or their weights change, nor with shifts of
    # the lambdas that are constant within a state.
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

function conditional_state_layout(plot::DiDIntPlot, treated, ccc, controls, weights)

    # One panel per treated state (the earliest treated first) and ccc option, a row of panels
    # per treated state: the lambdas of the treated state and of its control states, the latter
    # as one line per control state or as their average (see average_of_controls). A state
    # that is not yet treated is a control up to the period before its treatment.
    specs = plot.specs
    isempty(specs.treated_states) && error("facet = \"state\" needs the treated states: call didint_plot() with 'treated_states' or 'gvar'.")
    options = selected_ccc(unique(plot.conditional.ccc), ccc)
    treatment_period = Dict(zip(specs.treated_states, specs.treatment_periods))
    treated = by_treatment_time(plot, treated)

    series = empty_series()
    vlines = DataFrame(panel = Int[], x = Int[])
    panel_titles = String[]
    # the lambdas, numbers of observations, periods and control states under each option
    by_option = Dict(map(options) do option
        rows = plot.conditional[plot.conditional.ccc .== option, :]
        option => (lambda = Dict(zip(zip(rows.state, rows.period), rows.lambda)), n = Dict(zip(zip(rows.state, rows.period), rows.n)),
                   periods = sort(unique(rows.period)), controls_at = control_states_of(plot, unique(rows.state)))
    end)
    for s in treated, option in options
        (; lambda, n, periods, controls_at) = by_option[option]
        controls_of_s(period) = controls_at(s, period)
        push!(panel_titles, grid_panel_title(plot, s, option, length(options)))
        panel = length(panel_titles)
        if controls == "average"
            mean_periods, mean_values = average_of_controls(lambda, n, controls_of_s, periods, weights)
            add_line!(series, panel, "Control states", control_mean_label(weights), :control_mean, 0, mean_periods, mean_values)
        else
            control_periods = Dict{String, Vector{Int}}()
            for period in periods, c in controls_of_s(period)
                haskey(lambda, (c, period)) && push!(get!(control_periods, c, Int[]), period)
            end
            for c in sorted_states(collect(keys(control_periods)))
                add_line!(series, panel, c, "Control states", :control, 0, control_periods[c],
                          [lambda[(c, period)] for period in control_periods[c]])
            end
        end
        own_periods = [period for period in periods if haskey(lambda, (s, period))]
        add_line!(series, panel, s, "Treated state", :treated, 1, own_periods, [lambda[(s, period)] for period in own_periods])
        push!(vlines, (panel = panel, x = treatment_period[s]))
    end
    subtitle = length(options) == 1 ? ccc_subtitle(only(options)) : ""
    return PlotLayout("conditional", conditional_title(plot, options), subtitle, "Time", outcome_label(plot, options), "",
                      panel_titles, series, vlines, period_ticks(plot, series.x), [:vline => "First treated period"],
                      grid_columns(options))
end


###########################################################################################
# Section 6: Layouts of changes: conditional2, integrated and integrated2                 #
###########################################################################################

function with_base_period(periods, values, base)

    # The points of a line of changes since the period `base`, with the point of `base` itself,
    # where the change is zero
    position = searchsortedfirst(periods, base)
    return vcat(periods[1:position - 1], base, periods[position:end]), vcat(values[1:position - 1], 0.0, values[position:end])
end

function changes_layout(plot::DiDIntPlot, plot_type, treated, ccc, controls, weights)

    # One panel per treated state s, the earliest treated first: the change in the outcome of s
    # next to the changes of its control states, from the table of the plot type.
    #   "conditional2": the change in lambda since the period before the treatment of s, under
    #                   the ccc options other than "int"; with several options, one panel per
    #                   treated state and option, a row of panels per treated state
    #   "integrated":   the change in the predicted outcome since the previous period, at the
    #                   covariate means of s (Karim & Webb 2025, eqs. 36-38)
    #   "integrated2":  the same since the period before the treatment of s
    # In the two plots of changes since the period before treatment, the change is zero in
    # that period for every state, and the treated state's change minus the average of its
    # controls is ATT(s,t) in every treated period t and a placebo estimate before treatment.
    # The controls are drawn as their average in each period, weighted by their numbers of
    # observations in the two periods (the weights of didint()) or equally, or as one line per
    # control state.
    rows = plot_table(plot, plot_type)
    long = plot_type != "integrated"
    outcome = plot.specs.outcome
    if plot_type == "conditional2"
        options = selected_ccc(unique(rows.ccc), ccc)
        treated_column = :treated_state
        title = "Conditional parallel trends: changes since the pre-treatment period"
        subtitle = length(options) == 1 ? ccc_subtitle(only(options)) : ""
        ylabel = outcome_label(plot, options; prefix = "Change in ")
    else
        if !isnothing(ccc) && unique(lowercase.(ccc isa AbstractVector ? ccc : [ccc])) != ["int"]
            error("The integrated parallel trends plot is that of ccc = \"int\"; 'ccc' cannot select another option.")
        end
        options = [nothing]
        treated_column = :treated_state_covars
        title = "Integrated parallel trends: changes since the " * (long ? "pre-treatment period" : "previous period")
        subtitle, ylabel = "", "Change in predicted " * outcome
    end

    series = empty_series()
    vlines = DataFrame(panel = Int[], x = Int[])
    panel_titles = String[]
    for s in by_treatment_time(plot, treated), option in options
        option_rows = isnothing(option) ? rows : rows[rows.ccc .== option, :]
        state_rows = sort(option_rows[option_rows[!, treated_column] .== s, :], :period)
        own = state_rows[state_rows.state .== s, :]
        others = state_rows[state_rows.state .!= s, :]
        nrow(own) == 0 && continue
        push!(panel_titles, grid_panel_title(plot, s, option, length(options)))
        panel = length(panel_titles)
        line(periods, values) = long ? with_base_period(Int.(periods), Float64.(values), first(own.base_period)) : (periods, values)
        if controls == "average"
            periods = sort(unique(others.period))
            means = [(in_period = others.period .== period;
                      weighted_mean(others.delta[in_period], weights ? Float64.(others.n[in_period]) : ones(count(in_period))))
                     for period in periods]
            isempty(periods) || add_line!(series, panel, "Control states", control_mean_label(weights), :control_mean, 0, line(periods, means)...)
        else
            for c in sorted_states(unique(others.state))
                control_rows = others[others.state .== c, :]
                add_line!(series, panel, c, "Control states", :control, 0, line(control_rows.period, control_rows.delta)...)
            end
        end
        add_line!(series, panel, s, "Treated state", :treated, 1, line(own.period, own.delta)...)
        push!(vlines, (panel = panel, x = first(own.treat_period)))
    end
    isempty(panel_titles) && error("No data of plot_type = \"$plot_type\" for the selected treated states.")
    return PlotLayout(plot_type, title, subtitle, "Time", ylabel, "", panel_titles, series, vlines,
                      period_ticks(plot, series.x), [:vline => "First treated period"],
                      grid_columns(options))
end
