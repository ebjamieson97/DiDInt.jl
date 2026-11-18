###########################################################################################
# DiDIntMakieExt.jl: drawing a DiDIntPlot with Makie                                      #
#                                                                                         #
#   Section 1: plot(), the method for a DiDIntPlot                                        #
#   Section 2: Colours and theme                                                          #
#   Section 3: Drawing a layout: panels, lines, reference lines and legend                #
#   Section 4: Axis ticks                                                                 #
#                                                                                         #
# This extension is loaded when Makie is (through a backend such as CairoMakie). What a   #
# plot shows is decided in src/plots.jl (plot_layout); here it is only drawn.             #
###########################################################################################

module DiDIntMakieExt

using DiDInt
using DiDInt: DiDIntPlot, PlotLayout, plot_layout
using Makie


###########################################################################################
# Section 1: plot(), the method for a DiDIntPlot                                          #
###########################################################################################

"""
    plot(x::DiDIntPlot; plot_type, treated_states, ccc, window,
         facet, controls, weights, groupmin, ncols, size, title, theme)

Draw one of the plots of the data prepared by `didint_plot()` (event study, conditional parallel trends or integrated
parallel trends), as a Makie `Figure`. The arguments are described under [`DiDIntPlot`](@ref).
"""
function Makie.plot(x::DiDIntPlot;
                    plot_type::Union{AbstractString, Nothing} = nothing,
                    treated_states = nothing,
                    ccc::Union{AbstractString, AbstractVector{<:AbstractString}, Nothing} = nothing,
                    facet::AbstractString = "ccc",
                    controls::AbstractString = "average",
                    weights::Union{Bool, Nothing} = nothing,
                    window = nothing,
                    groupmin::Integer = 3,
                    ncols::Union{Integer, Nothing} = nothing,
                    size::Union{Tuple{<:Real, <:Real}, Nothing} = nothing,
                    title::Union{AbstractString, Nothing} = nothing,
                    theme::Makie.Theme = didint_theme())

    layout = plot_layout(x; plot_type = plot_type, treated_states = treated_states, ccc = ccc, facet = facet, controls = controls,
                         weights = weights, window = window, groupmin = groupmin, title = title)
    if !isnothing(ncols) && ncols < 1
        error("'ncols' must be at least 1.")
    end
    return Makie.with_theme(() -> draw_layout(layout, ncols, size), theme)
end


###########################################################################################
# Section 2: Colours and theme                                                            #
###########################################################################################

# Colours of the lines that are told apart by colour (treated states, estimates), in a fixed
# order; of the event study estimates and their confidence band; of the control states' lines
# and of their average; of the vertical line at the first treated period; and of text, axes
# and grid
const SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948",
                       "#8c5a2b", "#12a5c9", "#9a9c00", "#a23b8f", "#0d366b", "#7ac74f", "#7b1f1f", "#00695c"]
const SERIES_LINESTYLES = [:solid, :dash, :dashdot]
const ESTIMATE_COLOR = "#1f1f1d"
const ESTIMATE_BAND_COLOR = "#2a78d6"
const CONTROL_COLOR = "#b9b8b0"
const CONTROL_MEAN_COLOR = "#52514e"
const TREATMENT_LINE = (color = "#1f1f1d", linestyle = :dash, linewidth = 1.7)
const TEXT_COLOR = "#0b0b0b"
const SECONDARY_TEXT_COLOR = "#52514e"
const AXIS_COLOR = "#898781"
const GRID_COLOR = "#e6e5df"
const STRIP_COLOR = "#ecebe6"

# The colour of the coloured line number `index`, and its line style: when there are more
# lines than colours, the colours are used again with another style
series_color(index) = SERIES_COLORS[mod1(index, length(SERIES_COLORS))]
series_linestyle(index) = SERIES_LINESTYLES[mod1(cld(index, length(SERIES_COLORS)), length(SERIES_LINESTYLES))]

# The default look: white panels with a light grid and a frame, and a shaded strip with the
# title above each panel
function didint_theme()
    return Makie.Theme(
        fontsize = 13,
        textcolor = TEXT_COLOR,
        Axis = (backgroundcolor = :white,
                xgridcolor = GRID_COLOR, ygridcolor = GRID_COLOR, xgridwidth = 1, ygridwidth = 1,
                topspinecolor = AXIS_COLOR, bottomspinecolor = AXIS_COLOR, leftspinecolor = AXIS_COLOR, rightspinecolor = AXIS_COLOR,
                xtickcolor = AXIS_COLOR, ytickcolor = AXIS_COLOR, xticksize = 4, yticksize = 4,
                xticklabelcolor = SECONDARY_TEXT_COLOR, yticklabelcolor = SECONDARY_TEXT_COLOR,
                xticklabelsize = 12, yticklabelsize = 12),
        Legend = (framevisible = false, labelsize = 12, titlesize = 12, patchsize = (26, 14), padding = (0, 0, 0, 0)))
end


###########################################################################################
# Section 3: Drawing a layout: panels, lines, reference lines and legend                  #
###########################################################################################

# Lines with at most this many points are drawn with markers
const MAX_POINTS_WITH_MARKERS = 30

# Size of the title, and the width of a character of it relative to that size (an estimate
# for the bold font, used to make the title fit the figure)
const TITLE_FONTSIZE = 16
const TITLE_CHARACTER_WIDTH = 0.6

# A legend with more entries than this is not drawn, and the number of entries on a row of it
const MAX_LEGEND_ENTRIES = 40
const LEGEND_ENTRIES_PER_ROW = 7

# Estimated width of an entry of the legend: its mark and gaps, and each character of its label
const LEGEND_MARK_WIDTH = 42
const LEGEND_CHARACTER_WIDTH = 6.0

function draw_layout(layout::PlotLayout, ncols, size)

    # The figure: a title (with a second line under it), a grid of panels that share their axes (one panel per element of
    # `layout.panels`, filled row by row), the axis labels and a legend below the panels
    n_panels = length(layout.panels)
    ncols = isnothing(ncols) ? (n_panels <= 3 ? n_panels : ceil(Int, sqrt(n_panels))) : min(ncols, n_panels)
    nrows = cld(n_panels, ncols)
    with_markers = draws_markers(layout)
    elements, labels = legend_contents(layout, with_markers)
    has_legend = 2 <= length(labels) <= MAX_LEGEND_ENTRIES
    # By default the figure is wide enough for its panels and for its title; in a figure of a
    # given size, a title that is too long is drawn smaller
    # The legend has as many rows as its entries need to fit the width of the figure.
    title_width = TITLE_CHARACTER_WIDTH * TITLE_FONTSIZE * length(layout.title)
    width = isnothing(size) ? max(560, 300 * ncols + 80, ceil(Int, title_width) + 60) : size[1]
    legend_rows = legend_rows_for(labels, layout.legend_title, width - 40)
    if isnothing(size)
        size = (width, 250 * nrows + 90 + (isempty(layout.title) ? 0 : 30) + (isempty(layout.subtitle) ? 0 : 20) +
                       (has_legend ? 18 + 22 * legend_rows : 0))
    end
    title_fontsize = title_width <= size[1] - 60 ? TITLE_FONTSIZE : max(9, TITLE_FONTSIZE * (size[1] - 60) / title_width)

    figure = Makie.Figure(; size = size)
    panels = figure[2, 2] = Makie.GridLayout()
    if !isempty(layout.title) || !isempty(layout.subtitle)
        # the title and the line under it, centred above the panels
        header = figure[1, 2] = Makie.GridLayout()
        isempty(layout.title) || Makie.Label(header[1, 1], layout.title; font = :bold, fontsize = title_fontsize, tellwidth = false)
        isempty(layout.subtitle) || Makie.Label(header[2, 1], layout.subtitle; fontsize = 13, color = SECONDARY_TEXT_COLOR, tellwidth = false)
        Makie.rowgap!(header, 4)
    end
    Makie.Label(figure[2, 1], layout.ylabel; rotation = pi / 2, tellheight = false)
    Makie.Label(figure[3, 2], layout.xlabel; tellwidth = false)

    panel_width = (size[1] - 80) / ncols
    xticks = axis_ticks(layout, panel_width)
    axes = Makie.Axis[]
    for (panel, panel_title) in enumerate(layout.panels)
        row, column = fldmod1(panel, ncols)
        axis = Makie.Axis(panels[row, column]; xticks = xticks)
        push!(axes, axis)

        # The title of the panel, in a shaded strip above it
        Makie.Box(panels[row, column, Makie.Top()]; color = STRIP_COLOR, strokecolor = AXIS_COLOR, strokewidth = 1)
        Makie.Label(panels[row, column, Makie.Top()], panel_title; font = :bold, fontsize = 12, padding = (4, 4, 5, 5))

        draw_panel!(axis, layout, panel, with_markers)

        # Tick labels on the outside of the grid only: at the left of the first column, and
        # under the panels that have no panel below them
        column == 1 || Makie.hideydecorations!(axis; grid = false, ticks = false)
        panel + ncols <= n_panels && Makie.hidexdecorations!(axis; grid = false, ticks = false)
    end
    Makie.linkaxes!(axes...)
    Makie.colgap!(panels, 12)
    Makie.rowgap!(panels, 12)

    if has_legend
        if isempty(layout.legend_title)
            Makie.Legend(figure[4, 2], elements, labels; orientation = :horizontal, nbanks = legend_rows, tellwidth = false)
        else
            Makie.Legend(figure[4, 2], elements, labels, layout.legend_title; orientation = :horizontal, titleposition = :left,
                         nbanks = legend_rows, tellwidth = false)
        end
    end
    return figure
end

function legend_rows_for(labels, title, width)

    # The number of rows of the legend: the smallest number, with at most
    # LEGEND_ENTRIES_PER_ROW entries on a row, for which the legend is estimated to fit in
    # `width`. The entries fill the legend column by column.
    entry_width(label) = LEGEND_MARK_WIDTH + LEGEND_CHARACTER_WIDTH * length(label)
    for rows in cld(length(labels), LEGEND_ENTRIES_PER_ROW):length(labels)
        columns = Iterators.partition(labels, rows)
        sum(maximum(entry_width, column) for column in columns) + LEGEND_CHARACTER_WIDTH * length(title) <= width && return rows
    end
    return max(length(labels), 1)
end

function draws_markers(layout::PlotLayout)

    # Markers are drawn on the lines of estimates, and on those of a treated state and of the
    # average control state when a panel compares the two, unless the lines have many points
    points = Dict{Tuple{Int, String}, Int}()
    for row in eachrow(layout.series)
        row.role in (:estimate, :treated, :control_mean) || continue
        points[(row.panel, row.group)] = get(points, (row.panel, row.group), 0) + 1
    end
    isempty(points) && return false
    lines_in_panel = Dict{Int, Int}()
    for (panel, _) in keys(points)
        lines_in_panel[panel] = get(lines_in_panel, panel, 0) + 1
    end
    return maximum(values(lines_in_panel)) <= 2 && maximum(values(points)) <= MAX_POINTS_WITH_MARKERS
end

function draw_panel!(axis, layout::PlotLayout, panel, with_markers)

    # The lines of one panel, from back to front: the reference lines, the control states,
    # their average, then the treated states or the estimates
    for x in layout.vlines.x[layout.vlines.panel .== panel]
        Makie.vlines!(axis, x; TREATMENT_LINE...)
    end
    series = layout.series[layout.series.panel .== panel, :]
    for role in (:control, :state, :control_mean, :estimate, :treated)
        role_rows = series[series.role .== role, :]
        for group in unique(role_rows.group)
            rows = role_rows[role_rows.group .== group, :]
            draw_line!(axis, role, first(rows.color), rows.x, rows.y, rows.lower, rows.upper, with_markers)
        end
    end
    return axis
end

function draw_line!(axis, role, color, x, y, lower, upper, with_markers)

    # One line, drawn according to what it is
    if role == :control
        Makie.lines!(axis, x, y; color = (CONTROL_COLOR, 0.8), linewidth = 1)
    elseif role == :state
        Makie.lines!(axis, x, y; color = (series_color(color), 0.6), linewidth = 1.2)
    elseif role == :control_mean
        Makie.lines!(axis, x, y; color = CONTROL_MEAN_COLOR, linewidth = 2)
        with_markers && Makie.scatter!(axis, x, y; color = CONTROL_MEAN_COLOR, markersize = 7)
    elseif role == :treated
        Makie.lines!(axis, x, y; color = series_color(color), linestyle = series_linestyle(color), linewidth = 2)
        with_markers && Makie.scatter!(axis, x, y; color = series_color(color), markersize = 7)
    else
        # An estimate with its confidence interval, a black line in a light band: a band over
        # each run of periods that have an interval, and a square instead of a circle where
        # there is none
        has_interval = .!ismissing.(lower) .& .!ismissing.(upper)
        for run in runs_of(has_interval)
            Makie.band!(axis, x[run], Float64.(lower[run]), Float64.(upper[run]); color = (ESTIMATE_BAND_COLOR, 0.2))
        end
        Makie.lines!(axis, x, y; color = ESTIMATE_COLOR, linewidth = 2)
        any(has_interval) && Makie.scatter!(axis, x[has_interval], y[has_interval]; color = ESTIMATE_COLOR, markersize = 8)
        all(has_interval) || Makie.scatter!(axis, x[.!has_interval], y[.!has_interval]; color = ESTIMATE_COLOR, markersize = 8,
                                            marker = :rect)
    end
end

function runs_of(flags)

    # The ranges of consecutive true entries
    runs = UnitRange{Int}[]
    start = 0
    for (i, flag) in enumerate(flags)
        if flag && start == 0
            start = i
        elseif !flag && start != 0
            push!(runs, start:i - 1)
            start = 0
        end
    end
    start != 0 && push!(runs, start:length(flags))
    return runs
end

function legend_contents(layout::PlotLayout, with_markers)

    # The entries of the legend, as their elements and labels: one per label of the lines, the
    # coloured lines first (in the order of their colours), then the average of the control
    # states, then the control states; then the other marks (confidence band, estimates
    # without an interval, vertical reference lines)
    series = unique(layout.series[!, [:label, :role, :color]])
    order = Dict(:estimate => 1, :treated => 1, :state => 1, :control_mean => 2, :control => 3)
    series = unique(series[sortperm([(order[row.role], row.color) for row in eachrow(series)]), :], :label)
    elements = Any[legend_element(row, with_markers) for row in eachrow(series)]
    labels = String.(series.label)
    for (kind, label) in layout.legend_extras
        if kind == :band
            push!(elements, Makie.PolyElement(color = (ESTIMATE_BAND_COLOR, 0.2), strokewidth = 0))
        elseif kind == :no_interval
            push!(elements, Makie.MarkerElement(color = ESTIMATE_COLOR, marker = :rect, markersize = 8))
        else
            push!(elements, Makie.LineElement(; TREATMENT_LINE...))
        end
        push!(labels, label)
    end
    return elements, labels
end

function legend_element(entry, with_markers)

    # The line (and marker) that stands for an entry in the legend
    if entry.role == :control
        return Makie.LineElement(color = CONTROL_COLOR, linewidth = 1.5)
    end
    color = entry.role == :control_mean ? CONTROL_MEAN_COLOR : entry.role == :estimate ? ESTIMATE_COLOR : series_color(entry.color)
    line = Makie.LineElement(color = color, linestyle = entry.role == :treated ? series_linestyle(entry.color) : :solid, linewidth = 2)
    marker = Makie.MarkerElement(color = color, marker = :circle, markersize = 7)
    return entry.role == :estimate || (with_markers && entry.role != :state) ? [line, marker] : line
end


###########################################################################################
# Section 4: Axis ticks                                                                   #
###########################################################################################

function axis_ticks(layout::PlotLayout, panel_width)

    # Ticks of the horizontal axis that fit in a panel of the given width. For a period axis:
    # the labels of the dates, at every k-th period. For an event study: whole numbers of
    # periods, at multiples of a round step.
    if !isnothing(layout.xticks)
        positions, labels = layout.xticks
        isempty(positions) && return Makie.automatic
        label_width = 7.5 * (maximum(length, labels) + 2)
        step = max(1, cld(length(positions), max(2, floor(Int, panel_width / label_width))))
        keep = 1:step:length(positions)
        return (positions[keep], labels[keep])
    end
    first_x, last_x = extrema(layout.series.x)
    max_ticks = max(2, floor(Int, panel_width / 40))
    step = first(s for s in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 10^9) if (last_x - first_x) / s + 1 <= max_ticks)
    return collect(cld(first_x, step) * step:step:last_x)
end

end
