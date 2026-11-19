# Plotting Examples

`didint_plot()` prepares the data for parallel trends plots and event study plots, and
`plot` draws them once a [Makie](https://docs.makie.org) backend such as CairoMakie is
loaded. This page shows each kind of plot and its options, with the `merit.csv` data of the
[Examples](@ref) page. The formulas behind the plots are in the
[Plotting with `didint_plot()`](@ref) section.

```@setup plots
using CairoMakie
```

```@example plots
using DiDInt, CSV, DataFrames
using CairoMakie

merit = CSV.read(joinpath(@__DIR__, "../../test/data/merit.csv"), DataFrame);
treated_states = [34, 57, 58, 59, 61, 64, 71, 72, 85, 88];
treated_times  = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000];
nothing # hide
```

One call to `didint_plot()` prepares the data for every kind of plot:

```@example plots
plot_data = didint_plot(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    covariates = ["asian", "black", "male"]
)
```

The `plot_type` argument of `plot` chooses the plot:

| `plot_type`      | Plot                                                                                        |
|:-----------------|:--------------------------------------------------------------------------------------------|
| `"event"`        | Event study, one panel per `ccc` option                                                      |
| `"conditional"`  | Conditional parallel trends: the outcome means by period, per `ccc` option or per treated state |
| `"conditional2"` | Conditional parallel trends: changes since the pre-treatment period, per treated state        |
| `"integrated"`   | Integrated parallel trends: changes since the previous period, per treated state              |
| `"integrated2"`  | Integrated parallel trends: changes since the pre-treatment period, per treated state         |

In every plot a dashed vertical line marks the first treated period, the `treated_states`
argument selects the treated states to show, and panels of treated states are ordered
from the earliest treated state to the latest. The pre-treatment period of a treated state
is the period before its first treated period.

The data of each plot are in the table of the same name (`plot_data.conditional`,
`plot_data.integrated2`, and so on). To prepare only some of them, pass their names to the
`plot_type` argument of `didint_plot()`, for example
`plot_type = ["event", "integrated2"]`.

---

## Event Study Plots

The event study plot shows the mean of the treated states' covariate-adjusted outcomes by
period relative to treatment, with its confidence interval, in one panel per `ccc` option.

```@example plots
figure = plot(plot_data; plot_type = "event")
save("plotting_event.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Event study, one panel per ccc option](plotting_event.png)

A confidence interval is drawn where at least three treated states are behind an estimate;
the estimates without one are marked by squares. The `groupmin` argument changes that
number. The `ccc` argument selects the options to show and `window` the periods:

```@example plots
figure = plot(plot_data; plot_type = "event", ccc = ["none", "int"], window = (-6, 4))
save("plotting_event_window.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Event study under two ccc options, six periods before to four after treatment](plotting_event_window.png)

With `treated_states`, the means and intervals are computed for the selected states only,
and `weights = false` weights the treated states equally instead of by their numbers of
observations:

```@example plots
figure = plot(plot_data; plot_type = "event", ccc = "state",
              treated_states = [57, 59, 72, 85], weights = false, groupmin = 2)
save("plotting_event_states.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Event study of four treated states, weighted equally](plotting_event_states.png)

---

## Conditional Parallel Trends Plots

The conditional parallel trends plots are those of the `ccc` options other than `"int"`.
`plot_type = "conditional"`, the default plot, shows the covariate-adjusted mean outcome
of each state and period (the `lambda` column of `plot_data.conditional`).
`plot_type = "conditional2"` shows changes in those means since the pre-treatment period.

### Outcome means, one panel per `ccc` option

By default there is one panel per `ccc` option, including `"none"` for the unadjusted
means, with a line for every state. Each treated state has its own colour and the
never-treated states are grey.

```@example plots
figure = plot(plot_data)
save("plotting_conditional.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends, one panel per ccc option](plotting_conditional.png)

The `ccc` and `treated_states` arguments select the options and the treated states to show.
A treated state keeps its colour whichever states are shown.

```@example plots
figure = plot(plot_data; treated_states = 72)
save("plotting_conditional_one.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends of treated state 72](plotting_conditional_one.png)

Three treated states under two of the options:

```@example plots
figure = plot(plot_data; ccc = ["none", "state"],
              treated_states = [71, 58, 64])
save("plotting_conditional_some.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends of three treated states under two ccc options](plotting_conditional_some.png)

### Outcome means, one panel per treated state

With `facet = "state"` there is one panel per treated state and `ccc` option, a row of
panels per treated state; `ccc` selects the options, here a single one. Each panel compares
the treated state with the average of its control states. Parallel trends under that option
means that the two lines move in parallel before the dashed line.

```@example plots
figure = plot(plot_data; facet = "state", ccc = "state")
save("plotting_conditional_state.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends, one panel per treated state](plotting_conditional_state.png)

The average weights the control states by their numbers of observations; `weights = false`
weights them equally. `controls = "separate"` draws a line for each control state instead
of their average:

```@example plots
figure = plot(plot_data; facet = "state", ccc = "state", controls = "separate",
              treated_states = [71, 58])
save("plotting_conditional_state_separate.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends of two treated states with each control state](plotting_conditional_state_separate.png)

`ncols` sets the number of columns of panels, here for three selected states:

```@example plots
figure = plot(plot_data; facet = "state", ccc = "hom", weights = false,
              treated_states = [71, 64, 72], ncols = 1, size = (640, 780))
save("plotting_conditional_state_column.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends of three treated states in one column](plotting_conditional_state_column.png)

With several `ccc` options (or without `ccc`, which shows them all), each treated state has
a row of panels, one per option:

```@example plots
figure = plot(plot_data; facet = "state", ccc = ["hom", "state", "add"],
              treated_states = [71, 72])
save("plotting_conditional_state_ccc.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends of two treated states under three ccc options](plotting_conditional_state_ccc.png)

### Changes since the pre-treatment period

`didint()` compares the change in a treated state's adjusted mean outcome since its
pre-treatment period with the average change of its control states over the same periods.
`plot_type = "conditional2"` draws exactly those changes, one panel per treated state and
`ccc` option. Both lines are zero in the pre-treatment period. Here under a single option:

```@example plots
figure = plot(plot_data; plot_type = "conditional2", ccc = "state")
save("plotting_conditional2.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Conditional parallel trends as changes since the pre-treatment period](plotting_conditional2.png)

- In a treated period, the gap between the two lines is the ATT that `didint()` estimates
  for that state and period.
- Before treatment, the gap is a placebo estimate: the same comparison where no effect is
  expected yet. Parallel trends requires these gaps to be close to zero. Gaps that drift
  steadily towards or away from zero as treatment approaches suggest that the treated
  state was on a different trend, and a jump just before treatment suggests anticipation.

The gaps before and after treatment are in the same units, so the placebo gaps show how
large a violation of parallel trends would be next to the estimated effect.

Without `ccc`, every option is shown, with a row of panels per treated state and a column
per option. For one treated state this compares the options side by side:

```@example plots
figure = plot(plot_data; plot_type = "conditional2", treated_states = 72)
save("plotting_conditional2_ccc.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Changes since the pre-treatment period of treated state 72 under each ccc option](plotting_conditional2_ccc.png)

`controls = "separate"` draws each control state, here for one treated state:

```@example plots
figure = plot(plot_data; plot_type = "conditional2", ccc = "state",
              controls = "separate", treated_states = 72)
save("plotting_conditional2_separate.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Changes since the pre-treatment period of treated state 72 and of each of its control states](plotting_conditional2_separate.png)

### Not-yet-treated states as controls

The control states of a treated state are the never-treated states. Calling
`didint_plot()` with `notyet = true` adds the states of later cohorts while they are not
yet treated, as `didint()` does with the same argument. A not-yet-treated state is a
control up to the period before its own treatment.

```@example plots
plot_notyet = didint_plot(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    covariates = ["asian", "black", "male"],
    notyet     = true
);
nothing # hide
```

```@example plots
figure = plot(plot_notyet; plot_type = "conditional2", ccc = "state",
              controls = "separate", treated_states = [71, 64])
save("plotting_conditional_notyet.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Changes since the pre-treatment period with not-yet-treated control states](plotting_conditional_notyet.png)

---

## Integrated Parallel Trends Plots

The integrated parallel trends plots are those of `ccc = "int"`. There is one panel per
treated state, in which every state's outcome is predicted at the treated state's
covariate means.

### Changes since the previous period

`plot_type = "integrated"` draws the plot of Karim & Webb (2025): the change in the
treated state's predicted outcome since the previous period (the `delta` column of
`plot_data.integrated`) next to the average change of its control states. If the two lines
stay close to each other in the periods before the dashed line, integrated parallel trends
is plausible for that state. In the period at the dashed line, the gap between them is the
ATT that `didint()` estimates for that state's first treated period.

```@example plots
figure = plot(plot_data; plot_type = "integrated")
save("plotting_integrated.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Integrated parallel trends, changes since the previous period](plotting_integrated.png)

As for the conditional plots, `weights = false` weights the control states equally,
`controls = "separate"` draws each control state, and `notyet = true` in `didint_plot()`
adds the not-yet-treated states to the controls:

```@example plots
figure = plot(plot_data; plot_type = "integrated", controls = "separate",
              treated_states = [71, 58])
save("plotting_integrated_separate.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Integrated parallel trends of two treated states with each control state](plotting_integrated_separate.png)

### Changes since the pre-treatment period

`plot_type = "integrated2"` draws the same comparison with every change measured from the
pre-treatment period of the treated state, as `didint()` measures it. Both lines are zero
in that period.

```@example plots
figure = plot(plot_data; plot_type = "integrated2")
save("plotting_integrated2.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![Integrated parallel trends, changes since the pre-treatment period](plotting_integrated2.png)

As in `"conditional2"`, the gap between the two lines is the ATT that `didint()`
estimates in a treated period, and a placebo estimate before treatment. A trend that
differs between the treated state and its controls accumulates in this plot, so it is
easier to see than in the changes from one period to the next; on the other hand every
point of a line shares the noise of the pre-treatment period. The data behind this plot
are in `plot_data.integrated2`:

```@example plots
first(plot_data.integrated2, 4)
```

---

## Saving and Customising a Plot

`plot` returns a Makie `Figure`, which can be saved with `save` or modified like any other
Makie figure. The `title`, `size` and `ncols` arguments change its appearance
(`title = ""` leaves the title out), and `theme` takes a Makie `Theme` that replaces the
default look.

```@example plots
figure = plot(plot_data; plot_type = "integrated2", treated_states = [72, 57],
              title = "Merit aid and college attendance", size = (760, 420))
save("plotting_custom.png", figure; px_per_unit = 1.5) # hide
nothing # hide
```

![An integrated parallel trends plot with a title and size of its own](plotting_custom.png)

```julia
save("trends.png", figure)
```
