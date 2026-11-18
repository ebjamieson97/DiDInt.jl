# Functions

## `didint()` - Estimate ATT

The function `didint()` estimates the ATT while adjusting for covariates that may vary by state, time, or both.

```@docs
didint
```

## `DiDIntResult` - Results of `didint()`

`didint()` returns its results as a `DiDIntResult`, which prints a summary when displayed and stores the estimates as tables.

```@docs
DiDIntResult
```

The estimates, their standard errors and their labels can also be extracted as vectors:

```@docs
coef(::DiDIntResult)
stderror(::DiDIntResult)
coefnames(::DiDIntResult)
```

## `didint_plot()` - Parallel Trends and Event Study Plots

The `didint_plot()` function prepares the data for parallel trends plots and event study plots.

```@docs
didint_plot
```

## `DiDIntPlot` - Plot Data and Plots

`didint_plot()` returns its data as a `DiDIntPlot`, which prints a summary when displayed, stores the data of each kind of plot as a table, and is drawn with `plot` once a [Makie](https://docs.makie.org) backend such as CairoMakie is loaded.

```@docs
DiDIntPlot
```
