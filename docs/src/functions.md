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

## `didint_plot()` - Prepare Data for Visualization

The `didint_plot()` function produces a dataset in a long format that can easily be used for plotting parallel trends or event study plots.

```@docs
didint_plot
```
