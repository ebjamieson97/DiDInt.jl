# DiDInt.jl

[![CI](https://github.com/ebjamieson97/DiDInt.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ebjamieson97/DiDInt.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ebjamieson97/DiDInt.jl/graph/badge.svg?token=7Y95RD0P9N)](https://codecov.io/gh/ebjamieson97/DiDInt.jl)
[![Docs Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ebjamieson97.github.io/DiDInt.jl/stable)

Integrated difference-in-differences!

**DiDInt.jl** introduces two functions:
- `didint()`, which is used for estimating the average effect of treatment on the treated (ATT) while accounting for various violations of the common causal covariates (CCC) assumption
- `didint_plot()` which produces datasets that can easily be used to make integrated parallel trends plots, as well as conditional parallel trends plots that help determine which CCC violation should be accounted for when using `didint()`. `didint_plot` can also produce datasets that can easily be used to make event study plots.

To learn more about the common causal covariates assumption and integrated difference-in-differences please see [Karim & Webb (2025)](https://arxiv.org/abs/2412.14447).

A Stata wrapper for this package is available here: **[didintjl](https://github.com/ebjamieson97/didintjl)**.

An R wrapper for this package is available here: **[didintrjl](https://github.com/ebjamieson97/didintrjl)**.

See more detailed documentation on the [documentation page](https://ebjamieson97.github.io/DiDInt.jl/). Changes between versions are listed in the [changelog](CHANGELOG.md).

## Installation 
```
using Pkg
Pkg.add("DiDInt")
```

## `didint()` - Estimate ATT

The function `didint()` estimates the ATT while adjusting for covariates that may vary by state, time, or both. It works for both staggered and common adoption, and reports heteroskedasticity-robust and jackknife standard errors along with randomization inference p-values.

```julia
using DiDInt
using CSV, DataFrames # To read in data

# Define treated states and their respective treatment times
treated_states = ["34", "57", "58", "59", "61", "64", "71", "72", "85", "88"]
treated_times = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000]

# Load data (test/data/merit.csv in this repository)
test_data = CSV.read("merit.csv", DataFrame)
test_data.state = string.(test_data.state)

# Run DID-INT model
result = didint("coll", "state", "year", test_data,
                treatment_times = treated_times,
                treated_states = treated_states,
                seed = 1234, ccc = "state",
                covariates = [:male, :asian, :black],
                agg = "cohort", nperm = 399)

# DID-INT results
#   Outcome:          coll
#   Model:            Region-Varying DID-INT (ccc = "state")
#   Covariates:       male, asian, black
#   Aggregation:      "cohort" (by treatment time)
#   Weighting:        "both"
#   Controls:         never-treated states
#   Periods:          1 year, from 1989 to 2000
#   Standard errors:  HC1 (se) and jackknife (jknife_se)
#   Permutations:     399 (seed 1234)
#
# Aggregate ATT
#  estimate       att        se          pval         jknife_se  jknife_pval  ri_pval
# ────────────────────────────────────────────────────────────────────────────────────
#  Aggregate ATT  0.0462414  0.00769881  0.000959237  0.0123516   0.00046902   0.0975
#
# Sub-aggregate ATTs by treatment time (7)
#  treatment_time  att          se          pval       jknife_se        jknife_pval      ri_pval  weight
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────
#  1991-01-01       0.0553354   0.0222024   0.0342914  missing          missing           0.4575  0.201796
#  1993-01-01       0.0504489   0.0161905   0.0169391  missing          missing           0.5175  0.191535
#  1996-01-01       0.0424116   0.0255707   0.172536   missing          missing           0.5725  0.0756734
#  1997-01-01       0.0606879   0.0236147   0.037012         0.0328756        0.0708186   0.28    0.321077
#  1998-01-01       0.0398137   0.0406649   0.372511         0.0813398        0.626647    0.5825  0.108593
#  1999-01-01       0.00173622  0.00350574  0.707256   missing          missing           0.985   0.0354853
#  2000-01-01      -0.0253317   0.0333203   0.586179         0.0971445        0.795346    0.755   0.0658401

# The estimates are stored as tables...
result.agg          # the aggregate ATT with its standard errors and p-values
result.sub          # the sub-aggregate ATTs
result.att_st       # the ATT of each treated state in each period

# ...are available as vectors...
coef(result)                    # [0.0462414]
stderror(result)                # [0.00769881]
coef(result; level = :sub)      # the seven cohort ATTs

# ...and can be collected in one wide table
DataFrame(result)
```

`didint()` returns a `DiDIntResult`, which prints the summary shown above.

Other aggregations are available through the `agg` argument (`"cohort"`, `"simple"`, `"state"`, `"sgt"`, `"time"`, `"calendar"`), as are separate ATTs by policy when treated states adopt different policies (`"policy"`, `"policy_gt"`, `"policy_cohort"`, `"policy_calendar"`, `"policy_time"`). See the [documentation page](https://ebjamieson97.github.io/DiDInt.jl/) for details.

## `didint_plot()` - Prepare Data for Visualization

The `didint_plot()` function produces a dataset in a long format that can easily be used for plotting integrated parallel trends (the default), conditional parallel trends (`integrated = false`) or event study plots (`event = true`).

```julia

# Generate data for an integrated parallel trends plot: for each treated state, the change
# since the previous period in the predicted outcome of that state and of its comparison
# states, all evaluated at the treated state's covariate means
plot_integrated = didint_plot("coll", "state", "year", test_data,
                              treatment_times = treated_times,
                              treated_states = treated_states,
                              covariates = ["asian", "black", "male"])
plot_integrated[1:4, :]
# 4×9 DataFrame
#  Row │ treated_state  state   time    delta       ccc     period  treat_period  start_date  period_length
#      │ String         String  String  Float64     String  Int64   Int64         String      String
# ─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ 34             34      1990    -0.0791523  int          1            11  1989        1 year
#    2 │ 34             11      1990    -0.125941   int          1            11  1989        1 year
#    3 │ 34             12      1990     0.0223095  int          1            11  1989        1 year
#    4 │ 34             13      1990    -0.160495   int          1            11  1989        1 year


# Generate data for a conditional parallel trends plot: the covariate-adjusted mean
# outcome of every state and period, for each CCC specification
plot_parallel = didint_plot("coll", "state", "year", test_data,
                            treatment_times = treated_times,
                            treated_states = treated_states,
                            integrated = false, covariates = ["asian", "black", "male"])
plot_parallel[1:4, :]
# 4×8 DataFrame
#  Row │ state    time     lambda    ccc      period  start_date  treat_period  period_length
#      │ String?  String?  Float64?  String?  Int64?  String?     Int64?        String
# ─────┼──────────────────────────────────────────────────────────────────────────────────────
#    1 │ 11       1989     0.467895  hom           0  1989             missing  1 year
#    2 │ 11       1990     0.316912  hom           1  1989             missing  1 year
#    3 │ 11       1991     0.445936  hom           2  1989             missing  1 year
#    4 │ 11       1992     0.579997  hom           3  1989             missing  1 year


# Generate data for an event study plot
plot_event = didint_plot("coll", "state", "year", test_data,
                               treatment_times = treated_times, 
                               treated_states = treated_states,
                               event = true, covariates = ["asian", "black", "male"])
plot_event[1:4, :]
# 4×8 DataFrame
#  Row │ ccc     time_since_treatment  y         se         ci_lower   ci_upper  ngroup  period_length
#      │ String  Int64                 Float64   Float64?   Float64?   Float64?  Int64   String
# ─────┼───────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ hom                      -11  0.528292  0.0964386  -0.697077  1.75366        2  1 year
#    2 │ hom                      -10  0.447675  0.0554866   0.208936  0.686415       3  1 year
#    3 │ hom                       -9  0.470806  0.0223596   0.408725  0.532886       5  1 year
#    4 │ hom                       -8  0.452972  0.0248907   0.392067  0.513878       7  1 year
```

## Citations
- Karim & Webb (2025). ["Good Controls Gone Bad: Difference-in-Differences with Covariates".](https://arxiv.org/abs/2412.14447)
- MacKinnon & Webb (2020). ["Randomization inference for difference-in-differences with few treated clusters".](https://doi.org/10.1016/j.jeconom.2020.04.024)
