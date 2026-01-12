# DiDInt.jl

[![CI](https://github.com/ebjamieson97/DiDInt.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ebjamieson97/DiDInt.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ebjamieson97/DiDInt.jl/graph/badge.svg?token=7Y95RD0P9N)](https://codecov.io/gh/ebjamieson97/DiDInt.jl)
[![Docs Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ebjamieson97.github.io/DiDInt.jl/dev)

Intersection difference-in-differences!

**DiDInt.jl** introduces two functions:
- `didint()`, which is used for estimating the average effect of treatment on the treated (ATT) while accounting for various violations of the common causal covariates (CCC) assumption
- `didint_plot()` which produces datasets that can easily be used to make parallel trends plots that help determine which CCC violation should be accounted for when using `didint()`. `didint_plot` can also produce datasets that can easily be used to make event study plots.

To learn more about the common causal covariates assumption and intersection difference-in-differences please see [Karim & Webb (2025)](https://arxiv.org/abs/2412.14447).

A Stata wrapper for this package is available here: **[didintjl](https://github.com/ebjamieson97/didintjl)**.

See more detailed documentation on the [documentation page](https://ebjamieson97.github.io/DiDInt.jl/).

## Installation 
```
using Pkg
Pkg.add(url="https://github.com/ebjamieson97/DiDInt.jl")
```

## `didint()` - Estimate ATT

The function `didint()` estimates the ATT while adjusting for covariates that may vary by state, time, or both.

```julia
using DiDInt
using CSV, DataFrames # To read in data

# Define treated states and their respective treatment times
treated_states = ["34", "57", "58", "59", "61", "64", "71", "72", "85", "88"]
treated_times = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000]

# Load data
test_data = CSV.read("merit.csv", DataFrame)

# Run DID-INT model
result = DiDInt.didint("coll", "state", "year", test_data,
                       treatment_times = treated_times,
                       treated_states = treated_states, 
                       seed = 1234, ccc = "state",
                       covariates = [:male, :asian, :black],
                       agg = "cohort", nperm = 399)

# Show aggregate level results
result[:, 3:9]
        
# 7×7 DataFrame
#  Row │ agg_att          se_agg_att       pval_agg_att      jknifese_agg_att  jknifepval_agg_att  ri_pval_agg_att  nperm     
#      │ Float64?         Float64?         Float64?          Float64?          Float64?            Float64?         Float64?
# ─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │       0.0545816        0.0135013        0.00678187         0.0172599          0.00512895        0.0401003      399.0
#    2 │ missing          missing          missing            missing            missing           missing          missing
#    3 │ missing          missing          missing            missing            missing           missing          missing
#    4 │ missing          missing          missing            missing            missing           missing          missing
#    5 │ missing          missing          missing            missing            missing           missing          missing
#    6 │ missing          missing          missing            missing            missing           missing          missing
#    7 │ missing          missing          missing            missing            missing           missing          missing

# Show sub-aggregate level results
result[:, vcat([1,2], 10:end)]

# 7×11 DataFrame
#  Row │ treatment_time  att_cohort  se_att_cohort  pval_att_cohort  jknifese_att_cohort  jknifepval_att_cohort  ri_pval_att_cohort  weights    period  start_date  end_date 
#      │ Dates.Date      Float64?    Float64?       Float64?         Float64?             Float64?               Float64?            Float64    String  String      String
# ─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ 1991-01-01       0.0740159      0.0247873      0.00349694       missing                 missing                   0.255639  0.201796   1 year  1989        2000
#    2 │ 1993-01-01       0.0845589      0.0210913      0.000129329      missing                 missing                   0.215539  0.191535   1 year  1989        2000
#    3 │ 1996-01-01       0.0333588      0.0308039      0.283737         missing                 missing                   0.684211  0.0756734  1 year  1989        2000
#    4 │ 1997-01-01       0.0654413      0.0289131      0.028375               0.0389711               0.109479            0.288221  0.321077   1 year  1989        2000
#    5 │ 1998-01-01       0.0444878      0.0490408      0.370708               0.0900249               0.62685             0.528822  0.108593   1 year  1989        2000
#    6 │ 1999-01-01      -0.0211062      0.0207727      0.321739         missing                 missing                   0.842105  0.0354853  1 year  1989        2000
#    7 │ 2000-01-01      -0.0633142      0.0997027      0.539666               0.0971782               0.52251             0.401003  0.0658401  1 year  1989        2000

```

## `didint_plot()` - Prepare Data for Visualization

The `didint_plot()` function produces a dataset in a long format that can easily be used for plotting parallel trends or event study plots.

```julia

# Generate data for a parallel trends plot
plot_parallel = DiDInt.didint_plot("coll", "state", "year", test_data,
                                   treatment_times = treated_times, 
                                   treated_states = treated_states,
                                   event = false, covariates = ["asian", "black", "male"])
plot_parallel[1:4, :]
# 4×8 DataFrame
#  Row │ state    time     lambda    ccc      period  start_date  treat_period  period_length 
#      │ String?  String?  Float64?  String?  Int64?  String?     Int64?        String
# ─────┼──────────────────────────────────────────────────────────────────────────────────────
#    1 │ 12       1989     0.503763  hom           0  1989             missing  1 year
#    2 │ 12       1990     0.483555  hom           1  1989             missing  1 year
#    3 │ 12       1991     0.539888  hom           2  1989             missing  1 year
#    4 │ 12       1992     0.226398  hom           3  1989             missing  1 year


# Generate data for an event study plot
plot_event = DiDInt.didint_plot("coll", "state", "year", test_data,
                               treatment_times = treated_times, 
                               treated_states = treated_states,
                               event = true, covariates = ["asian", "black", "male"])
plot_event[1:4, :]
# 4×8 DataFrame
#  Row │ ccc     time_since_treatment  y         se         ci_lower    ci_upper  ngroup  period_length 
#      │ String  Int64                 Float64   Float64?   Float64?    Float64?  Int64   String
# ─────┼────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ hom                      -11  0.528292  0.243      -2.55931    3.6159         2  1 year
#    2 │ hom                      -10  0.447675  0.116099   -0.0518564  0.947207       3  1 year
#    3 │ hom                       -9  0.470806  0.0278348   0.393524   0.548087       5  1 year
#    4 │ hom                       -8  0.452972  0.030141    0.37922    0.526725       7  1 year
```

## Citations
- Karim & Webb (2025). ["Good Controls Gone Bad: Difference-in-Differences with Covariates".](https://arxiv.org/abs/2412.14447)
- MacKinnon & Webb (2020). ["Randomization inference for difference-in-differences with few treated clusters".](https://doi.org/10.1016/j.jeconom.2020.04.024)