# Examples

The examples below use the `merit.csv` dataset included in the package's test suite
(`test/data/merit.csv`). It contains individual-level survey data on college attendance
(`coll`), state identifiers (`state`), survey year (`year`), and demographic covariates
(`male`, `asian`, `black`). The dataset covers U.S. states over the period 1989–2000 and
is used to study the effect of merit scholarships on college attendance.

```@example didint
using DiDInt, CSV, DataFrames

merit = CSV.read(joinpath(@__DIR__, "../../test/data/merit.csv"), DataFrame);
first(merit, 5)
```

---

## Staggered Adoption

This example illustrates a staggered adoption setting in which ten states adopted
merit-aid programmes at different points between 1991 and 2000.

### Running `didint()`

```@example didint
# Note that treated_states & state column in the merit data could alternatively
# be string vectors, here they just happen to be numeric codes instead of state names
treated_states = [34, 57, 58, 59, 61, 64, 71, 72, 85, 88];
treated_times  = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000];

result = DiDInt.didint(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    seed       = 1234,
    ccc        = "state",
    covariates = [:male, :asian, :black],
    agg        = "cohort",
    nperm      = 399
);
```

Inspect the aggregate-level ATT and inference columns:

```@example didint
result[:, 3:9]
```

And the cohort-level ATT and inference columns:

```@example didint
result[:, vcat([1, 2], 10:ncol(result))]
```

### Running `didint_plot()`

`didint_plot()` can produce data for two types of plots depending on the `event` argument.

By default, the outputted dataset will include the parallel trends (or event study) data for each `ccc` option (`"hom"`,
`"time"`, `"state"`, `"add"`, `"int"`) as well as the data estimated without covariates (`"none"`). You can specify which
`ccc` options you want to plot by inputting a vector of strings to the `ccc` argument.

#### Parallel trends plot data

Setting `event = false` returns a long-format DataFrame of residualised outcome means
by state and period, suitable for plotting parallel trends across CCC specifications.

```@example didint
plot_parallel = DiDInt.didint_plot(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    event      = false,
    covariates = ["asian", "black", "male"]
);

plot_parallel[1:4, :]
```

#### Event study plot data

Setting `event = true` returns a DataFrame of treated-state means by period relative
to treatment, along with standard errors and confidence interval bounds, suitable for
plotting an event study.

```@example didint
plot_event = DiDInt.didint_plot(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    event      = true,
    covariates = ["asian", "black", "male"]
);

plot_event[1:4, :]
```

The `lambda` column in the parallel trends output contains the residualised outcome
for each state-period cell. The `ccc` column indicates which model was used (`"hom"`,
`"time"`, `"state"`, `"add"`, or `"int"`). In the event study output,
`time_since_treatment = 0` marks the treatment period; negative values are
pre-treatment and positive values are post-treatment.

---

## Common Adoption

This example illustrates a **common adoption** setting: a single treated state
(state `"71"`, treated in 1991) compared against a single untreated control state
(state `"73"`). The data are filtered to just these two states before calling `didint()`.

```@example didint
merit_common = filter(row -> row.state ∈ [71, 73], merit);

result_common = DiDInt.didint(
    "coll", "state", "year", merit_common;
    treated_states  = [71],
    treatment_times = [1991],
    seed       = 1234,
    ccc        = "hom",
    covariates = [:male, :asian, :black],
    agg        = "cohort"
);
```

```@example didint
result_common
```


### Trying a different CCC assumption

The `ccc` argument controls the assumption made about how covariate effects vary across
states and time. Swapping `ccc = "hom"` for any of `"time"`, `"state"`, `"add"`, or
`"int"` changes whether covariate effects are allowed to vary by time, by state,
additively by both, or interactively by both, respectively. For example:

```@example didint
result_common_int = DiDInt.didint(
    "coll", "state", "year", merit_common;
    treated_states  = [71],
    treatment_times = [1991],
    seed       = 1234,
    ccc        = "int",
    covariates = [:male, :asian, :black],
    agg        = "cohort"
);
```

```@example didint
result_common_int
```

See the [Functions](@ref) for function syntax and the [Under the Hood](@ref) page for more details.