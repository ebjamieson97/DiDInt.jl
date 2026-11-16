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

result = didint(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    seed       = 1234,
    ccc        = "state",
    covariates = [:male, :asian, :black],
    agg        = "cohort",
    nperm      = 399
)
```

`didint()` returns a `DiDIntResult`. Displaying it, as above, prints a summary of the
specification, the aggregate ATT, and the sub-aggregate ATTs (here, one per treatment
time). Jackknife standard errors are `missing` for cohorts with a single treated state
(see the [Jackknife](@ref) section).

### Working with the results

The estimates are stored as tables in the fields of the result. The aggregate ATT with
its standard errors and p-values:

```@example didint
result.agg
```

The sub-aggregate ATTs, with the weight of each in the aggregate ATT:

```@example didint
result.sub
```

The ATT for each treated state in each period from its treatment time on, which is what
every other estimate is computed from:

```@example didint
first(result.att_st, 5)
```

The estimates and standard errors are also available as vectors:

```@example didint
coef(result), stderror(result), stderror(result; type = :jackknife)
```

```@example didint
coefnames(result; level = :sub), coef(result; level = :sub)
```

Finally, every result can be collected in a single wide table with `DataFrame(result)`
(the same table is passed to the Stata and R wrappers), and a result can be written
straight to a file, for example with `CSV.write("results.csv", result)`:

```@example didint
names(DataFrame(result))
```

### Other aggregations

The `agg` argument controls how the ATTs for each treated state and period are grouped
into sub-aggregate ATTs (see the [Aggregation](@ref) section). For example, `agg = "time"`
groups them by the number of periods since treatment:

```@example didint
result_time = didint(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    seed       = 1234,
    ccc        = "state",
    covariates = [:male, :asian, :black],
    agg        = "time",
    nperm      = 399
);

result_time.sub
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
plot_parallel = didint_plot(
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
plot_event = didint_plot(
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
`"time"`, `"state"`, `"add"`, `"int"`, or `"none"` for no covariates). In the event study output,
`time_since_treatment = 0` marks the treatment period; negative values are
pre-treatment and positive values are post-treatment.

---

## Separate ATTs by Policy

When treated states adopt different policies, the policy aggregation options report an
ATT for each policy. The merit data does not record different types of merit-aid
programmes, so purely for illustration the ten treated states are split here into two
made-up policy types, `"A"` and `"B"`, using the `treated_policies` argument (the `pvar`
argument can be used instead when the policy is stored in a column of the data).

```@example didint
treated_policies = ["A", "A", "B", "A", "B", "A", "B", "B", "A", "B"];

result_policy = didint(
    "coll", "state", "year", merit;
    treatment_times  = treated_times,
    treated_states   = treated_states,
    treated_policies = treated_policies,
    seed       = 1234,
    ccc        = "state",
    covariates = [:male, :asian, :black],
    agg        = "policy",
    policy_aggregation = true,
    nperm      = 399
)
```

Setting `policy_aggregation = true` also reports the pooled ATT across policies, here
weighting the two policies equally (the `policy_weight` column).
See the [Policies](@ref) section for the other policy aggregation options and for the
weighting of policies in the pooled ATT.

---

## Common Adoption

When every treated state is treated at the same time, `didint()` is used in exactly the
same way, and every aggregation option is available. This example keeps the two states
that adopted merit-aid programmes in 1998 along with the never-treated states. A single
treatment time applies to every treated state.

```@example didint
merit_common = filter(row -> row.state ∈ [57, 72] || !(row.state ∈ treated_states), merit);

result_common = didint(
    "coll", "state", "year", merit_common;
    treated_states  = [57, 72],
    treatment_times = 1998,
    seed       = 1234,
    ccc        = "hom",
    covariates = [:male, :asian, :black],
    agg        = "time",
    nperm      = 399
)
```

### Trying a different CCC assumption

The `ccc` argument controls the assumption made about how covariate effects vary across
states and time. Swapping `ccc = "hom"` for any of `"time"`, `"state"`, `"add"`, or
`"int"` changes whether covariate effects are allowed to vary by time, by state,
additively by both, or interactively by both, respectively. For example:

```@example didint
result_common_int = didint(
    "coll", "state", "year", merit_common;
    treated_states  = [57, 72],
    treatment_times = 1998,
    seed       = 1234,
    ccc        = "int",
    covariates = [:male, :asian, :black],
    agg        = "time",
    nperm      = 399
);

result_common_int.sub
```

---

## A Single Control State

This example compares a single treated state (state `71`, treated in 1991) against a
single untreated control state (state `73`). With `agg = "sgt"`, an ATT is reported for
state `71` in each period from 1991 onwards. With only one control state, the usual
standard error of each of these ATTs cannot be computed; setting `edgecase = true`
recovers it from the first-stage estimates instead.

```@example didint
merit_pair = filter(row -> row.state ∈ [71, 73], merit);

result_pair = didint(
    "coll", "state", "year", merit_pair;
    treated_states  = [71],
    treatment_times = [1991],
    seed       = 1234,
    ccc        = "hom",
    covariates = [:male, :asian, :black],
    agg        = "sgt",
    edgecase   = true
);

select(result_pair.sub, :state, :t, :att, :se)
```

No p-value accompanies these standard errors, as the regression behind each ATT has no
residual degrees of freedom. With one treated state and one control state there is also nothing for the jackknife to
compute and only one alternative treatment assignment for the randomization inference, so
`didint()` issues warnings to that effect. For more information on these standard errors,
including a caveat for panel data, see the section
[Computation of Edge Case Standard Errors](@ref) on the "Details" page of the documentation site.

See the [Functions](@ref) for function syntax and the [Under the Hood](@ref) page for more details.
