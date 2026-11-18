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

`didint_plot()` prepares the data for parallel trends plots and event study plots, and
returns them as a `DiDIntPlot`. One call prepares the data for every kind of plot:

```@example didint
plot_data = didint_plot(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    covariates = ["asian", "black", "male"]
)
```

With a [Makie](https://docs.makie.org) backend loaded, `plot(plot_data; plot_type = ...)`
draws each of these: the [Plotting Examples](@ref) page shows every plot, and the
[Plotting with `didint_plot()`](@ref) section gives the formulas behind them. The data of
each kind of plot are in the table of the same name. The `plot_type` argument of
`didint_plot()` selects the kinds to prepare, and `ccc` the options.

#### Conditional parallel trends data

The `conditional` table has the mean outcome of each state and period, residualised by
the covariates as each `ccc` option other than `"int"` specifies (`lambda`), and the number
of observations behind it (`n`). `"none"` is the option without covariates.

```@example didint
plot_data.conditional[1:4, :]
```

How the series of a treated state changes over time relative to the series of its
controls is exactly the comparison that `didint()` bases its estimates on under these
options. The *level* of a series rests on a convention (the adjusted means are evaluated
at covariates equal to zero), but another choice of covariate values would shift the
series together and leave their parallelism unaffected.

The `conditional2` table has those changes themselves: for each treated state
(`treated_state`), the change in `lambda` of that state and of each of its control states
(`state`) since the period before its treatment (`base_period`).

```@example didint
plot_data.conditional2[1:4, :]
```

#### Integrated parallel trends data

Under `ccc = "int"` another choice of covariate values would change each point of a
series by its own amount, so the levels cannot be compared. The `integrated` table
therefore has, for each treated state (`treated_state_covars`), the change in the
predicted outcome since the previous period (`delta`) of that state and of each of its
control states (`state`), with every state evaluated at the covariate means of
`treated_state_covars`. `treat_period` is the period in which that treated state is first
treated, and `n` the number of observations of `state` in the two periods.

```@example didint
plot_data.integrated[1:4, :]
```

The `integrated2` table has the same changes measured from the period before the treated
state's treatment (`base_period`). In both tables the control states are the never-treated
states; `notyet = true` adds the states of later cohorts while they are not yet treated.

#### Event study data

The `event` table has the treated-state means by period relative to treatment, with
standard errors and confidence interval bounds.

```@example didint
plot_data.event[1:4, :]
```

The `ccc` column indicates which model was used (`"hom"`, `"time"`, `"state"`, `"add"`,
`"int"`, or `"none"` for no covariates), and `time_since_treatment = 0` marks the
treatment period; negative values are pre-treatment and positive values are
post-treatment.

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

### Cohort ATTs With and Without Policies

`agg = "cohort"` reports one ATT per treatment time and ignores the policies, whereas
`agg = "policy_cohort"` reports one ATT per treatment time *within each policy*. With the
made-up policies above, the two treated states of the 1998 cohort have different policies,
as do the two treated states of the 2000 cohort; the two states of the 1997 cohort both
have policy `"A"`, and every other cohort consists of a single state.

Under `agg = "cohort"` there are seven cohort ATTs, whether or not policies are given:

```@example didint
result_cohort = didint(
    "coll", "state", "year", merit;
    treatment_times  = treated_times,
    treated_states   = treated_states,
    treated_policies = treated_policies,
    seed       = 1234,
    ccc        = "state",
    covariates = [:male, :asian, :black],
    agg        = "cohort",
    nperm      = 399
);

select(result_cohort.sub, :treatment_time, :att, :se, :jknife_se, :weight)
```

Under `agg = "policy_cohort"` the cohorts are formed within each policy:

```@example didint
result_policy_cohort = didint(
    "coll", "state", "year", merit;
    treatment_times  = treated_times,
    treated_states   = treated_states,
    treated_policies = treated_policies,
    seed       = 1234,
    ccc        = "state",
    covariates = [:male, :asian, :black],
    agg        = "policy_cohort",
    nperm      = 399
)
```

Comparing the two sets of results:

- The 1998 and 2000 cohorts, which contain both policies, are each split in two: one ATT
  for the state with policy `"A"` and one for the state with policy `"B"`. Each of these
  is now based on a single treated state, so it has no jackknife standard error. The
  split can matter: the 2000 cohort has an ATT of about -0.025 as a whole, which is made
  up of about -0.046 under policy `"A"` and 0.070 under policy `"B"`.
- The cohorts that contain a single policy (1991, 1993, 1996, 1997 and 1999) appear once,
  under that policy, with the same ATT and standard errors as under `agg = "cohort"`.
  Their randomization inference p-values differ slightly, because for the policy
  aggregation options a state's policy is reassigned together with its treatment time.
- The `weight` column now gives the weight of each cohort ATT *within its policy*, and the
  ATT of each policy (the first table of the output) is the weighted mean of the cohort
  ATTs of that policy. There is no single aggregate ATT; setting
  `policy_aggregation = true` adds the pooled ATT across policies, as in the example above.

---

## The Aggregate ATT Across Aggregation Options

With the default weighting (`weighting = "both"`), and also with `weighting = "att"`, the
aggregate ATT does not depend on the `agg` option (see the
[Aggregate ATT Equivalency](@ref) section). Running the staggered adoption example above
with every (non-policy) `agg` option and the default weighting:

```@example didint
aggs = ["cohort", "simple", "state", "sgt", "time", "calendar"]

results = [didint(
    "coll", "state", "year", merit;
    treatment_times = treated_times,
    treated_states  = treated_states,
    seed       = 1234,
    ccc        = "state",
    covariates = [:male, :asian, :black],
    agg        = agg,
    nperm      = 399
) for agg in aggs]

# The aggregate ATT and its inference, one row per agg option
aggregate(result) = result.agg[1, :]
DataFrame(
    agg         = aggs,
    agg_att     = [aggregate(result).att for result in results],
    se          = [aggregate(result).se for result in results],
    pval        = [aggregate(result).pval for result in results],
    jknife_se   = [aggregate(result).jknife_se for result in results],
    jknife_pval = [aggregate(result).jknife_pval for result in results],
    ri_pval     = [aggregate(result).ri_pval for result in results]
)
```

The aggregate ATT (`agg_att`) is identical in every row. Its HC standard error and p-value
(`se`, `pval`) differ, as they are computed from the spread of the sub-aggregate ATTs of
each option. The jackknife standard error and p-value do not differ: each removal of a state
recomputes the same aggregate ATT, whatever the grouping. Neither does the randomization
inference p-value, except under `"state"` and `"sgt"`, whose randomization inference draws
an additional stand-in state for each treated state (see the
[Randomization Inference](@ref) section).

---

## Common Adoption

When every treated state is treated at the same time, `didint()` is used in exactly the
same way, and every aggregation option is available. This example keeps the two states
that adopted merit-aid programmes in 1998 along with the never-treated states. A single
treatment time applies to every treated state.

```@example didint
merit_common = filter(
    row -> row.state ∈ [57, 72] || !(row.state ∈ treated_states),
    merit
);

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
