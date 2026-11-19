# Under the Hood

Most pertinent details of **DiDInt.jl** can be found in the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447), however, some details that might be of particular interest to users of the **DiDInt.jl** package are elaborated upon here.

## Five Step Estimation Procedure

The DID-INT estimator is computed in five steps:

1. The means of the outcome of interest (conditional on covariates) are computed for each period for each state (denoted as $\lambda$ in the DID-INT paper).
2. Within state differences between these $\lambda$ values (long differences) are computed for treated states and for their controls.
3. An ATT is computed for each treated state $s$ in each period $t$ from its treatment time on, denoted $\mathrm{ATT}(s,t)$: the long difference of the treated state minus the (weighted) mean of the long differences of its controls.
4. The (weighted) means of the $\mathrm{ATT}(s,t)$ values within groups give the subaggregate ATTs (more details on the choice of groups in the [Aggregation](@ref) section).
5. The (weighted) mean of the subaggregate ATTs gives the aggregate ATT.

!!! note "Pooling Together Different Policies"
    If the policies evaluated with a policy aggregation are relatively similar, an overall effect across them may be of interest. This optional sixth step is off by default and is turned on by setting `policy_aggregation = true`; `policy_weighting` then sets how the policies are weighted.

More explicitly:

Step 1. In principle, the following regression is estimated to obtain conditional outcome means $\hat{\lambda}_{s,t}$:

```math
Y_{i,s,t} = \sum_s \sum_t \lambda_{s,t} I(s,t) + f\!\left(X^k_{i,s,t}\right) + \epsilon_{i,s,t}
```

where $I(s,t)$ is a dummy variable equal to 1 if observation $i$ belongs to state $s$ in period $t$, and $f\!\left(X^k_{i,s,t}\right)$ is a function of covariates whose form depends on the common causal covariates (CCC) violation that is being accounted for (see [Common Causal Covariates and Model Specification](@ref)). The estimated $\hat{\lambda}_{s,t}$ values are the covariate-adjusted means of the outcome for state $s$ in period $t$.

!!! note "Implementation"
    In practice, this regression is not estimated directly on the full dataset as this would be computationally prohibitive for datasets with many states, periods, and covariates. Rather, the Frisch-Waugh-Lovell theorem is applied in the case of homogeneous DID-INT, region-varying DID-INT, time-varying DID-INT, and two one-way DID-INT to recover covariate coefficients which can then be used in conjunction with the outcome and covariate means per $(s,t)$ group to retrieve the $\hat{\lambda}_{s,t}$ values. For two-way DID-INT, the regression is run separately for each $(s,t)$ group. The resulting $\hat{\lambda}_{s,t}$ values are numerically equivalent to those from the combined regression shown above. More details on estimating $\hat{\lambda}_{s,t}$ can be found in the [Estimating Lambda](@ref) section.

Step 2. For each treated state $s$ belonging to treatment cohort $g$ (the states first treated in period $g$), long differences are computed between each period $t \geq g$ and the period immediately prior to treatment, $r = g - 1$:

```math
\widehat{\Delta}_{s,t} = \left(\hat{\lambda}_{s,t} - \hat{\lambda}_{s,r}\right) + \left(\hat{\gamma}_{s,t} - \hat{\gamma}_{s,r}\right)'\,\overline{X}_{s,t}
```

The analogous difference is computed for each control state $c$ of the treated state $s$:

```math
\widehat{\Delta}_{c \to s,t} = \left(\hat{\lambda}_{c,t} - \hat{\lambda}_{c,r}\right) + \left(\hat{\gamma}_{c,t} - \hat{\gamma}_{c,r}\right)'\,\overline{X}_{s,t}
```

Here $\hat{\gamma}_{s,t}$ are the covariate coefficients that apply to state $s$ in period $t$ under the chosen CCC violation, and $\overline{X}_{s,t}$ are the means of the covariates of the *treated* state $s$ in period $t$. These are the integrated long differences of the DID-INT paper (equations 46 and 47): the change in the covariate-adjusted outcome is evaluated over the covariate distribution of the treated state, which, because $f$ is linear in the covariates, only requires that state's covariate means. More on when the covariate term matters in the [Integrated Long Differences](@ref) section.

Step 3. For each treated state $s$ and each period $t$ from its treatment time on, the ATT is the long difference of $s$ minus the weighted mean of the long differences of its controls (equation 48 of the DID-INT paper):

```math
\widehat{\mathrm{ATT}}(s,t) = \widehat{\Delta}_{s,t} - \frac{\sum_{c \in C(s,t)} v_{c}\, \widehat{\Delta}_{c \to s,t}}{\sum_{c \in C(s,t)} v_{c}}
```

where $C(s,t)$ is the set of controls: the never-treated states and, optionally, the states not yet treated in period $t$ (see [Not Yet Treated Cells](@ref)). In this step, weights $v_c$ are used (by default) according to the number of observations associated with each control long difference (see [Weighting](@ref)).

!!! note "Implementation"
    Each $\widehat{\mathrm{ATT}}(s,t)$ is computed as the coefficient $\beta$ in the weighted regression of the long differences of state $s$ and of its controls on an intercept and a treatment indicator, $\widehat{\Delta} = \alpha + \beta\, d + \varepsilon$, where $d = 1$ for the treated state and 0 for the controls. This coefficient equals the expression above and the regression also yields its standard error.

Step 4. Each subaggregate ATT is a weighted mean of the $\widehat{\mathrm{ATT}}(s,t)$ values that belong to its group $j$:

```math
\widehat{\mathrm{ATT}}_{j} = \frac{\sum_{(s,t) \in j} u_{s,t}\, \widehat{\mathrm{ATT}}(s,t)}{\sum_{(s,t) \in j} u_{s,t}}
```

By default the weight $u_{s,t}$ is the number of observations of treated state $s$ in period $t$ (see [Weighting](@ref)).

!!! note "Implementation"
    Each subaggregate ATT is computed as the coefficient $\beta_j$ in the weighted regression of the $\widehat{\mathrm{ATT}}(s,t)$ values of group $j$ on a constant, $\widehat{\mathrm{ATT}}(s,t) = \beta_j + \varepsilon_{s,t}$, with weights $u_{s,t}$. This coefficient equals the weighted mean above, and the regression also yields the standard error of the subaggregate ATT: the heteroskedasticity-consistent standard error of $\beta_j$, which reflects how much the $\widehat{\mathrm{ATT}}(s,t)$ values of the group vary around their weighted mean. It does not use the standard errors of the individual $\widehat{\mathrm{ATT}}(s,t)$ values from Step 3. A subaggregate ATT that consists of a single $\widehat{\mathrm{ATT}}(s,t)$ has no variation to compute a standard error from, and reports the standard error of that $\widehat{\mathrm{ATT}}(s,t)$ from Step 3 instead.

Step 5. The aggregate ATT is a weighted mean of the subaggregate ATTs:

```math
\widehat{\mathrm{ATT}} = \sum^{J}_{j = 1} w_{j}\, \widehat{\mathrm{ATT}}_{j}
```

where $J$ denotes the total number of subaggregate ATTs and the weights, $w_j$, are scaled such that $\sum^{J}_{j = 1} w_{j}\, = 1$. By default $w_j$ is proportional to the number of treated observations behind subaggregate ATT $j$ (see [Weighting](@ref)).

!!! note "Implementation"
    The aggregate ATT is computed as the coefficient in the weighted regression of the subaggregate ATTs on a constant with weights $w_j$, which also yields its standard error. An aggregate ATT with a single subaggregate ATT reports the standard error of that subaggregate ATT. See the [Standard Errors and P-Values](@ref) section.

!!! note "Common Adoption"
    When every treated state is treated at the same time (common adoption) the procedure is exactly the same, with a single treatment cohort. $\hat{\lambda}_{s,t}$ is estimated for every period, the long differences are taken relative to the period immediately before the common treatment time, and every aggregation option is available.

## Common Causal Covariates and Model Specification

A more robust discussion of the common causal covariates assumption can be found in the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447). This section serves only as a quick reference for the different forms of DID-INT. The following table shows the mapping between the DID-INT variation, the functional form of $f\!\left(X^k_{i,s,t}\right)$, and the string value passed to the `ccc` argument in the **DiDInt.jl** package in order to use that form of DID-INT.

| DID-INT Variation            | Functional Form                                                                                                                          | `ccc`     |
| :--------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------- | :---------- |
| Homogeneous DID-INT          | $f(X_{i,s,t}) = \sum_{k=1}^{K} \gamma^k X^k_{i,s,t}$                                                                                   | `"hom"`   |
| Region-Varying DID-INT       | $f(X_{i,s,t}) = \sum_{s=1}^{S} \sum_{k=1}^{K} \gamma^k_s I(s) X^k_{i,s,t}$                                                             | `"state"` |
| Time-Varying DID-INT         | $f(X_{i,s,t}) = \sum_{t=1}^{T} \sum_{k=1}^{K} \gamma^k_t I(t) X^k_{i,s,t}$                                                             | `"time"`  |
| Two-Way DID-INT              | $f(X_{i,s,t}) = \sum_{s=1}^{S} \sum_{t=1}^{T} \sum_{k=1}^{K} \gamma^k_{s,t} I(s)I(t) X^k_{i,s,t}$                                      | `"int"`   |
| Two One-Way DID-INT          | $f(X_{i,s,t}) = \sum_{s=1}^{S} \sum_{k=1}^{K} \gamma^k_s I(s) X^k_{i,s,t} + \sum_{t=1}^{T} \sum_{k=1}^{K} \gamma^k_t I(t) X^k_{i,s,t}$ | `"add"`   |

## Integrated Long Differences

The covariate term in the long differences of Step 2 of the [Five Step Estimation Procedure](@ref), $\left(\hat{\gamma}_{\cdot,t} - \hat{\gamma}_{\cdot,r}\right)'\,\overline{X}_{s,t}$, depends on the DID-INT variation:

- Homogeneous and region-varying DID-INT: the covariate coefficients do not change over time, so the term is zero in every long difference.
- Time-varying and two one-way DID-INT: the term is not zero, but it is the same for the treated state and for each of its controls (they share the time-specific coefficients, and both are evaluated at $\overline{X}_{s,t}$), so it cancels in $\widehat{\mathrm{ATT}}(s,t)$.
- Two-way DID-INT: the coefficients are specific to each state and period, and the term does not cancel.

Consequently the covariate term changes the estimates only when `ccc = "int"` is used with covariates. For every other `ccc` option, and whenever no covariates are specified, $\widehat{\mathrm{ATT}}(s,t)$ only depends on the $\hat{\lambda}_{s,t}$ values:

```math
\widehat{\mathrm{ATT}}(s,t) = \left(\hat{\lambda}_{s,t} - \hat{\lambda}_{s,r}\right) - \frac{\sum_{c \in C(s,t)} v_{c} \left(\hat{\lambda}_{c,t} - \hat{\lambda}_{c,r}\right)}{\sum_{c \in C(s,t)} v_{c}}
```

## Aggregation

In Step 4 of the [Five Step Estimation Procedure](@ref), a choice needs to be made in terms of which $\widehat{\mathrm{ATT}}(s,t)$ values should be averaged together. For example, we could average the $\widehat{\mathrm{ATT}}(s,t)$ values of the states treated at treatment time $g$, which would yield the subaggregate ATT for $g$. Repeating this for each treatment time would give us the set of subaggregate ATTs based on `"cohort"` aggregation. Conversely, we could average the $\widehat{\mathrm{ATT}}(s,t)$ values of a particular treated state, $s^*$, in order to retrieve the subaggregate ATTs per treated state. The following table gives an overview of the aggregation groupings and the string value passed to the `agg` argument in the **DiDInt.jl** package. Every option is available for both staggered adoption and common adoption scenarios.

| Grouping                                                         | $\widehat{\mathrm{ATT}}(s,t)$ values averaged in each subaggregate ATT        | `agg`                |
| :--------------------------------------------------------------- | :---------------------------------------------------------------------------- | :------------------- |
| By treatment time, $g$                                           | Those of the states treated at $g$, in every period                           | `"cohort"` (default) |
| By treated state, $s^*$                                          | Those of state $s^*$, in every period                                         | `"state"`            |
| By treatment time, $g$, and period, $t$                          | Those of the states treated at $g$, in period $t$                             | `"simple"`           |
| By treated state, $s^*$, and period, $t$                         | None: each $\widehat{\mathrm{ATT}}(s,t)$ is reported                          | `"sgt"`              |
| By periods since treatment, $p$                                  | Those with $t - g = p$, across treatment times                                | `"time"`             |
| By period, $t$                                                   | Those in period $t$, across treatment times                                   | `"calendar"`         |
| By policy (and within policy)                                    | See the [Policies](@ref) section                                              | `"policy"`, `"policy_gt"`, `"policy_cohort"`, `"policy_calendar"`, `"policy_time"` |

Whichever option is chosen, the aggregate ATT (`agg_att`) is always reported: it is the weighted mean of the subaggregate ATTs of that option. Because it is a mean of subaggregates, the aggregate ATT can differ across `agg` options. It is the same for every `agg` option when `weighting` is `"both"` (the default) or `"att"` (see [Aggregate ATT Equivalency](@ref)).

## Weighting

An estimate in **DiDInt.jl** is built from three successive averages, and each of them can be weighted by numbers of observations:

1. the long differences of the control states are averaged within each $\widehat{\mathrm{ATT}}(s,t)$ (Step 3 of the [Five Step Estimation Procedure](@ref));
2. the $\widehat{\mathrm{ATT}}(s,t)$ values are averaged within each subaggregate ATT (Step 4);
3. the subaggregate ATTs are averaged within the aggregate ATT (Step 5).

The `weighting` argument decides which averages are weighted: `"diff"` weights the first, `"att"` the second and third, `"both"` (the default) all three, and `"none"` none.

Throughout, $n_{s,t}$ denotes the number of observations of state $s$ in period $t$, and $r$ the period immediately before the treatment time of the treated state in question.

### First Average: Control States Within an ATT(s,t)

Each $\widehat{\mathrm{ATT}}(s,t)$ compares one treated state with its control states $C(s,t)$:

```math
\widehat{\mathrm{ATT}}(s,t) = \widehat{\Delta}_{s,t} - \frac{\sum_{c \in C(s,t)} v_{c}\, \widehat{\Delta}_{c \to s,t}}{\sum_{c \in C(s,t)} v_{c}}
```

with the weight of control state $c$ given by:

```math
v_c = \begin{cases} n_{c,t} + n_{c,r} & \text{if } \texttt{weighting} \in \{\texttt{"both"}, \texttt{"diff"}\} \\ 1 & \text{otherwise} \end{cases}
```

The weights only apply to the control states: there is a single treated state, so nothing is averaged on the treated side. A weighted control state counts in proportion to the number of observations behind its long difference, which are its observations in period $t$ and in period $r$.

### Second Average: ATT(s,t) Values Within a Subaggregate ATT

Each subaggregate ATT averages the $\widehat{\mathrm{ATT}}(s,t)$ values of its group $j$ (see [Aggregation](@ref)):

```math
\widehat{\mathrm{ATT}}_{j} = \frac{\sum_{(s,t) \in j} u_{s,t}\, \widehat{\mathrm{ATT}}(s,t)}{\sum_{(s,t) \in j} u_{s,t}}
```

The weight of the $\widehat{\mathrm{ATT}}(s,t)$ value of treated state $s$ in period $t$ is given by:

```math
u_{s,t} = \begin{cases} n_{s,t} & \text{if } \texttt{weighting} \in \{\texttt{"both"}, \texttt{"att"}\} \\ 1 & \text{otherwise} \end{cases}
```

that is, the number of post-treatment observations of state $s$ in period $t$.

Under `"sgt"` aggregation every subaggregate ATT is a single $\widehat{\mathrm{ATT}}(s,t)$, so this average does not occur.

### Third Average: Subaggregate ATTs Within the Aggregate ATT

The aggregate ATT averages the $J$ subaggregate ATTs:

```math
\widehat{\mathrm{ATT}} = \sum_{j=1}^{J} w_j\, \widehat{\mathrm{ATT}}_j
```

with the weight of subaggregate ATT $j$ given by:

```math
w_j = \begin{cases} \dfrac{N_j}{\sum_{l=1}^{J} N_l} & \text{if } \texttt{weighting} \in \{\texttt{"both"}, \texttt{"att"}\} \\[10pt] \dfrac{1}{J} & \text{otherwise} \end{cases}
```

A weighted subaggregate ATT counts in proportion to $N_j = \sum_{(s,t) \in j} n_{s,t}$, its number of observations from treated states in post-treatment periods. The weights $w_j$ are reported alongside the subaggregate ATTs (the `weight` column of the `sub` table of the results). For the policy aggregation options, the ATT of each policy is formed from the subaggregate ATTs of that policy in the same way, and the pooled ATT across policies has its own weights (see [Policies](@ref)).

#### Aggregate ATT Equivalency

With `weighting = "both"` or `"att"`, the aggregate ATT does not depend on the grouping: since each subaggregate ATT is $\sum_{(s,t) \in j} n_{s,t}\, \widehat{\mathrm{ATT}}(s,t) / N_j$, its weighted mean is $\sum n_{s,t}\, \widehat{\mathrm{ATT}}(s,t) / \sum n_{s,t}$ over all $\widehat{\mathrm{ATT}}(s,t)$ values, whichever `agg` is chosen. The HC standard error and p-value can still differ, since they are computed from different subaggregates; the jackknife standard error and p-value cannot, and neither can the randomization inference p-value except under `"state"` and `"sgt"`. See [this example](@ref "The Aggregate ATT Across Aggregation Options").

## Policies

When treated states adopt different policies, **DiDInt.jl** can report separate ATTs for each policy. The policy of each treated state is given either by the `pvar` argument (the name of a column holding each state's policy, with missing values, `0`, `"0"` or `""` for never-treated states) or by the `treated_policies` argument (a vector with one policy per element of `treated_states`). If neither is specified, every treated state is taken to have the same policy.

Policies do not change which states serve as controls: the controls of a treated state are the never-treated states and, if `notyet = true`, the states not yet treated in that period, whatever their (future) policy. Policies therefore only matter for the policy aggregation options, which form the groups of the [Aggregation](@ref) section within each policy:

| Grouping                                                  | `agg`               |
| :-------------------------------------------------------- | :------------------ |
| By policy                                                 | `"policy"`          |
| By policy, treatment time, and period                     | `"policy_gt"`       |
| By policy and treatment time                              | `"policy_cohort"`   |
| By policy and period                                      | `"policy_calendar"` |
| By policy and periods since treatment                     | `"policy_time"`     |

For these options the ATT of each policy, $\widehat{\mathrm{ATT}}_p$, is the weighted mean of the subaggregate ATTs of that policy, in the same way that the aggregate ATT is the weighted mean of all subaggregate ATTs for the other aggregation options. There is no `agg_att` column for the policy aggregation options.

Setting `policy_aggregation = true` additionally reports a pooled ATT across policies in the `att_pool` column:

```math
\widehat{\mathrm{ATT}}_{pool} = \sum_{p} \omega_p\, \widehat{\mathrm{ATT}}_p, \qquad \sum_p \omega_p = 1
```

with the weight of each policy, reported in the `policy_weights` column, chosen through the `policy_weighting` argument:

| `policy_weighting`  | $\omega_p$ is proportional to...                              | The pooled ATT is the effect for...   |
| :------------------ | :------------------------------------------------------------ | :------------------------------------ |
| `"equal"` (default) | 1                                                             | The typical policy.                   |
| `"treated_obs"`     | The number of treated observations under policy $p$.          | The typical treated observation.      |
| `"treated_cells"`   | The number of $\widehat{\mathrm{ATT}}(s,t)$ values of policy $p$. | The typical treated state-period. |
| `"treated_states"`  | The number of treated states with policy $p$.                 | The typical treated state.            |

The heteroskedasticity-consistent standard error of the pooled ATT comes from the weighted regression of the $\widehat{\mathrm{ATT}}_p$ values on a constant, with degrees of freedom equal to the number of policies minus one (so two policies leave a single degree of freedom). Jackknife and randomization inference results are also reported for each $\widehat{\mathrm{ATT}}_p$ and for the pooled ATT; during the randomization inference procedure a state's policy and treatment time are reassigned together.

## Not Yet Treated Cells

Optionally (not by default), long differences from treated states, prior to the state's treatment time, can be used as control long differences in the Step 3 computation of the $\widehat{\mathrm{ATT}}(s,t)$ values. For example, in a scenario with two treated states and 10 periods, State A being treated at period 2 and State B at period 8, then control long differences can be calculated for State B for periods 2-7 (before B's treatment) to match a subset of the long differences calculated for State A.

Using not-yet-treated long differences as controls is available by setting the `notyet` argument to `true`. A state that is treated in or before period $t$ is never used as a control in period $t$.

## Ensuring Proper Comparisons

In order to best handle data with missing values or inconsistent time period coverage across states, **DiDInt.jl** only computes $\widehat{\mathrm{ATT}}(s,t)$ when the treated state $s$ is observed in both period $t$ and the period immediately before its treatment time, and when at least one control state is observed in both of those periods. Control states that are missing one of the two periods are left out of that $\widehat{\mathrm{ATT}}(s,t)$, and treated state-period pairs without any valid control are dropped. Because every subaggregate and aggregate ATT is a mean of $\widehat{\mathrm{ATT}}(s,t)$ values, each of which compares a treated state with controls in the same pair of periods, no estimate ever compares long differences taken over different periods.

## Period Grid Construction & Date Matching Procedure

**DiDInt.jl** can handle data at various temporal frequencies including daily, weekly, monthly, and yearly observations. The package constructs a period grid and matches observations to discrete periods on that grid.

### Automatic Period Detection

By default, if `freq` is not specified, **DiDInt.jl** automatically detects the appropriate period length from the data by examining the maximum distance between temporally adjacent time observations and identifies whether the data is best represented at a yearly, monthly, weekly, or daily frequency (or some combination thereof) and constructs period lengths accordingly. The automatic detection accounts for leap years.

### Manual Period Specification

Users can manually specify the period length using the `freq` and `freq_multiplier` arguments. For example:

- `freq = "year"` and `freq_multiplier = 1` creates annual periods
- `freq = "month"` and `freq_multiplier = 3` creates quarterly periods
- `freq = "week"` and `freq_multiplier = 2` creates bi-weekly periods

The period grid is constructed starting from `start_date` and ending at `end_date`, creating evenly-spaced periods of the specified length. If the `start_date` and `end_date` arguments were not specified, then the earliest and latest dates in the data are used as the `start_date` and `end_date`, respectively. If `start_date` and `end_date` were specified, then any data falling outside the `start_date` to `end_date` range is dropped from the analysis.

### Date Matching Procedure

Once the period grid is established, each observation in the data is matched to a specific period. The matching procedure follows these rules:

1. Basic matching: An observation at time $t$ is matched to the latest period grid date that is less than or equal to $t$.
2. Treatment boundary adjustment: If there exists a treatment time $T^*$ between two adjacent periods on the grid, observations are matched as follows:

   - If $t < T^*$: the observation is matched to the earlier of the two periods
   - If $t \geq T^*$: the observation is matched to the later of the two periods

   This adjustment ensures that observations are not incorrectly assigned to pre-treatment periods when they occur after treatment begins within that period.
3. Treatment time matching: Treatment times themselves are matched to the first grid period date that is greater than or equal to the treatment time, ensuring treatments are associated with the correct period on the grid.

### Treatment Times Outside the Data

A treated state is part of the analysis as a treated state only if its treatment time falls within the data used, that is, after the first date and no later than the last date (after any restriction through `start_date` and `end_date`). Otherwise, with a warning:

- A state treated at or before the first date is treated throughout the data. It has no pre-treatment period, so no ATT can be computed for it, and it is not a valid control either. Its observations are dropped by `didint()`.
- A state treated after the last date is untreated throughout the data, and is used as a control like a never-treated state.

Each treated state can only have one treatment time: listing a state twice, or a `gvar` column that takes more than one value within a state, is an error.


## Standard Errors and P-Values

**DiDInt.jl** reports three types of p-values and two sets of standard errors, for every subaggregate ATT and for the aggregate ATT (and, for the policy aggregation options, for the ATT of each policy and for the pooled ATT).

- Standard errors using heteroskedasticity-consistent covariance estimators (HCCME): HC0, HC1 (default), HC2, HC3, or HC4 and p-values via a two-tailed test using a $t$-distribution with $n - k$ degrees of freedom, where $n$ is the number of observations in the relevant regression and $k$ is the number of estimated regression parameters. The desired HCCME can be selected via the `hc` argument.
- Jackknife standard errors and p-values, described in more detail in the [Jackknife](@ref) section.
- Randomization inference (RI) p-values, via a randomization procedure described in the [Randomization Inference](@ref) section.

There are some circumstances in which the jackknife standard errors and p-values or the randomization inference p-values cannot be calculated. These situations are noted in the [Jackknife](@ref) and [Randomization Inference](@ref) sections, respectively.

The standard errors using HCCMEs come from the following regressions:

- For each $\widehat{\mathrm{ATT}}(s,t)$: the Step 3 regression of the long differences of the treated state and of its controls on an intercept and a treatment indicator ($n$ is one plus the number of controls and $k = 2$).
- Subaggregate ATTs: the weighted regression of the $\widehat{\mathrm{ATT}}(s,t)$ values of the group on a constant, using the weights $u_{s,t}$ ($n$ is the number of $\widehat{\mathrm{ATT}}(s,t)$ values in the group and $k = 1$). A subaggregate ATT that consists of a single $\widehat{\mathrm{ATT}}(s,t)$ reports the standard error and p-value of that $\widehat{\mathrm{ATT}}(s,t)$.
- Aggregate ATT: the weighted regression of the subaggregate ATTs on a constant, using the weights $w_j$ ($n$ is the number of subaggregate ATTs and $k = 1$). If there is a single subaggregate ATT, the aggregate ATT reports its standard error and p-value.

!!! note "Caveats for the HCCME Standard Errors"
    The regressions on a constant treat the estimates being averaged as independent. $\widehat{\mathrm{ATT}}(s,t)$ values from the same period share control states, so these standard errors can understate the uncertainty of subaggregate and aggregate ATTs. The jackknife and the randomization inference procedure do not rely on that independence.

    The single treated state in the Step 3 regression has a leverage of exactly one. HC2, HC3, and HC4 divide by one minus the leverage, so the HC2, HC3, and HC4 standard errors of an $\widehat{\mathrm{ATT}}(s,t)$ (reported under `"sgt"` aggregation, and for any subaggregate ATT that consists of a single $\widehat{\mathrm{ATT}}(s,t)$) are `NaN`. HC0 and HC1 are unaffected.

When an $\widehat{\mathrm{ATT}}(s,t)$ has only one control state, the Step 3 regression has two observations and two parameters: the model is saturated and the regression-based standard error is degenerate. In such cases, it is still possible to recover the heteroskedastic robust standard error using the variances and covariances of the $\hat{\lambda}_{s,t}$ terms. Given one treated state $s$ and one control state $c$, the variance of $\widehat{\mathrm{ATT}}(s,t)$ (without the covariate term of the [Integrated Long Differences](@ref)) can be decomposed into variance and covariance terms of the $\hat{\lambda}_{s,t}$ values:

```math
\begin{aligned}
\widehat{\mathrm{ATT}}(s,t)
&=
(\hat{\lambda}_{s,t} - \hat{\lambda}_{s,r})
-
(\hat{\lambda}_{c,t} - \hat{\lambda}_{c,r})
\\[6pt]
\mathrm{Var}(\widehat{\mathrm{ATT}}(s,t))
&=
\mathrm{Var}(\hat{\lambda}_{s,t} - \hat{\lambda}_{s,r})
+
\mathrm{Var}(\hat{\lambda}_{c,t} - \hat{\lambda}_{c,r})
\\
&\quad
-2\,\mathrm{Cov}\!\left(
\hat{\lambda}_{s,t} - \hat{\lambda}_{s,r},
\hat{\lambda}_{c,t} - \hat{\lambda}_{c,r}
\right)
\\[6pt]
\mathrm{Var}(\widehat{\mathrm{ATT}}(s,t))
&=
\mathrm{Var}(\hat{\lambda}_{s,r})
+
\mathrm{Var}(\hat{\lambda}_{s,t})
+
\mathrm{Var}(\hat{\lambda}_{c,r})
+
\mathrm{Var}(\hat{\lambda}_{c,t})
\\
&\quad
-2\,\mathrm{Cov}(\hat{\lambda}_{s,r}, \hat{\lambda}_{s,t})
-2\,\mathrm{Cov}(\hat{\lambda}_{c,r}, \hat{\lambda}_{c,t})
\\
&\quad
-2\Big(
\mathrm{Cov}(\hat{\lambda}_{s,t}, \hat{\lambda}_{c,t})
-
\mathrm{Cov}(\hat{\lambda}_{s,t}, \hat{\lambda}_{c,r})
-
\mathrm{Cov}(\hat{\lambda}_{s,r}, \hat{\lambda}_{c,t})
+
\mathrm{Cov}(\hat{\lambda}_{s,r}, \hat{\lambda}_{c,r})
\Big)
\end{aligned}
```

Further details on recovering the standard errors in such an edge case can be found in the section [Computation of Edge Case Standard Errors](@ref).

## Jackknife

There are two variations to the **DiDInt.jl** jackknife procedure. The default option is to calculate jackknife estimates of the ATTs while dropping one state at a time directly from Step 3 in the [Five Step Estimation Procedure](@ref), and the $\widehat{\mathrm{ATT}}(s,t)$ values, the subaggregate ATTs, and the aggregate ATT are recomputed from the $\hat{\lambda}_{s,t}$ values of the remaining states. This is referred to as the "fast jackknife". The alternative method, the "true jackknife", calculates jackknife estimates of the ATTs while dropping one state at a time starting from Step 1. The true jackknife can be used in lieu of the fast jackknife by setting the `truejack` argument to `true`.

Note that for the `ccc` options of `"state"` and `"int"` that the true jackknife is always equal to the fast jackknife as the covariate effects are local to the state or state-time, respectively, so that removing a state from the dataset does not affect the estimated covariate effects at other states and thus does not affect the estimated $\hat{\lambda}_{s,t}$ values at other states. For `"hom"`, `"time"` and `"add"`, the covariate effects are estimated across states, so removing a state changes the $\hat{\lambda}_{s,t}$ of the others, and the true jackknife must start over from Step 1.

The jackknife standard error is computed as:

```math
\widehat{\mathrm{se}}_{\mathrm{jack}} = \sqrt{\frac{N-1}{N} \sum_{i=1}^{N}\left(\hat{\beta}_{(-i)} - \overline{\beta}\right)^2}, \qquad \overline{\beta} = \frac{1}{N} \sum_{i=1}^{N} \hat{\beta}_{(-i)}
```

where $N$ is the number of states whose removal leaves the estimate defined, $\hat{\beta}_{(-i)}$ is the ATT estimate with state $i$ removed, and $\overline{\beta}$ is the mean of those estimates. This is the cluster jackknife variance estimator centered at the mean of the delete-one estimates ($\mathrm{CV}_{3J}$ in [MacKinnon, Nielsen, and Webb (2023)](https://doi.org/10.1016/j.jeconom.2022.04.001); equation 87 of the DID-INT paper). Jackknife p-values are computed via a two-tailed test of the original ATT estimate using a $t$-distribution with $N - 1$ degrees of freedom.

A minimum of two treated states and two control states are required to compute jackknife standard errors and p-values; if this condition is not met, these quantities are returned as `missing`. The same minimum applies to each estimate separately: a subaggregate ATT (or the ATT of a policy) only has a jackknife standard error if the $\widehat{\mathrm{ATT}}(s,t)$ values it averages involve at least two treated states and at least two control states, and if it can be computed with any one state removed. Otherwise its jackknife standard error and p-value are returned as `missing`.

The states that are removed in turn are those that contribute at least one long difference to the estimates; a state that never does (for example one that is only observed before every treatment time) is not counted. If removing a state leaves a subaggregate ATT without any treated state (for example the only state of a treatment cohort), that subaggregate ATT has no jackknife standard error (it is `missing`), and the aggregate ATT for that removal is the mean over the remaining subaggregate ATTs. The jackknife standard error of the aggregate ATT is still reported. The same applies to the pooled ATT across policies when removing a state leaves a policy without any treated state.

Further, note that for the aggregation methods of `"state"` and `"sgt"`, where each subaggregate ATT involves only one treated state, that the jackknife standard errors and jackknife p-values are never estimated at the subaggregate level: removing that one treated state leaves nothing to estimate. This limitation does not prevent the jackknife standard error and jackknife p-value from being calculated for the *aggregate* ATT under `"state"` or `"sgt"` aggregation, so long as there are at least two treated and two control states.

## Randomization Inference

The randomization procedure used in **DiDInt.jl** follows the procedure described in [MacKinnon and Webb (2020)](https://doi.org/10.1016%2Fj.jeconom.2020.04.024), although for the `"state"` and `"sgt"` aggregation options, there is an additional randomization step.

The **DiDInt.jl** implementation of the randomization inference procedure is as follows:

1. The total number of unique alternative treatment assignments is calculated by $\frac{N!}{(N-n)! \prod_{m \in M} n_m!} - 1$, where $N$ is the total number of states, $n$ is the number of treated states, $M$ is the set of unique treatment times, and $n_m$ is the number of states treated at time $m$. (For the policy aggregation options, $M$ is the set of unique combinations of treatment time and policy.) The number of permutations considered during the procedure, $nperm^{*}$, is the lesser of this total and the user-supplied `nperm` value.
2. If the total does not exceed `nperm`, every alternative assignment is used. Otherwise, $nperm^{*}$ distinct alternative assignments are drawn at random, each being equally likely. The assignments are drawn without repetition and never include the actual assignment; the `seed` argument makes the draws reproducible.
3. Subaggregate and aggregate ATTs are then calculated for each alternative assignment, as if the states that received the treatment times in that assignment had been the treated states. The $\hat{\lambda}_{s,t}$ values do not depend on which states are treated, so Step 1 is not repeated, and the estimates for all alternative assignments are computed together (see [Fast Computation of Randomization Inference and Jackknife](@ref)). The distribution of estimates resulting from the alternative assignments is compared against the estimates computed when using the actual assignment of treatment times:

```math
p_{RI} = \frac{1 + \sum_{r=1}^{nperm^{*}} \mathbf{1}\left(|ATT_{RI,r}| \geq |ATT_{Actual}|\right)}{nperm^{*} + 1}
```

This is the more conservative of the two p-values discussed in [MacKinnon and Webb (2020)](https://doi.org/10.1016%2Fj.jeconom.2020.04.024). An RI p-value is reported for every subaggregate ATT, each compared against the estimates for the same group under the alternative assignments, and for the aggregate ATT.

!!! note "Drawing Alternative Assignments"
    The method for drawing alternative assignments is essentially a combination of two procedures: Floyd's algorithm for sampling distinct integers, and the problem of finding the $k$-th permutation of a lexicographically ordered set. The latter is used in its multiset version, which is adjusted for values that occur more than once in the set, since several states can share a treatment time and several states are never treated.

    Floyd's algorithm ([Bentley and Floyd, 1987](https://doi.org/10.1145/30401.315746)) draws distinct whole numbers at random, without replacement. Each number $k$ is mapped to the $k$-th assignment in a lexicographically ordered list of every distinct treatment assignment, using the procedure for finding the $k$-th permutation of an ordered multiset; the list itself is never created. The actual assignment is excluded by computing its position in that list and adding 1 to every drawn number at or above that position.

Alternative assignments are evaluated in batches of 500, so memory use does not grow with `nperm`.

For the `"state"` and `"sgt"` options, each treated state's treatment time is assigned to a randomly drawn state in each alternative assignment, and that state's estimate is used in its place.

### Kernel-Smoothed P-Values

With few alternative assignments, the p-value above can only take a few distinct values. When the number of permutations used, $nperm^{*}$, is smaller than the `ri_kernel_cutoff` argument (99 by default), **DiDInt.jl** instead reports the kernel-smoothed p-value of [Racine and MacKinnon (2007)](https://doi.org/10.1016/j.csda.2006.11.013), which replaces the indicator function with a smooth function:

```math
p_{RI} = 1 - \frac{1}{nperm^{*}} \sum_{r=1}^{nperm^{*}} \Phi\!\left(\frac{|ATT_{Actual}| - |ATT_{RI,r}|}{h}\right), \qquad h = 1.575\; \hat{\sigma}_{RI}\; (nperm^{*})^{-4/9}
```

where $\Phi$ is the standard normal cumulative distribution function and $\hat{\sigma}_{RI}$ is the sample standard deviation of the $ATT_{RI,r}$ values. This applies whether $nperm^{*}$ is small because few alternative assignments exist or because a small `nperm` was requested. Kernel smoothing can be turned off by setting `ri_kernel_cutoff` to `nothing`, `false`, or `0`. If there is only one alternative assignment, or if the $ATT_{RI,r}$ values do not vary, the unsmoothed p-value is reported.

### Ties

Although every alternative assignment differs from the actual assignment, an estimate under an alternative assignment can be exactly equal to the actual estimate. Assignments are unique as a whole, but a given subaggregate ATT only depends on part of the assignment: for example, the subaggregate ATT for treatment time $g$ and period $t$ only depends on which states are treated at $g$ and which states are controls in period $t$. An alternative assignment that reshuffles other states leaves that subaggregate ATT unchanged. Such ties are most common when `notyet = true` for early treatment cohorts with few states, as the set of not-yet-treated controls in an early period changes little across assignments. The aggregate ATT depends on the entire assignment and almost never ties.

Ties are counted in favour of the null hypothesis (they satisfy $|ATT_{RI,r}| \geq |ATT_{Actual}|$). Because two estimates that are equal in exact arithmetic can differ by rounding error, estimates within a small relative tolerance (about $1.5 \times 10^{-8}$) of the actual estimate are treated as ties.

This applies to the unsmoothed p-value. With [Kernel-Smoothed P-Values](@ref), an estimate equal to the actual estimate contributes $\Phi(0) = 0.5$ to the sum, so a tie counts as half of an exceedance rather than a whole one.

### When Randomization Inference Is Not Performed

It should also be noted that the randomization inference procedure will not run unless it is possible to calculate the long differences for each $(g,t)$ used in the computation of the actual ATT values for *every* state, treated or not: under an alternative assignment any state can take the place of a treated state. If a state is not observed in one of the periods involved (which may occur in datasets with missing values or inconsistent time period coverage across states), the randomization inference p-values are returned as `missing`, with a warning. This safeguard is to ensure that the randomization procedure only runs if it is possible to make a proper comparison between the actual $ATT$ values and the $ATT_{RI}$ values.

## Estimating Lambda

As mentioned previously, estimating the regression model shown in Step 1 would be computationally infeasible, especially for a dataset with many states, time periods, and covariates. For example, the number of columns in the design matrix for a two-way DID-INT model would be $n_{s} \cdot n_{t} + n_{s} \cdot n_{t} \cdot n_{k}$ where $n_{s}$ is the number of states, $n_{t}$ is the number of time periods, and $n_{k}$ is the number of covariates. Consequently, **DiDInt.jl** takes three different approaches for estimating the $\hat{\lambda}_{s,t}$ values depending on the specific DID-INT variation. For each approach, the resulting $\hat{\lambda}_{s,t}$ values are numerically equivalent to those from the combined regression formulation shown in Step 1.

In each case the estimated covariate coefficients and the covariate means of each $(s,t)$ cell are retained, as they are needed for the [Integrated Long Differences](@ref).

### Two-Way DID-INT

Two-way DID-INT is the simplest case. As the effects of the covariates are local to each $(s,t)$ group, separate regressions of the form:

```math
Y_{i,s,t} = \lambda_{s,t} + \sum_{k=1}^{K} \gamma^k_{s,t} X^k_{i,s,t} + \epsilon_{i,s,t}
```

are estimated for each $(s,t)$ cell, where $\hat{\lambda}_{s,t}$ is simply the estimated intercept for the local regression. 

Covariates that are redundant in the regression of a cell are left out of it and have a coefficient of zero there (see [Rank Deficiency and Lambda Identification](@ref)).

### Region-Varying, Time-Varying, and Homogeneous DID-INT

For region-varying, time-varying, and homogeneous DID-INT, the Frisch-Waugh-Lovell (FWL) theorem is applied. The procedure is the same across all three cases, differing only in the grouping over which the covariate coefficients $\hat{\beta}$ are estimated:

1. Both the outcome and covariates are demeaned within each $(s,t)$ cell - which is equivalent to residualizing the outcome and covariates by the $I(s,t)$ cell dummies.
2. The covariate coefficients $\hat{\beta}$ are estimated by regressing the demeaned outcome on the demeaned covariates (without an intercept), where the regression is run:
   - Globally across all observations for homogeneous DID-INT
   - Within each time period $t$ for time-varying DID-INT
   - Within each state $s$ for region-varying DID-INT

   Covariates that are redundant in a regression are left out of it and have a coefficient of zero there (see [Rank Deficiency and Lambda Identification](@ref)).
3. The $\hat{\lambda}_{s,t}$ values are then recovered as:

```math
\hat{\lambda}_{s,t} = \overline{Y}_{s,t} - \hat{\beta}'\overline{X}_{s,t}
```

where $\overline{Y}_{s,t}$ and $\overline{X}_{s,t}$ are the cell means of the outcome and covariates for state $s$ in period $t$, respectively. This follows from taking expectations of the regression model conditional on $(s,t)$:

```math
\mathbb{E}[Y_{i,s,t} \mid s,t] = \lambda_{s,t} + \beta'\mathbb{E}[X_{i,s,t} \mid s,t]
```

which gives $\overline{Y}_{s,t} = \hat{\lambda}_{s,t} + \hat{\beta}'\overline{X}_{s,t}$, and rearranging yields the recovery formula above.

### Two One-Way DID-INT

Two one-way DID-INT follows the same FWL-based procedure as the preceding three cases, differing only in how step 2 is structured to reflect the additive nature of two one-way DID-INT:

1. Both the outcome and covariates are demeaned within each $(s,t)$ cell - again, equivalent to residualizing by the $I(s,t)$ cell dummies.
2. Because two one-way DID-INT implies that covariate effects are additively separable into a state-specific component, $\hat{\beta}_s$, and a time-specific component, $\hat{\beta}_t$, both sets of coefficients must be estimated jointly. This is done by regressing the demeaned outcome $\tilde{Y}_{i,s,t}$ (where $\tilde{Y}_{i,s,t} = Y_{i,s,t} - \overline{Y}_{s,t}$) on a structured sparse design matrix $Z$ of dimension $n \times (n_s + n_t)K$. The matrix $Z$ consists of $n_s + n_t$ blocks of $K$ columns each - one block per state, followed by one block per time period. For observation $i$ in cell $(s,t)$, the demeaned covariates $\tilde{X}_{i,s,t}$ are placed in the $s$-th state block and the $t$-th time block, with all other entries in that row set to zero. The regression model used to recover $\hat{\beta}_s$ and $\hat{\beta}_t$ is given by: $\tilde{Y}_{i,s,t} = \hat{\beta}'Z_{i,s,t} + \varepsilon_{i,s,t}$, where the first $n_{s} \cdot K$ entries of $\hat{\beta}$ correspond to the $\hat{\beta}_s$ coefficients, and the following $n_{t} \cdot K$ entries of $\hat{\beta}$ correspond to the $\hat{\beta}_t$ coefficients.
3. The $\hat{\lambda}_{s,t}$ values are then recovered as:

$\hat{\lambda}_{s,t} = \overline{Y}_{s,t} - \hat{\beta}_s'\overline{X}_{s,t} - \hat{\beta}_t'\overline{X}_{s,t}$

### Rank Deficiency and Lambda Identification

When the covariates are linearly dependent after demeaning within cells (a cell being an $(s,t)$ group), the $\hat{\lambda}_{s,t}$ values are not uniquely determined by the data. The ATT is unaffected only when the resulting ambiguity in $\hat{\lambda}_{s,t}$ is additive in state and period, because such shifts cancel in the long differences; collinearity among covariates in general does not have this form.

Under `"hom"`, `"state"`, `"time"` and `"int"`, **DiDInt.jl** therefore drops redundant covariates within each block of the data in which the covariate coefficients are estimated (all data, each state, each period, or each cell). It goes through the covariates in the order given in `covariates` and drops any covariate that is, up to a relative tolerance, a linear combination of those already kept: a covariate is dropped when less than $10^{-7}$ of it (measured by its norm, after demeaning within cells) is left once the part explained by the covariates kept before it is removed. The covariates are never reordered, so the estimates are unique, unchanged when covariates are rescaled, and identical across procedures. A covariate that does not vary within any cell of a block is dropped as well.

Collinearity is checked on the covariates centered within each cell, because the regression has an intercept: a covariate can be collinear with the intercept (for example, constant in the cell) without being collinear with the other covariates. Under `"int"`, another case is a categorical covariate whose reference category does not occur in a cell: the cell's dummy variables add up to one, and the dummy variable listed last is dropped there. A cell with $n$ observations keeps at most $n - 1$ covariates.

Dropping collinear covariates fixes a convention: the data of a cell cannot tell their coefficients apart, so **DiDInt.jl** drops them by a fixed rule that depends only on the order of the `covariates` argument and sets their coefficients to zero. Under `ccc = "int"` this can change an $\widehat{\mathrm{ATT}}(s,t)$, but only when a cell's coefficients are evaluated at covariate values outside that cell's data, which happens for the treated state's base period and for the control cells, all evaluated at $\overline{X}_{s,t}$. **DiDInt.jl** warns in exactly that case, giving the number of $\widehat{\mathrm{ATT}}(s,t)$ values affected and the covariates responsible. Under `"hom"`, `"state"` and `"time"` the same rule is used, and a warning names every covariate dropped for being collinear, since the estimates can depend on which one is dropped (which depends on the order of `covariates`); a covariate dropped only because it does not vary within any cell of a block is not reported.

Under `"add"`, each cell's covariate coefficients are $\gamma_s + \gamma_t$, so adding a constant to every state's coefficients and subtracting it from every period's leaves every $\hat{\lambda}_{s,t}$ and ATT unchanged; this rank deficiency is always present and harmless. **DiDInt.jl** solves the regression with Julia's `\` operator and does not drop covariates. Further rank deficiencies leave the $\hat{\lambda}_{s,t}$ not unique. Some still leave the ATTs unique, for example a covariate that is constant over time within a state, whose effect shifts every $\hat{\lambda}_{s,t}$ of that state by the same amount and cancels in its long differences. Others can make the ATTs depend on the solution. **DiDInt.jl** checks every such direction and warns when the ATTs are affected. Edge case standard errors are not available under `"add"`.

The warning reports the rank of the design matrix $Z$ next to the rank that the split between states and periods leaves. With $K$ covariates, $n_s$ states and $n_t$ periods, $Z$ has $K (n_s + n_t)$ columns: one coefficient per covariate for each state and for each period. The split costs one rank per covariate, so the highest rank that $Z$ can have is

```math
K \left(n_s + n_t - C\right)
```

where $C$ is the number of separate groups that the panel falls into, a group being a set of states and periods that are linked by observations ($C = 1$ unless some states are only observed in periods in which none of the other states are observed). Full column rank is never reached, and a further rank deficiency exists when the rank is below this value. For example, with 20 states, 12 periods and 3 covariates, $Z$ has $3 \times 32 = 96$ columns and the expected rank is $3 \times (20 + 12 - 1) = 93$. If one of the covariates does not vary within any cell of one particular state, that state's coefficient for that covariate is not determined by the data and the rank is 92. If an ATT is affected (the covariate's value in that state changes over time), the warning reads: *"Under ccc = "add" the covariate slopes are less well determined than usual: the first-stage regression has rank 92 where 93 is expected. ..."*

## Computation of Edge Case Standard Errors

In general, the heteroskedasticity-robust standard error of an $\widehat{\mathrm{ATT}}(s,t)$ is retrieved directly from the Step 3 regression. However, as noted in the [Standard Errors and P-Values](@ref) section, when the treated state has only one control state in that period the regression has two observations and two parameters, the model is saturated, and the regression-based standard error is degenerate. In such cases, the variance of the ATT can instead be recovered from the variances and covariances of the $\hat{\lambda}_{s,t}$ values, as shown in the decomposition in the aforementioned section. This is referred to as the edge case.

The edge case is a property of an individual $\widehat{\mathrm{ATT}}(s,t)$, whatever the aggregation method and whether adoption is common or staggered. The edge case standard errors therefore appear in the results:

- under `"sgt"` aggregation, for every $\widehat{\mathrm{ATT}}(s,t)$ with a single control state;
- for any subaggregate ATT that consists of a single $\widehat{\mathrm{ATT}}(s,t)$ with a single control state (and for the aggregate ATT if that is the only subaggregate ATT).

Subaggregate ATTs that average two or more $\widehat{\mathrm{ATT}}(s,t)$ values obtain their standard error from the spread of those values, as described in the [Standard Errors and P-Values](@ref) section, and do not involve the edge case. The jackknife and the randomization inference procedure only use point estimates and are never affected. An edge case standard error is reported without a p-value, as the saturated regression leaves no residual degrees of freedom for a $t$-distribution.

Recovering the variance of the ATT from the variances and covariances of the $\hat{\lambda}_{s,t}$ values necessitates computing the variance-covariance matrix of the $\hat{\lambda}_{s,t}$ values, $\widehat{\mathrm{Var}}(\hat{\lambda})$. This is a more computationally expensive procedure than the standard lambda estimation procedure, so by default **DiDInt.jl** does not compute it. Computation of edge case standard errors can be toggled on by setting `edgecase = true`. When edge case computation is toggled on, the set of required $\hat{\lambda}_{s,t}$ values needed to compute ATTs is determined prior to any computations, and those involved in the saturated Step 3 regression(s) are identified. Only the necessary blocks of data are sent to the edge case computation procedure, where both the $\hat{\lambda}_{s,t}$ values and their variances/covariances are recovered. The separation of the data into blocks sent to the edge case procedure, and the rest of the data which is fed into the standard [Estimating Lambda](@ref) procedure is determined by the $(s,t)$ groups to which the implicated $\hat{\lambda}_{s,t}$ belong to as well as the specific DID-INT variation. For example, time-varying DID-INT requires that all observations in group $t$ are grouped together, so any observation belonging to a $t$ period implicated in a saturated Step 3 regression is sent to the edge case procedure. Parallel logic applies to the other DID-INT variations, however it should be noted that the edge case procedure is not an option for two one-way DID-INT as this model is always rank deficient and $\widehat{\mathrm{Var}}(\hat{\lambda})$ is not recoverable. 

!!! note "Consistency with Lambda Estimation"
    The regular lambda estimation procedure and the edge case procedure use exactly the same covariates in every block. Covariates that do not vary within any $(s,t)$ cell of a block are absorbed by the cell dummies and dropped. Under `"hom"`, `"state"`, `"time"` and `"int"`, covariates that are collinear with other covariates in the block are then dropped by the rule described in [Rank Deficiency and Lambda Identification](@ref). The covariates kept are determined once per block and used by both procedures, so the two return identical $\hat{\lambda}_{s,t}$ values and ATT estimates.

### Computing $\widehat{\mathrm{Var}}(\hat{\lambda})$

The variance-covariance matrix $\widehat{\mathrm{Var}}(\hat{\lambda})$ is an $n_\lambda \times n_\lambda$ matrix, where $n_\lambda = n_s \cdot n_t$ is the total number of $\hat{\lambda}_{s,t}$ values. For homogeneous DID-INT it is computed by following a single regression, and for all other DID-INT variations it is computed in blocks following regressions on subsets of the data. In each block regression, the raw outcome $Y_{i,s,t}$ is regressed on the $I(s,t)$ cell dummies of the block and on the original (non-demeaned) covariates that the block keeps (see [Rank Deficiency and Lambda Identification](@ref)), after which the covariance matrix of the $\hat{\lambda}_{s,t}$ values from that regression is recovered using the chosen HCCME.

The specific regression used to compute $\widehat{\mathrm{Var}}(\hat{\lambda})$ mirrors the CCC structure of the chosen DID-INT variant, which also determines which of the terms of the variance decomposition in the [Standard Errors and P-Values](@ref) section can be nonzero:

- Two-way (`"int"`): a separate regression is run for each $(s,t)$ cell, with the cell dummy and covariates as regressors. The HCCME for each cell's dummy coefficient is placed into the corresponding diagonal entry of $\widehat{\mathrm{Var}}(\hat{\lambda})$. Off-diagonal entries between different cells are zero since the regressions are run separately, so only the four variance terms of the decomposition remain.
- Region-varying (`"state"`): a separate regression is run for each state $s$, pooling all time periods within that state. The HCCME for the cell dummies within each state block is placed into the corresponding entries of $\widehat{\mathrm{Var}}(\hat{\lambda})$ (including off-diagonal covariance terms). The covariances between the two periods of the same state are nonzero; those between the treated state and the control state are zero.
- Time-varying (`"time"`): analogous to region-varying, with separate regressions run for each time period $t$. The covariances between the treated state and the control state in the same period are nonzero; all others are zero.
- Homogeneous (`"hom"`): a single global regression is run across all observations. The HCCME for all cell dummies from the regression is exactly $\widehat{\mathrm{Var}}(\hat{\lambda})$, and every term of the decomposition can be nonzero.

Under HC1, each block regression's HCCME is scaled by $n_b / (n_b - k_b)$, where $n_b$ and $k_b$ are that regression's number of observations and number of parameters (cell dummies plus covariates), not those of the full sample.

Because the block regressions only include the covariates that the block keeps, their design matrices have full column rank, and rank deficiency is not a reason for a `missing` edge case standard error. An edge case standard error is `missing` in the following cases only:

- one of the block regressions involved has no residual degrees of freedom (it has as many parameters as observations);
- the estimated variance of one of the four $\hat{\lambda}_{s,t}$ values is exactly zero (for example a binary outcome that takes a single value in one of the cells);
- under HC2, HC3, and HC4, a block regression includes a cell with a single observation (that observation has a leverage of one).

### Edge Case Standard Errors for Two-Way DID-INT with Covariates

For homogeneous, region-varying, and time-varying DID-INT, the covariate term of the [Integrated Long Differences](@ref) cancels, $\widehat{\mathrm{ATT}}(s,t)$ is a combination of four $\hat{\lambda}_{s,t}$, and the decomposition above gives its variance.

Under `ccc = "int"` with covariates, the covariate term of the [Integrated Long Differences](@ref) does not cancel. Four cells enter the estimate, $(s,t)$, $(s,r)$, $(c,t)$ and $(c,r)$, each evaluated at the covariate means of the treated cell, $\bar{x} = \overline{X}_{s,t}$:

```math
\widehat{\mathrm{ATT}}(s,t) = \left(\hat{\lambda}_{s,t} + \hat{\gamma}_{s,t}'\bar{x}\right) - \left(\hat{\lambda}_{s,r} + \hat{\gamma}_{s,r}'\bar{x}\right) - \left(\hat{\lambda}_{c,t} + \hat{\gamma}_{c,t}'\bar{x}\right) + \left(\hat{\lambda}_{c,r} + \hat{\gamma}_{c,r}'\bar{x}\right)
```

Conditional on the observed covariates, the four cells' regressions are independent (see [Independence Across Cells](@ref)), so, with $k$ running over the four cells,

```math
\begin{aligned}
\mathrm{Var}\!\left(\widehat{\mathrm{ATT}}(s,t) \mid X\right)
&= \sum_k \mathrm{Var}\!\left(\hat{\lambda}_k + \hat{\gamma}_k'\bar{x} \mid X\right) \\
&= \sum_k \left[\mathrm{Var}(\hat{\lambda}_k) + 2\,\bar{x}'\,\mathrm{Cov}(\hat{\gamma}_k, \hat{\lambda}_k) + \bar{x}'\,\mathrm{Var}(\hat{\gamma}_k)\,\bar{x}\right]
\end{aligned}
```

with the variances and covariances estimated from each cell's regression. This treats $\overline{X}_{s,t}$ as fixed and so omits its sampling variability, a reasonable simplification for an edge case that does not report p-values. Without covariates it reduces to $\sum_k \mathrm{Var}(\hat{\lambda}_k)$.

### Independence Across Cells

!!! warning "Panel Data"
    Every edge case standard error assumes that observations in different $(s,t)$ cells are independent of each other. This holds for repeated cross-sections, where each period samples different individuals. It does not hold for panel data, where the same individuals are observed in several periods: the estimates for a state's periods are then correlated, and the formulas above omit that covariance. With panel data the edge case standard errors are generally biased.

## Plotting with `didint_plot()`

`didint_plot()` prepares the data for five kinds of plots, which `plot` draws; the [Plotting Examples](@ref) page shows each of them. All are built from the first-stage estimates $\hat{\lambda}_{s,t}$ of Step 1 of the [Five Step Estimation Procedure](@ref), estimated on all states exactly as in `didint()`. The `plot_type` argument of `didint_plot()` chooses the kinds of data to prepare (by default all of them), and each kind is returned as a table of the same name.

Below, $g_s$ is the first treated period of a treated state $s$, and $r_s = g_s - 1$ its pre-treatment period. The control states $C_{s,t}$ of $s$ in period $t$ are the never-treated states and, with `notyet = true`, the states of later cohorts that are not yet treated in $t$, as in `didint()`.

### Conditional Parallel Trends

The `"conditional"` data are the $\hat{\lambda}_{s,t}$ values themselves, for every state and period, under every `ccc` option other than `"int"`:

```math
\hat{\lambda}_{s,t} = \overline{Y}_{s,t} - \hat{\gamma}_{s,t}'\,\overline{X}_{s,t}
```

where $\hat{\gamma}_{s,t}$ are the covariate coefficients that the `ccc` option implies for the cell (`"none"` uses $\overline{Y}_{s,t}$ itself). Under these options `didint()` compares exactly the changes over time in these series, so series that move in parallel before treatment support conditional parallel trends under that option, and comparing the options shows which form of covariate adjustment makes the trends parallel.

The plot has one panel per `ccc` option, with the series of every state: each treated state in its own colour, the never-treated states in grey. With `facet = "state"` it instead has one panel per treated state $s$ and `ccc` option, which compares the series of $s$ with the average of its control states (or, with `controls = "separate"`, with each of them). Because `didint()` compares the *change* of a treated state with the average of its controls' changes, the line of the average control state is built from changes. It starts at the average of the controls' $\hat{\lambda}$ in the first period, and then

```math
\overline{\lambda}^{C}_{s,t} = \overline{\lambda}^{C}_{s,t-1} + \sum_{c \in C_{s,t}} w_{c,t} \left(\hat{\lambda}_{c,t} - \hat{\lambda}_{c,t-1}\right), \qquad w_{c,t} = \frac{n_{c,t} + n_{c,t-1}}{\sum_{c' \in C_{s,t}} \left(n_{c',t} + n_{c',t-1}\right)}
```

with $w_{c,t} = 1 / |C_{s,t}|$ under `weights = false`. With the same control states and weights in every period this is the average of the controls' $\hat{\lambda}_{c,t}$. An average of those levels would move whenever a not-yet-treated state leaves the controls or the numbers of observations change, even if every control state were parallel to the treated state; this line does not.

### Changes Since the Pre-Treatment Period

The `"conditional2"` data are the differences from which `didint()` computes its estimates under the `ccc` options other than `"int"`, for every period instead of the treated periods only. For a treated state $s$ and a state $j$ (the state $s$ itself or one of its control states),

```math
\widehat{\delta}_{j,t} = \hat{\lambda}_{j,t} - \hat{\lambda}_{j,r_s}
```

for every period $t \neq r_s$; in period $r_s$ the change is zero. The plot has one panel per treated state and `ccc` option, which compares $\widehat{\delta}_{s,t}$ with the average over its control states,

```math
\overline{\delta}_{s,t} = \sum_{c \in C_{s,t}} w_{c,t}\, \widehat{\delta}_{c,t}, \qquad w_{c,t} = \frac{n_{c,t} + n_{c,r_s}}{\sum_{c' \in C_{s,t}} \left(n_{c',t} + n_{c',r_s}\right)}
```

(or, with `controls = "separate"`, with each $\widehat{\delta}_{c,t}$). These are the differences and the weights of Steps 2 and 3 of the [Five Step Estimation Procedure](@ref), so the gap between the two lines is the estimate itself:

```math
\widehat{\delta}_{s,t} - \overline{\delta}_{s,t} = \widehat{ATT}(s,t) \qquad \text{for } t \geq g_s
```

For $t < r_s$ the gap is a placebo estimate (see [Placebo Estimates Before Treatment](@ref)).

### Integrated Parallel Trends

Under `"int"` the $\hat{\lambda}_{s,t}$ values depend on how the covariates are coded, so the `"integrated"` data are instead those of the integrated parallel trends plot of the DID-INT paper (equations 36 to 38). The outcome is regressed on the covariates separately in every $(s,t)$ cell, which gives $\hat{\lambda}_{s,t}$ and $\hat{\gamma}_{s,t}$. For a treated state $s$ and a state $j$, the change in the predicted outcome of $j$ between periods $t-1$ and $t$, evaluated at the covariate means of $s$ in period $t$, is

```math
\widehat{\delta}_{j \to s,t} = \left(\hat{\lambda}_{j,t} - \hat{\lambda}_{j,t-1}\right) + \left(\hat{\gamma}_{j,t} - \hat{\gamma}_{j,t-1}\right)'\,\overline{X}_{s,t}
```

It is computed for every treated state $s$ (the `treated_state_covars` column), for $j = s$ and for each control state of $s$. The `n` column holds $n_{j,t} + n_{j,t-1}$.

These are the [Integrated Long Differences](@ref) of `didint()` under `ccc = "int"`, taken from one period to the next rather than from the pre-treatment period; in the first treated period the two coincide. The plot has one panel per treated state, which compares $\widehat{\delta}_{s \to s,t}$ with the average of $\widehat{\delta}_{c \to s,t}$ over its control states, weighted by $n_{c,t} + n_{c,t-1}$ (or with each of them). If the two lines stay close to each other before treatment, integrated parallel trends is plausible for that state. In the first treated period the gap between them is the estimate of `didint()` for that state and period, $\widehat{ATT}(s, g_s)$.

The `"integrated2"` data measure the same changes from the pre-treatment period:

```math
\widehat{\delta}^{\,r}_{j \to s,t} = \left(\hat{\lambda}_{j,t} - \hat{\lambda}_{j,r_s}\right) + \left(\hat{\gamma}_{j,t} - \hat{\gamma}_{j,r_s}\right)'\,\overline{X}_{s,t}
```

For $t \geq g_s$ these are the [Integrated Long Differences](@ref) of `didint()` themselves; the data have them for the earlier periods as well, and they are zero in period $r_s$. The controls are averaged with the weights $n_{c,t} + n_{c,r_s}$. As in `"conditional2"`, the gap between the two lines is $\widehat{ATT}(s,t)$ in every treated period and a placebo estimate before treatment.

### Placebo Estimates Before Treatment

In the plots of changes since the pre-treatment period (`"conditional2"` and `"integrated2"`), the gap between the treated state and the average of its controls in a period before treatment is the same difference-in-differences that `didint()` computes for a treated period, applied where no treatment has taken place. Under parallel trends, and without anticipation of the treatment, it is zero in expectation. These gaps are the counterpart of the pre-treatment estimates of Callaway & Sant'Anna (2021) with a fixed base period, and the gaps before treatment in the `"integrated"` plot are the counterpart of their estimates with a base period that moves with $t$.

- Gaps before treatment that scatter around zero and are small next to the gaps after treatment support parallel trends.
- Gaps that drift steadily towards or away from zero as treatment approaches suggest that the treated state was on a different trend, which the gaps after treatment then partly continue.
- A jump in the last periods before treatment suggests anticipation.

The gaps before and after treatment are in the same units, so they show how large a violation would be relative to the estimated effect. A trend that differs between a treated state and its controls accumulates in the changes since the pre-treatment period, where the changes from one period to the next show it as a small gap in each period; in return, every point of a line shares the sampling error of period $r_s$. The gap is zero by construction in period $r_s$, and the plots do not show the sampling error of the placebo estimates.

### Event Study

For the `"event"` data, the $\hat{\lambda}_{s,t}$ values of the treated states are averaged by period relative to treatment, $e = t - g_s$, for each requested `ccc` option:

```math
\overline{\lambda}_{e} = \sum_{s} \omega_{s,e}\, \hat{\lambda}_{s,\,g_s + e}, \qquad \omega_{s,e} = \frac{n_{s,\,g_s+e}}{\sum_{s'} n_{s',\,g_{s'}+e}}
```

where the sums run over the treated states observed $e$ periods from treatment (`weights = false` uses equal weights). The standard error and confidence interval of each average come from the weighted regression of those $\hat{\lambda}_{s,t}$ values on a constant. The plot has one panel per `ccc` option and shows the average covariate-adjusted outcome of the treated states around treatment; intervals are drawn where at least `groupmin` (by default 3) treated states are behind the average. It contains no comparison states, so it describes the path of the treated states and does not estimate dynamic treatment effects, which `didint()` reports with `agg = "time"`.

### Drawing the Plots

`plot` draws a `DiDIntPlot` once a [Makie](https://docs.makie.org) backend such as CairoMakie is loaded. Its `plot_type` argument chooses among `"event"`, `"conditional"`, `"conditional2"`, `"integrated"` and `"integrated2"`, and `treated_states` selects the treated states to show. A dashed vertical line marks the first treated period, and the panels of treated states are ordered from the earliest treated state to the latest.

In the panels of a treated state, `weights = true` (the default) gives the control states the weights that `didint()` gives them under `weighting = "both"` or `"diff"`, and `weights = false` the equal weights of `"att"` or `"none"`.

## Fast Computation of Randomization Inference and Jackknife

The randomization inference procedure and the fast jackknife recompute every estimate many times: once for each alternative assignment of the treatment times, or once for each removed state. **DiDInt.jl** computes the estimates of all of these cases together. The $\hat{\lambda}_{s,t}$ values do not depend on which states are treated, so Step 1 is done once, and every later step is a weighted mean.

### Fast Randomization Inference

The alternative assignments, drawn as described in the [Randomization Inference](@ref) section, are stored in a matrix $L$ with one row per assignment and one column per state; each entry is the treatment time the assignment gives that state, or 0 if untreated. For each cell $(g, t)$ of a treatment time $g$ and a period $t \geq g$, with `d`, `v`, `u` and `n` holding one value per state:

```julia
# with notyet = true: controls = (L .== 0) .| (L .> t)
controls = (L .== 0)
treated  = (L .== g)

# control mean of every assignment, ATT(s,t) of every state in every assignment
a   = (controls * (v .* d)) ./ (controls * v)
att = d' .- a

numerator   .+= (treated .* att) * u
denominator .+= treated * u
nobs        .+= treated * n
```

Here $d_j = \hat{\lambda}_{j,t} - \hat{\lambda}_{j,r}$ is the long difference of state $j$ (with $r = g - 1$), $v_j$ its weight as a control, $u_j$ its weight as a treated state, and $n_j = n_{j,t}$ (see [Weighting](@ref)). Entry $(b, j)$ of `att` is $\widehat{\mathrm{ATT}}_b(j,t)$, the $\widehat{\mathrm{ATT}}(s,t)$ of state $j$ in assignment $b$. The three sums are kept separately for each subaggregate; after the last cell, `numerator ./ denominator` gives every subaggregate ATT of every assignment, and the aggregate ATT of each assignment is their mean, weighted by `nobs` under `"both"` and `"att"` and equally otherwise. States without a long difference in the cell are left out of both the controls and the treated states, and an assignment without any control for a cell has no $\widehat{\mathrm{ATT}}(s,t)$ in that cell.

Under `ccc = "int"` with covariates, each long difference also includes its covariate term, $(\hat{\gamma}_{j,t} - \hat{\gamma}_{j,r})'\,\overline{X}_{i,t}$, evaluated at the covariate means of the state $i$ treated in that case, using the same first-stage estimates in every case.

### Fast Jackknife

The fast jackknife uses a similar computation with one row per removed state. Every row of $L$ is the actual assignment, and row $h$ leaves state $h$ out of both the controls and the treated states:

```julia
L    = repeat(permutedims(actual), S)   # S rows, one per removed state
keep = .!Matrix{Bool}(I, S, S)          # keep[h, j] is false when j == h
controls = (L .== 0) .& keep
treated  = (L .== g) .& keep
```

The rest is unchanged, and row $h$ of the result is the estimate without state $h$.
