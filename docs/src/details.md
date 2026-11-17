# Under the Hood

Most pertinent details of **DiDInt.jl** can be found in the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447), however, some details that might be of particular interest to users of the **DiDInt.jl** package are elaborated upon here.

## Four Step Estimation Procedure

The **DID-INT** estimator is computed in four steps:

1. The means of the outcome of interest (conditional on covariates) are computed for each period for each state (denoted as $\lambda$ in the DID-INT paper).
2. Within state differences between these $\lambda$ values (long differences) are computed for treated states and for their controls.
3. An ATT is computed for each treated state $s$ in each period $t$ from its treatment time on, denoted $\mathrm{ATT}(s,t)$: the long difference of the treated state minus the (weighted) mean of the long differences of its controls.
4. The (weighted) means of the $\mathrm{ATT}(s,t)$ values within groups give the subaggregate ATTs (more details on the choice of groups in the [Aggregation](@ref) section), and the (weighted) mean of the subaggregate ATTs gives the aggregate ATT.

More explicitly:

**Step 1.** In principle, the following regression is estimated to obtain conditional outcome means $\hat{\lambda}_{s,t}$:

```math
Y_{i,s,t} = \sum_s \sum_t \lambda_{s,t} I(s,t) + f\!\left(X^k_{i,s,t}\right) + \epsilon_{i,s,t}
```

where $I(s,t)$ is a dummy variable equal to 1 if observation $i$ belongs to state $s$ in period $t$, and $f\!\left(X^k_{i,s,t}\right)$ is a function of covariates whose form depends on the common causal covariates (CCC) violation that is being accounted for (see [Common Causal Covariates and Model Specification](@ref)). The estimated $\hat{\lambda}_{s,t}$ values are the covariate-adjusted means of the outcome for state $s$ in period $t$.

!!! note "Implementation"
    In practice, this regression is not estimated directly on the full dataset as this would be computationally prohibitive for datasets with many states, periods, and covariates. Rather, the Frisch-Waugh-Lovell theorem is applied in the case of homogeneous DID-INT, region-varying DID-INT, time-varying DID-INT, and two one-way DID-INT to recover covariate coefficients which can then be used in conjunction with the outcome and covariate means per $(s,t)$ group to retrieve the $\hat{\lambda}_{s,t}$ values. For two-way DID-INT, the regression is run separately for each $(s,t)$ group. The resulting $\hat{\lambda}_{s,t}$ values are numerically equivalent to those from the combined regression shown above. More details on estimating $\hat{\lambda}_{s,t}$ can be found in the [Estimating Lambda](@ref) section.

**Step 2.** For each treated state $s$ belonging to treatment cohort $g$ (the states first treated in period $g$), long differences are computed between each period $t \geq g$ and the period immediately prior to treatment, $r = g - 1$:

```math
\widehat{\Delta}_{s,t} = \left(\hat{\lambda}_{s,t} - \hat{\lambda}_{s,r}\right) + \left(\hat{\gamma}_{s,t} - \hat{\gamma}_{s,r}\right)'\,\overline{X}_{s,t}
```

The analogous difference is computed for each control state $c$ of the treated state $s$:

```math
\widehat{\Delta}_{c \to s,t} = \left(\hat{\lambda}_{c,t} - \hat{\lambda}_{c,r}\right) + \left(\hat{\gamma}_{c,t} - \hat{\gamma}_{c,r}\right)'\,\overline{X}_{s,t}
```

Here $\hat{\gamma}_{s,t}$ are the covariate coefficients that apply to state $s$ in period $t$ under the chosen CCC violation, and $\overline{X}_{s,t}$ are the means of the covariates of the *treated* state $s$ in period $t$. These are the integrated long differences of the DID-INT paper (equations 46 and 47): the change in the covariate-adjusted outcome is evaluated over the covariate distribution of the treated state, which, because $f$ is linear in the covariates, only requires that state's covariate means. More on when the covariate term matters in the [Integrated Long Differences](@ref) section.

**Step 3.** For each treated state $s$ and each period $t$ from its treatment time on, the ATT is the long difference of $s$ minus the weighted mean of the long differences of its controls (equation 48 of the DID-INT paper):

```math
\widehat{\mathrm{ATT}}(s,t) = \widehat{\Delta}_{s,t} - \frac{\sum_{c \in C(s,t)} v_{c}\, \widehat{\Delta}_{c \to s,t}}{\sum_{c \in C(s,t)} v_{c}}
```

where $C(s,t)$ is the set of controls: the never-treated states and, optionally, the states not yet treated in period $t$ (see [Not Yet Treated Cells](@ref)). In this step, weights $v_c$ are used (by default) according to the number of observations associated with each control long difference (see [Weighting](@ref)).

!!! note "Implementation"
    $\widehat{\mathrm{ATT}}(s,t)$ is computed as the coefficient $\beta$ in the weighted regression of the long differences of state $s$ and of its controls on an intercept and a treatment indicator, $\widehat{\Delta} = \alpha + \beta\, d + \varepsilon$, where $d = 1$ for the treated state and 0 for the controls. This coefficient equals the expression above and the regression also yields its standard error.

**Step 4.** Each subaggregate ATT is a weighted mean of the $\widehat{\mathrm{ATT}}(s,t)$ values that belong to its group $j$, and the aggregate ATT is a weighted mean of the subaggregate ATT estimates:

```math
\widehat{\mathrm{ATT}}_{j} = \frac{\sum_{(s,t) \in j} u_{s,t}\, \widehat{\mathrm{ATT}}(s,t)}{\sum_{(s,t) \in j} u_{s,t}}, \qquad
\widehat{\mathrm{ATT}} = \sum^{J}_{j = 1} w_{j}\, \widehat{\mathrm{ATT}}_{j}
```

where $J$ denotes the total number of subaggregate ATTs and the weights, $w_j$, are scaled such that $\sum^{J}_{j = 1} w_{j}\, = 1$. By default the weights used here reflect the number of treated observations associated with each $\widehat{\mathrm{ATT}}(s,t)$ and with each subaggregate ATT (see [Weighting](@ref)).

!!! note "Implementation"
    Each subaggregate ATT is computed as the coefficient $\beta_j$ in the weighted regression of the $\widehat{\mathrm{ATT}}(s,t)$ values of group $j$ on a constant, $\widehat{\mathrm{ATT}}(s,t) = \beta_j + \varepsilon_{s,t}$, with weights $u_{s,t}$. This coefficient equals the weighted mean above, and the regression also yields the standard error of the subaggregate ATT: the heteroskedasticity-consistent standard error of $\beta_j$, which reflects how much the $\widehat{\mathrm{ATT}}(s,t)$ values of the group vary around their weighted mean. It does not use the standard errors of the individual $\widehat{\mathrm{ATT}}(s,t)$ values from Step 3.

    The aggregate ATT is computed in the same way, as the coefficient in the weighted regression of the subaggregate ATTs on a constant with weights $w_j$, which also yields its standard error. A subaggregate ATT that consists of a single $\widehat{\mathrm{ATT}}(s,t)$ has no variation to compute a standard error from, and reports the standard error of that $\widehat{\mathrm{ATT}}(s,t)$ from Step 3 instead; likewise an aggregate ATT with a single subaggregate ATT reports the standard error of that subaggregate ATT. See the [Standard Errors and P-Values](@ref) section.

No regression pools long differences from different treatment cohorts or periods: every estimate is an explicit weighted mean of $\widehat{\mathrm{ATT}}(s,t)$ values, each of which only compares states within the same period.

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

The covariate term in the long differences of Step 2 of the [Four Step Estimation Procedure](@ref), $\left(\hat{\gamma}_{\cdot,t} - \hat{\gamma}_{\cdot,r}\right)'\,\overline{X}_{s,t}$, depends on the DID-INT variation:

- **Homogeneous** and **region-varying** DID-INT: the covariate coefficients do not change over time, so the term is zero in every long difference.
- **Time-varying** and **two one-way** DID-INT: the term is not zero, but it is the same for the treated state and for each of its controls (they share the time-specific coefficients, and both are evaluated at $\overline{X}_{s,t}$), so it cancels in $\widehat{\mathrm{ATT}}(s,t)$.
- **Two-way** DID-INT: the coefficients are specific to each state and period, and the term does not cancel.

Consequently the covariate term changes the estimates only when `ccc = "int"` is used with covariates. For every other `ccc` option, and whenever no covariates are specified, $\widehat{\mathrm{ATT}}(s,t)$ only depends on the $\hat{\lambda}_{s,t}$ values:

```math
\widehat{\mathrm{ATT}}(s,t) = \left(\hat{\lambda}_{s,t} - \hat{\lambda}_{s,r}\right) - \frac{\sum_{c \in C(s,t)} v_{c} \left(\hat{\lambda}_{c,t} - \hat{\lambda}_{c,r}\right)}{\sum_{c \in C(s,t)} v_{c}}
```

## Aggregation

In Step 4 of the [Four Step Estimation Procedure](@ref), a choice needs to be made in terms of which $\widehat{\mathrm{ATT}}(s,t)$ values should be averaged together. For example, we could average the $\widehat{\mathrm{ATT}}(s,t)$ values of the states treated at treatment time $g$, which would yield the subaggregate ATT for $g$. Repeating this for each treatment time would give us the set of subaggregate ATTs based on `"cohort"` aggregation. Conversely, we could average the $\widehat{\mathrm{ATT}}(s,t)$ values of a particular treated state, $s^*$, in order to retrieve the subaggregate ATTs per treated state. The following table gives an overview of the aggregation groupings and the string value passed to the `agg` argument in the **DiDInt.jl** package. Every option is available for both staggered adoption and common adoption scenarios.

| Grouping                                                         | $\widehat{\mathrm{ATT}}(s,t)$ values averaged in each subaggregate ATT        | `agg`                |
| :--------------------------------------------------------------- | :---------------------------------------------------------------------------- | :------------------- |
| By treatment time, $g$                                           | Those of the states treated at $g$, in every period                           | `"cohort"` (default) |
| By treated state, $s^*$                                          | Those of state $s^*$, in every period                                         | `"state"`            |
| By treatment time, $g$, and period, $t$                          | Those of the states treated at $g$, in period $t$                             | `"simple"`           |
| By treated state, $s^*$, and period, $t$                         | None: each $\widehat{\mathrm{ATT}}(s,t)$ is reported                          | `"sgt"`              |
| By periods since treatment, $p$                                  | Those with $t - g = p$, across treatment times                                | `"time"`             |
| By period, $t$                                                   | Those in period $t$, across treatment times                                   | `"calendar"`         |
| By policy (and within policy)                                    | See the [Policies](@ref) section                                              | `"policy"`, `"policy_gt"`, `"policy_cohort"`, `"policy_calendar"`, `"policy_time"` |

Whichever option is chosen, the aggregate ATT (`agg_att`) is always reported: it is the weighted mean of the subaggregate ATTs of that option. Because it is a mean of subaggregates, the aggregate ATT can differ across `agg` options. It is the same for every `agg` option when `weighting` is `"att"` or `"both"` and `treated_weight = "treated"` (the defaults), in which case it is the mean of the $\widehat{\mathrm{ATT}}(s,t)$ values weighted by their numbers of treated observations (see the example [The Aggregate ATT Across Aggregation Options](@ref)).

## Weighting

An estimate in **DiDInt.jl** is built from three successive averages, and each of them can be weighted by numbers of observations:

1. the long differences of the **control states** are averaged within each $\widehat{\mathrm{ATT}}(s,t)$ (Step 3 of the [Four Step Estimation Procedure](@ref));
2. the **$\widehat{\mathrm{ATT}}(s,t)$ values** are averaged within each subaggregate ATT (Step 4);
3. the **subaggregate ATTs** are averaged within the aggregate ATT (Step 4).

The `weighting` argument decides which of these averages are weighted, and the `treated_weight` argument decides which observations of a treated state count in the second average. By default (`weighting = "both"`, `treated_weight = "treated"`) all three averages are weighted. The table gives an overview, and the three averages are described in turn below.

| `weighting`        | 1. Controls within an $\widehat{\mathrm{ATT}}(s,t)$ | 2. $\widehat{\mathrm{ATT}}(s,t)$ values within a subaggregate ATT        | 3. Subaggregate ATTs within the aggregate ATT |
| :----------------- | :-------------------------------------------------- | :----------------------------------------------------------------------- | :-------------------------------------------- |
| `"both"` (default) | Weighted                                            | Weighted                                                                 | Weighted                                      |
| `"att"`            | Equal                                               | Weighted with `treated_weight = "treated"`, equal with `"population"`    | Weighted                                      |
| `"diff"`           | Weighted                                            | Equal with `treated_weight = "treated"`, weighted with `"population"`    | Equal                                         |
| `"none"`           | Equal                                               | Equal                                                                    | Equal                                         |

Throughout, $n_{s,t}$ denotes the number of observations of state $s$ in period $t$, and $r$ the period immediately before the treatment time of the treated state in question.

### First Average: Control States Within an ATT(s,t)

Each $\widehat{\mathrm{ATT}}(s,t)$ compares one treated state with its control states $C(s,t)$:

```math
\widehat{\mathrm{ATT}}(s,t) = \widehat{\Delta}_{s,t} - \frac{\sum_{c \in C(s,t)} v_{c}\, \widehat{\Delta}_{c \to s,t}}{\sum_{c \in C(s,t)} v_{c}}
```

with the weight of control state $c$ given by:

```math
v_c = \begin{cases} n_{c,t} + n_{c,r} & \text{if } \texttt{weighting} \text{ is } \texttt{"diff"} \text{ or } \texttt{"both"} \\ 1 & \text{otherwise} \end{cases}
```

The weights only apply to the control states: there is a single treated state, so nothing is averaged on the treated side. A weighted control state counts in proportion to the number of observations behind its long difference, which are its observations in period $t$ and in period $r$.

### Second Average: ATT(s,t) Values Within a Subaggregate ATT

Each subaggregate ATT averages the $\widehat{\mathrm{ATT}}(s,t)$ values of its group $j$ (see [Aggregation](@ref)). These are estimates for treated states, and the weights are numbers of observations of those treated states:

```math
\widehat{\mathrm{ATT}}_{j} = \frac{\sum_{(s,t) \in j} u_{s,t}\, \widehat{\mathrm{ATT}}(s,t)}{\sum_{(s,t) \in j} u_{s,t}}
```

with the weight $u_{s,t}$ given by:

| `weighting` | $u_{s,t}$ with `treated_weight = "treated"` (default) | $u_{s,t}$ with `treated_weight = "population"` |
| :---------- | :----------------------------------------------------- | :---------------------------------------------- |
| `"both"`    | $n_{s,t}$                                              | $n_{s,t} + n_{s,r}$                             |
| `"att"`     | $n_{s,t}$                                              | 1                                               |
| `"diff"`    | 1                                                      | $n_{s,t} + n_{s,r}$                             |
| `"none"`    | 1                                                      | 1                                               |

With `treated_weight = "treated"`, an $\widehat{\mathrm{ATT}}(s,t)$ is weighted by the treated state's observations in period $t$, so that a subaggregate ATT estimates the effect on the average treated observation of its group. With `treated_weight = "population"`, it is weighted like a control long difference: by the observations of both the period $t$ and the period before treatment that were used in its calculation. The two options give the same subaggregate ATT when it averages a single $\widehat{\mathrm{ATT}}(s,t)$, and can otherwise differ, including for a single treated state observed over several periods (the ratio $n_{s,r}/n_{s,t}$ changes from period to period).

Under `"sgt"` aggregation every subaggregate ATT is a single $\widehat{\mathrm{ATT}}(s,t)$, so this average does not occur.

### Third Average: Subaggregate ATTs Within the Aggregate ATT

The aggregate ATT averages the $J$ subaggregate ATTs:

```math
\widehat{\mathrm{ATT}} = \sum_{j=1}^{J} w_j\, \widehat{\mathrm{ATT}}_j
```

with the weight of subaggregate ATT $j$ given by:

```math
w_j = \begin{cases} \dfrac{N_j}{\sum_{l=1}^{J} N_l} & \text{if } \texttt{weighting} \text{ is } \texttt{"att"} \text{ or } \texttt{"both"} \\[10pt] \dfrac{1}{J} & \text{otherwise} \end{cases}
\qquad \text{where} \qquad
N_j = \sum_{(s,t) \in j} n_{s,t}
```

A weighted subaggregate ATT counts in proportion to $N_j$, its number of *treated* observations: the observations of the treated states in the periods of the $\widehat{\mathrm{ATT}}(s,t)$ values that it averages. This is the same whichever `treated_weight` option is used. The weights $w_j$ are reported alongside the subaggregate ATTs (the `weight` column of the `sub` table of the results).

For the policy aggregation options, the ATT of each policy is formed from the subaggregate ATTs of that policy in the same way, and the pooled ATT across policies has its own weights (see [Policies](@ref)).

Under the default options the second average uses $u_{s,t} = n_{s,t}$ and the third uses $N_j$, the sum of those same weights. The two averages then collapse into one:

```math
\widehat{\mathrm{ATT}} = \frac{\sum_{j} N_j\, \widehat{\mathrm{ATT}}_j}{\sum_{j} N_j} = \frac{\sum_{(s,t)} n_{s,t}\, \widehat{\mathrm{ATT}}(s,t)}{\sum_{(s,t)} n_{s,t}}
```

which does not depend on how the $\widehat{\mathrm{ATT}}(s,t)$ values were grouped. This is why the aggregate ATT is the same for every `agg` option when `weighting` is `"att"` or `"both"` and `treated_weight = "treated"`.

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

For these options the ATT of each policy, $\widehat{\mathrm{ATT}}_p$, is reported in the `att_p` column (on the first row of each policy) along with its standard error and p-values. It is the weighted mean of the subaggregate ATTs of that policy, in the same way that the aggregate ATT is the weighted mean of all subaggregate ATTs for the other aggregation options. There is no `agg_att` column for the policy aggregation options.

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

1. **Basic matching**: An observation at time $t$ is matched to the latest period grid date that is less than or equal to $t$.
2. **Treatment boundary adjustment**: If there exists a treatment time $T^*$ between two adjacent periods on the grid, observations are matched as follows:

   - If $t < T^*$: the observation is matched to the earlier of the two periods
   - If $t \geq T^*$: the observation is matched to the later of the two periods

   This adjustment ensures that observations are not incorrectly assigned to pre-treatment periods when they occur after treatment begins within that period.
3. **Treatment time matching**: Treatment times themselves are matched to the first grid period date that is greater than or equal to the treatment time, ensuring treatments are associated with the correct period on the grid.

### Treatment Times Outside the Data

A treated state is part of the analysis as a treated state only if its treatment time falls within the data used, that is, after the first date and no later than the last date (after any restriction through `start_date` and `end_date`). Otherwise, with a warning:

- A state treated **at or before the first date** is treated throughout the data. It has no pre-treatment period, so no ATT can be computed for it, and it is not a valid control either. Its observations are dropped by `didint()`.
- A state treated **after the last date** is untreated throughout the data, and is used as a control like a never-treated state.

Each treated state can only have one treatment time: listing a state twice, or a `gvar` column that takes more than one value within a state, is an error.


## Standard Errors and P-Values

**DiDInt.jl** reports three types of p-values and two sets of standard errors, for every subaggregate ATT and for the aggregate ATT (and, for the policy aggregation options, for the ATT of each policy and for the pooled ATT).

- Standard errors using heteroskedasticity-consistent covariance estimators (HCCME): HC0, HC1 (default), HC2, HC3, or HC4 and p-values via a two-tailed test using a $t$-distribution with $n - k$ degrees of freedom, where $n$ is the number of observations in the relevant regression and $k$ is the number of estimated regression parameters. The desired HCCME can be selected via the `hc` argument.
- Jackknife standard errors and p-values, described in more detail in the [Jackknife](@ref) section.
- Randomization inference (RI) p-values, via a randomization procedure described in the [Randomization Inference](@ref) section.

There are some circumstances in which the jackknife standard errors and p-values or the randomization inference p-values cannot be calculated. These situations are noted in the [Jackknife](@ref) and [Randomization Inference](@ref) sections, respectively.

The standard errors using HCCMEs come from the following regressions:

- **$\widehat{\mathrm{ATT}}(s,t)$**: the Step 3 regression of the long differences of the treated state and of its controls on an intercept and a treatment indicator ($n$ is one plus the number of controls and $k = 2$).
- **Subaggregate ATTs**: the weighted regression of the $\widehat{\mathrm{ATT}}(s,t)$ values of the group on a constant, using the weights $u_{s,t}$ ($n$ is the number of $\widehat{\mathrm{ATT}}(s,t)$ values in the group and $k = 1$). A subaggregate ATT that consists of a single $\widehat{\mathrm{ATT}}(s,t)$ reports the standard error and p-value of that $\widehat{\mathrm{ATT}}(s,t)$.
- **Aggregate ATT**: the weighted regression of the subaggregate ATTs on a constant, using the weights $w_j$ ($n$ is the number of subaggregate ATTs and $k = 1$). If there is a single subaggregate ATT, the aggregate ATT reports its standard error and p-value.

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

There are two variations to the **DiDInt.jl** jackknife procedure. The default option is to calculate jackknife estimates of the ATTs while dropping one state at a time directly from Step 3 in the [Four Step Estimation Procedure](@ref), and the $\widehat{\mathrm{ATT}}(s,t)$ values, the subaggregate ATTs, and the aggregate ATT are recomputed from the $\hat{\lambda}_{s,t}$ values of the remaining states. This is referred to as the "fast jackknife". The alternative method, the "true jackknife", calculates jackknife estimates of the ATTs while dropping one state at a time starting from Step 1. The true jackknife can be used in lieu of the fast jackknife by setting the `truejack` argument to `true`.

Note that for the `ccc` options of `"state"` and `"int"` that the true jackknife is always equal to the fast jackknife as the covariate effects are local to the state or state-time, respectively, so that removing a state from the dataset does not affect the estimated covariate effects at other states and thus does not affect the estimated $\hat{\lambda}_{s,t}$ values at other states.

The jackknife standard error is computed as:

```math
\widehat{\mathrm{se}}_{\mathrm{jack}} = \sqrt{\frac{N-1}{N} \sum_{i=1}^{N}\left(\hat{\beta}_{(-i)} - \hat{\beta}\right)^2}
```

where $N$ is the number of states whose removal leaves the estimate defined, $\hat{\beta}_{(-i)}$ is the ATT estimate with state $i$ removed, and $\hat{\beta}$ is the original ATT estimate. This is the cluster jackknife variance estimator centered at the original estimate ($\mathrm{CV}_3$ in [MacKinnon, Nielsen, and Webb (2023)](https://doi.org/10.1016/j.jeconom.2022.04.001)). Jackknife p-values are computed via a two-tailed test using a $t$-distribution with $N - 1$ degrees of freedom.

A minimum of two treated states and two control states are required to compute jackknife standard errors and p-values; if this condition is not met, these quantities are returned as `missing`. The same minimum applies to each estimate separately: a subaggregate ATT (or the ATT of a policy) only has a jackknife standard error if the $\widehat{\mathrm{ATT}}(s,t)$ values it averages involve at least two treated states and at least two control states, and if it can be computed with any one state removed. Otherwise its jackknife standard error and p-value are returned as `missing`.

The states that are removed in turn are those that contribute at least one long difference to the estimates; a state that never does (for example one that is only observed before every treatment time) is not counted. If removing a state leaves a subaggregate ATT without any treated state (for example the only state of a treatment cohort), that subaggregate ATT has no jackknife standard error (it is `missing`), and the aggregate ATT for that removal is the mean over the remaining subaggregate ATTs. The jackknife standard error of the aggregate ATT is still reported. The same applies to the pooled ATT across policies when removing a state leaves a policy without any treated state.

Further, note that for the aggregation methods of `"state"` and `"sgt"`, where each subaggregate ATT involves only one treated state, that the jackknife standard errors and jackknife p-values are never estimated at the subaggregate level: removing that one treated state leaves nothing to estimate. This limitation does **not** prevent the jackknife standard error and jackknife p-value from being calculated for the *aggregate* ATT under `"state"` or `"sgt"` aggregation, so long as there are at least two treated and two control states.

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
    The method for drawing alternative assignments is essentially a combination of two classic procedures: Floyd's algorithm for sampling distinct integers, and the problem of finding the $k$-th permutation of a lexicographically ordered set. The latter is used in its multiset version, which is adjusted for values that occur more than once in the set, since several states can share a treatment time and several states are never treated.

    **Assignments as numbered arrangements.** With the states in a fixed order, an assignment is an arrangement of labels across the states, where the label of a state is its treatment time or "never treated" (for the policy aggregation options, its combination of treatment time and policy). Every possible assignment is a rearrangement of the labels of the actual assignment. Swapping two states that carry the same label changes nothing, which is why the number of distinct arrangements is the multinomial coefficient in step 1 rather than $N!$. Listing the distinct arrangements in lexicographic (dictionary) order gives each of them a number: $0$ for the first, $1$ for the second, and so on.

    **From a number to its arrangement.** The arrangement with a given number $k$ can be built directly, one state at a time, without listing any other arrangement. For the first state, the candidate labels are tried in order. The number of arrangements that begin with the smallest label is again a multinomial coefficient (that of the labels left over); if $k$ is smaller than this count, the first state receives that label, and otherwise the count is subtracted from $k$ and the next label is tried. The same is then done for the second state with the labels that remain, and so on until every state has a label.

    **Drawing distinct numbers.** Drawing alternative assignments is thereby reduced to drawing distinct whole numbers. Floyd's algorithm ([Bentley and Floyd, 1987](https://doi.org/10.1145/30401.315746)) draws $nperm^{*}$ distinct numbers from a range such that every set of numbers is equally likely, using exactly $nperm^{*}$ random draws: no draw ever has to be discarded for repeating an earlier one. Each drawn number is then converted into its arrangement. This guarantees distinct, equally likely assignments, however many states and treatment times there are. (The number of arrangements can exceed what a 64-bit integer holds; larger integer types are used when needed.)

    **Leaving out the actual assignment.** The actual assignment is itself one of the arrangements, with some number $k_0$, and must not be drawn. If there are $A$ alternative assignments, there are $A + 1$ arrangements in total, numbered $0$ to $A$. The numbers are drawn from the shorter range $0$ to $A - 1$, and every drawn number that is greater than or equal to $k_0$ is then increased by one. The drawn numbers below $k_0$ stay as they are and those from $k_0$ on move up by one, so together they can land on every number from $0$ to $A$ except $k_0$, each with the same probability. When every alternative assignment is used (step 2), the numbers $0$ to $A - 1$ are simply all taken, with the same shift.

### Kernel-Smoothed P-Values

With few alternative assignments, the p-value above can only take a few distinct values. When the number of permutations used, $nperm^{*}$, is no larger than the `ri_kernel_cutoff` argument (99 by default), **DiDInt.jl** instead reports the kernel-smoothed p-value of [Racine and MacKinnon (2007)](https://doi.org/10.1016/j.csda.2006.11.013), which replaces the indicator function with a smooth function:

```math
p_{RI} = 1 - \frac{1}{nperm^{*}} \sum_{r=1}^{nperm^{*}} \Phi\!\left(\frac{|ATT_{Actual}| - |ATT_{RI,r}|}{h}\right), \qquad h = 1.575\; \hat{\sigma}_{RI}\; (nperm^{*})^{-4/9}
```

where $\Phi$ is the standard normal cumulative distribution function and $\hat{\sigma}_{RI}$ is the sample standard deviation of the $ATT_{RI,r}$ values. This applies whether $nperm^{*}$ is small because few alternative assignments exist or because a small `nperm` was requested. Kernel smoothing can be turned off by setting `ri_kernel_cutoff` to `nothing`, `false`, or `0`. If there is only one alternative assignment, or if the $ATT_{RI,r}$ values do not vary, the unsmoothed p-value is reported.

### Ties

Although every alternative assignment differs from the actual assignment, an estimate under an alternative assignment can be exactly equal to the actual estimate. Assignments are unique as a whole, but a given subaggregate ATT only depends on part of the assignment: for example, the subaggregate ATT for treatment time $g$ and period $t$ only depends on which states are treated at $g$ and which states are controls in period $t$. An alternative assignment that reshuffles other states leaves that subaggregate ATT unchanged. Such ties are most common when `notyet = true` for early treatment cohorts with few states, as the set of not-yet-treated controls in an early period changes little across assignments. The aggregate ATT depends on the entire assignment and almost never ties.

Ties are counted in favour of the null hypothesis (they satisfy $|ATT_{RI,r}| \geq |ATT_{Actual}|$). Because two estimates that are equal in exact arithmetic can differ by rounding error, estimates within a small relative tolerance (about $1.5 \times 10^{-8}$) of the actual estimate are treated as ties.

This applies to the unsmoothed p-value. With [Kernel-Smoothed P-Values](@ref), an estimate equal to the actual estimate contributes $\Phi(0) = 0.5$ to the sum, so a tie counts as half of an exceedance rather than a whole one.

### Aggregation by State

When the randomization inference procedure is performed while the aggregation is set to `"state"` or `"sgt"`, each subaggregate ATT involves a single treated state, and an alternative assignment does not determine which state should take its place when more than one state receives the same treatment time. There is therefore an additional level of randomization: for each actually treated state, one state is drawn at random from the set of states that the alternative assignment associates with the treatment time of the actually treated state, and is used in its place in every period. This is in order to preserve the structure of the `"state"` and `"sgt"` aggregation methods whereby only one treated state is considered at a time.

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

If, for any $(s,t)$ group, a covariate causes the design matrix to be rank-deficient, then that covariate is dropped from the regression for that $(s,t)$ group, and its coefficient is taken to be zero in that group when computing the [Integrated Long Differences](@ref). A covariate that does not vary within the group is dropped, as it cannot be told apart from the intercept. Among the covariates that do vary, the ones to drop are found by going through the covariates in the order in which they were given in the `covariates` argument (with the dummy variables of a categorical covariate in its place): a covariate is dropped if, after centering at the group means, it is a linear combination of the covariates kept before it. The covariates are not reordered in any way (no column pivoting is involved), so which of several collinear covariates is dropped only depends on the order of `covariates`, and not on the units in which the covariates are measured.

Centering matters because the regression has an intercept: covariates can be collinear with the intercept without being collinear with each other. The typical case is a categorical covariate whose reference category does not occur in the group, so that the group's dummy variables add up to one; the dummy variable listed last is then dropped in that group. The group contains no information about the absent category, and the estimates can then depend on which category was chosen as the reference. A group with $n$ observations can keep at most $n - 1$ covariates, so in very small groups the covariates listed last are dropped.

!!! note "Collinear Covariates"
    Dropping collinear covariates is a choice: when covariates are collinear within an $(s,t)$ group, the data of that group cannot tell their coefficients apart, and a different choice (dropping another one of the collinear covariates, or spreading the coefficient across them) fits the group's data equally well but can yield different ATT estimates. **DiDInt.jl** resolves this by a fixed rule that depends only on the order of the `covariates` argument, and issues a warning that names every covariate dropped for being collinear and the groups in which it was dropped. The same rule applies under `"hom"`, `"state"`, and `"time"`, in the part of the data over which the covariate coefficients are estimated; see [Rank Deficiency and Lambda Identification](@ref). No warning is given for covariates that are dropped only because they do not vary within a group.

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

When the covariates (demeaned within cells) are linearly dependent, the $\hat{\lambda}_{s,t}$ values are not uniquely determined by the data. The ATT is unaffected only when the resulting ambiguity in $\hat{\lambda}_{s,t}$ is additive in state and period, because such shifts cancel in the long differences; collinearity among covariates in general does not have this form.

**Homogeneous, region-varying, time-varying, and two-way DID-INT.** Under `"hom"`, `"state"`, `"time"` and `"int"`, **DiDInt.jl** therefore drops redundant covariates within each block of the data in which the covariate coefficients are estimated (all data, each state, each period, or each $(s,t)$ cell). It goes through the covariates in the order given in `covariates` and drops any covariate that is, up to a relative tolerance, a linear combination of those already kept: a covariate is dropped when less than $10^{-7}$ of it (measured by its norm, after demeaning within cells) is left once the part explained by the covariates kept before it is removed. The estimates are then unique, unchanged when covariates are rescaled, and identical across procedures. When a covariate that varies within cells is dropped, a warning names it: the estimates then depend on which covariate was dropped, which can be changed by reordering `covariates`. Dropping a covariate does not make the data informative about its coefficient; it fixes a convention. A covariate that does not vary within any cell of a block is absorbed by the $I(s,t)$ cell dummies and is dropped without a warning; a covariate measured at the state-period level is the typical case, and it has no effect on the estimates under any DID-INT variation.

**Two one-way DID-INT.** Under `"add"`, the design is always rank deficient: adding a constant to every state's slopes and subtracting it from every period's slopes leaves each cell's slope $\gamma_s + \gamma_t$, and therefore every $\hat{\lambda}_{s,t}$, unchanged. Any least-squares solution gives the same $\hat{\lambda}_{s,t}$ values and ATTs, so **DiDInt.jl** solves this regression with Julia's `\` operator and does not drop any covariate. It checks that this is the only rank deficiency and warns otherwise, for example when a covariate does not vary within any cell of some state or period, or when covariates are collinear within cells; the ATTs may then depend on the solution chosen. Edge case standard errors are not available under `"add"`.

The check compares the rank of the design matrix $Z$ with the rank that the split between states and periods leaves. With $K$ covariates that vary within cells, $n_s$ states and $n_t$ periods, $Z$ has $K (n_s + n_t)$ columns: one coefficient per covariate for each state and for each period. The split costs one rank per covariate, so the highest rank that $Z$ can have is

```math
K \left(n_s + n_t - C\right)
```

where $C$ is the number of separate groups that the panel falls into, a group being a set of states and periods that are linked by observations ($C = 1$ unless some states are only observed in periods in which none of the other states are observed). Full column rank is never reached, and the warning is given only when the rank is below this value. For example, with 20 states, 12 periods and 3 covariates, $Z$ has $3 \times 32 = 96$ columns and the expected rank is $3 \times (20 + 12 - 1) = 93$. If one of the covariates does not vary within any cell of one particular state, that state's coefficient for that covariate is not determined by the data, the rank is 92, and the warning reads: *"Under ccc = "add" the covariate slopes are less well determined than usual: the first-stage regression has rank 92 where 93 is expected. ..."*

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

- **Two-way** (`"int"`): a separate regression is run for each $(s,t)$ cell, with the cell dummy and covariates as regressors. The HCCME for each cell's dummy coefficient is placed into the corresponding diagonal entry of $\widehat{\mathrm{Var}}(\hat{\lambda})$. Off-diagonal entries between different cells are zero since the regressions are run separately, so only the four variance terms of the decomposition remain.
- **Region-varying** (`"state"`): a separate regression is run for each state $s$, pooling all time periods within that state. The HCCME for the cell dummies within each state block is placed into the corresponding entries of $\widehat{\mathrm{Var}}(\hat{\lambda})$ (including off-diagonal covariance terms). The covariances between the two periods of the same state are nonzero; those between the treated state and the control state are zero.
- **Time-varying** (`"time"`): analogous to region-varying, with separate regressions run for each time period $t$. The covariances between the treated state and the control state in the same period are nonzero; all others are zero.
- **Homogeneous** (`"hom"`): a single global regression is run across all observations. The HCCME for all cell dummies from the regression is exactly $\widehat{\mathrm{Var}}(\hat{\lambda})$, and every term of the decomposition can be nonzero.

Under HC1, each block regression's HCCME is scaled by $n_b / (n_b - k_b)$, where $n_b$ and $k_b$ are that regression's number of observations and number of parameters (cell dummies plus covariates), not those of the full sample.

Because the block regressions only include the covariates that the block keeps, their design matrices have full column rank, and rank deficiency is not a reason for a `missing` edge case standard error. An edge case standard error is `missing` in the following cases only:

- one of the block regressions involved has no residual degrees of freedom (it has as many parameters as observations);
- the estimated variance of one of the four $\hat{\lambda}_{s,t}$ values is exactly zero (for example a binary outcome that takes a single value in one of the cells);
- under HC2, HC3, and HC4, a block regression includes a cell with a single observation (that observation has a leverage of one).

### Edge Case Standard Errors for Two-Way DID-INT with Covariates

For homogeneous, region-varying, and time-varying DID-INT, the covariate term of the [Integrated Long Differences](@ref) cancels, $\widehat{\mathrm{ATT}}(s,t)$ is a combination of four $\hat{\lambda}_{s,t}$ values with fixed coefficients, and the decomposition above gives its variance. For two-way DID-INT with covariates the covariate term does not cancel. Since each cell's regression passes through the cell means, $\hat{\lambda}_{s,t} + \hat{\gamma}_{s,t}'\overline{X}_{s,t} = \overline{Y}_{s,t}$, and with a single control state $c$:

```math
\widehat{\mathrm{ATT}}(s,t) = \left(\overline{Y}_{s,t} - \hat{\delta}'\,\overline{X}_{s,t}\right) - \left(\hat{\lambda}_{s,r} + \hat{\lambda}_{c,t} - \hat{\lambda}_{c,r}\right), \qquad \hat{\delta} = \hat{\gamma}_{s,r} + \hat{\gamma}_{c,t} - \hat{\gamma}_{c,r}
```

The term $\hat{\delta}'\,\overline{X}_{s,t}$ is a product of two estimates: covariate coefficients estimated from three cells, and covariate means estimated from the treated state's cell in period $t$. **DiDInt.jl** treats $\overline{X}_{s,t}$ as an estimate (of the mean of the covariate distribution that the long differences are integrated over) rather than as fixed, and computes:

```math
\mathrm{Var}(\widehat{\mathrm{ATT}}(s,t)) = \mathrm{Var}(\overline{Z}_{s,t}) + \sum_{k \in \mathcal{K}} a' V_k\, a + \mathrm{tr}\!\left(V_{\delta}\, \Sigma\right)
```

The rest of this section first shows where the three terms come from and then defines each of them.

**Where the formula comes from.** Four cells are involved: the treated state's cell in period $t$, which is $(s,t)$, and the three other cells, collected in the set

```math
\mathcal{K} = \{(s,r),\ (c,t),\ (c,r)\}
```

Write each of the two estimates in the product as its true value plus an estimation error: $\hat{\delta} = \delta + e_{\delta}$ and $\overline{X}_{s,t} = \mu + e_{X}$, where $\mu$ is the mean of the covariates in cell $(s,t)$. Multiplying out,

```math
\hat{\delta}'\,\overline{X}_{s,t} = \delta'\,\overline{X}_{s,t} + e_{\delta}'\,\mu + e_{\delta}'\,e_{X}
```

Substituting this into the expression for $\widehat{\mathrm{ATT}}(s,t)$ and collecting the terms by where their estimation error comes from splits it into three parts:

```math
\widehat{\mathrm{ATT}}(s,t) =
\underbrace{\left(\overline{Y}_{s,t} - \delta'\,\overline{X}_{s,t}\right)}_{\text{(i)}}
- \underbrace{\left(\hat{\lambda}_{s,r} + \hat{\lambda}_{c,t} - \hat{\lambda}_{c,r} + e_{\delta}'\,\mu\right)}_{\text{(ii)}}
- \underbrace{e_{\delta}'\,e_{X}}_{\text{(iii)}}
```

1. Part (i) carries the treated cell's own sampling error: the means of the outcome and of the covariates over the observations of cell $(s,t)$, with the coefficients $\delta$ held at their true values.
2. Part (ii) carries the estimation error of the regressions of the three other cells. The error $e_{\delta}$ is made up of the errors in $\hat{\gamma}_{s,r}$, $\hat{\gamma}_{c,t}$ and $\hat{\gamma}_{c,r}$, so part (ii) is, up to a constant, a combination of three quantities, one per cell: the cell's estimated intercept plus its estimated covariate coefficients multiplied by $\mu$. That is each cell's regression line evaluated at the covariate means of the treated cell.
3. Part (iii) is the product of the two errors, the one in $\hat{\delta}$ and the one in $\overline{X}_{s,t}$.

The four cells are estimated by separate regressions on separate sets of observations, so their estimation errors are independent of one another (see [Independence Across Cells](@ref) for when this holds). Part (i) only involves cell $(s,t)$, and part (ii) only involves the three other cells, each through its own regression, so all of these errors are independent. The product in part (iii) has mean zero and is uncorrelated with the other two parts, because each of its two factors has mean zero and is independent of the other. The variance of $\widehat{\mathrm{ATT}}(s,t)$ is therefore the plain sum of the variances of the parts, with no covariance terms, which gives the three terms of the formula.

**The terms of the formula.** In the order in which they appear:

1. The first term, ``\mathrm{Var}(\overline{Z}_{s,t})``, is the variance of part (i): the sampling variation of the treated cell. Here ``\overline{Z}_{s,t}`` is the mean of ``Z_i = Y_i - \hat{\delta}' X_i`` over the observations ``i`` of cell ``(s,t)``, and ``\mathrm{Var}(\overline{Z}_{s,t})`` is the HCCME variance of that mean. Since ``\overline{Z}_{s,t} = \overline{Y}_{s,t} - \hat{\delta}'\,\overline{X}_{s,t}``, it accounts for the sampling variation in ``\overline{Y}_{s,t}`` and in ``\overline{X}_{s,t}``, and for their covariance. In this term ``\hat{\delta}`` is treated as fixed, because its estimation error is handled by the other two terms.

2. The second term, the sum of ``a' V_k\, a`` over the three cells ``k`` in ``\mathcal{K}``, is the variance of part (ii): the uncertainty in what each of the three other cells contributes to the counterfactual. Let ``\hat{\theta}_k = (\hat{\lambda}_k, \hat{\gamma}_k)`` collect the intercept and the covariate coefficients from the regression of cell ``k``, let ``V_k`` be their HCCME covariance matrix, and let ``a = (1, \overline{X}_{s,t})``, a one followed by the covariate means of the treated cell. Then ``\hat{\lambda}_k + \hat{\gamma}_k'\,\overline{X}_{s,t} = a' \hat{\theta}_k`` is the fitted value of the regression of cell ``k`` at the covariate means of the treated cell, and ``a' V_k\, a`` is the variance of that fitted value. The three cells enter ``\widehat{\mathrm{ATT}}(s,t)`` with different signs (plus or minus), but the signs do not matter for the variance, because the cells are independent: the variances simply add up.

3. The third term, ``\mathrm{tr}(V_{\delta}\, \Sigma)``, is the variance of part (iii): the product of the estimation errors in ``\hat{\delta}`` and in ``\overline{X}_{s,t}``. With ``K`` covariates, both matrices are ``K \times K``. The first, ``V_{\delta} = V_{\gamma,(s,r)} + V_{\gamma,(c,t)} + V_{\gamma,(c,r)}``, is the covariance matrix of ``\hat{\delta}``: the sum of the covariate-coefficient blocks of the three matrices ``V_k`` (each ``V_k`` without the row and column of its intercept). The second, ``\Sigma``, is the HCCME covariance matrix of the *mean* ``\overline{X}_{s,t}``, not of the covariates ``X_i`` of individual observations. Finally, ``\mathrm{tr}(A)`` denotes the trace of a square matrix ``A``, the sum of its diagonal entries. With a single covariate all of these are numbers, and the term is just ``\mathrm{Var}(\hat{\delta}) \cdot \mathrm{Var}(\overline{X}_{s,t})``.

The trace appears because the errors $e_{\delta}$ and $e_{X}$ are independent and both have mean zero, so that

```math
\mathrm{Var}\!\left(e_{\delta}'\,e_{X}\right) = \mathbb{E}\!\left[e_{\delta}'\,e_{X}\,e_{X}'\,e_{\delta}\right] = \mathbb{E}\!\left[e_{\delta}'\,\Sigma\,e_{\delta}\right] = \mathrm{tr}\!\left(V_{\delta}\,\Sigma\right)
```

How much the third term matters depends on the number of observations per cell. Both $\mathrm{Var}(\hat{\delta})$ and $\mathrm{Var}(\overline{X}_{s,t})$ shrink roughly in proportion to $1/n$ as the cells grow, as do the first two terms of the formula, so their product shrinks roughly like $1/n^2$. With about 100 observations per cell, the third term is therefore roughly 100 times smaller than the other terms. With about 5 observations per cell, it is only about 5 times smaller and can matter.

A covariate that is dropped from the regression of one of the cells has a coefficient of zero there, and zero rows and columns in the $V_k$ of that cell.

### Independence Across Cells

!!! warning "Panel Data"
    Every edge case standard error assumes that observations in different $(s,t)$ cells are independent of each other. This holds for repeated cross-sections, where each period samples different individuals. It does not hold for panel data, where the same individuals are observed in several periods: the estimates for a state's periods are then correlated, and the formulas above omit that covariance (under two-way DID-INT with covariates, this includes the covariance between $\overline{X}_{s,t}$ and the same state's estimates for the period before treatment). With panel data the edge case standard errors are generally biased.

## Plotting with `didint_plot()`

`didint_plot()` is intended as a simple visual tool. For each requested `ccc` option it plots the $\hat{\lambda}_{s,t}$ values from Step 1 of the [Four Step Estimation Procedure](@ref) (equation 45 of the DID-INT paper) by state and period, or, for event study plots, the averages of the treated states' $\hat{\lambda}_{s,t}$ values by periods since treatment. These show how each form of covariate adjustment affects the trends in the outcome. They are not the [Integrated Long Differences](@ref) used by `didint()`, so the plots are not meant to reproduce the ATT estimates exactly (the two differ when `ccc = "int"` is used with covariates).

## Fast Computation of Randomization Inference and Jackknife

The randomization inference procedure and the fast jackknife both recompute every estimate many times: once for each alternative assignment of the treatment times, or once for each removed state. **DiDInt.jl** does not loop over these cases and rerun the estimation each time. It computes the estimates for all of them together, which is why the procedures remain fast for large values of `nperm` and for many states. Three facts make this possible.

1. **The first step never has to be repeated.** The $\hat{\lambda}_{s,t}$ values (and, for two-way DID-INT with covariates, the covariate coefficients and covariate means) do not depend on which states are treated. They are computed once, in Step 1 of the [Four Step Estimation Procedure](@ref).
2. **No regression is needed.** $\widehat{\mathrm{ATT}}(s,t)$ is the long difference of the treated state minus the weighted mean of the long differences of its controls. This expression equals the Step 3 regression coefficient and can be evaluated directly.
3. **Every estimate is a weighted mean of $\widehat{\mathrm{ATT}}(s,t)$ values.** Once the $\widehat{\mathrm{ATT}}(s,t)$ values are available for every case, the subaggregate ATTs and the aggregate ATT follow from sums and divisions.

### What Is Computed Once

A "cell" below is a pair $(g,t)$ of a treatment time $g$ of the actual assignment and a period $t \geq g$, with $r = g - 1$ the period before that treatment time. For every cell and for every state $i$, treated or not, the following are computed once and stored in tables with one row per state and one column per cell:

- the long difference $d_i = \hat{\lambda}_{i,t} - \hat{\lambda}_{i,r}$;
- the weight $v_i$ that the state would have as a control, the weight $u_i$ that its $\widehat{\mathrm{ATT}}(s,t)$ would have within a subaggregate ATT, and its number of observations $n_{i,t}$ (see [Weighting](@ref));
- whether the state is observed in both periods, so that the long difference exists.

Any state can be treated at $g$ under an alternative assignment, which is why these are needed for every state and not only for the states that are actually treated at $g$.

### All Assignments at Once

The cases to evaluate are written as the rows of an assignment table $L$ with one row per case and one column per state: $L_{b,i}$ is the treatment time of state $i$ in case $b$, or $0$ if the state is never treated in that case.

- For **randomization inference**, each row is one alternative assignment.
- For the **fast jackknife**, every row is the actual assignment, there is one row per state, and row $b$ additionally switches state $b$ off, so that it is neither a treated state nor a control in that row.

For one cell $(g,t)$, let $C$ be the table of zeros and ones that marks the controls of every case ($C_{b,i} = 1$ if state $i$ is a control in case $b$: never treated or, with `notyet = true`, not yet treated in period $t$), and let $d$ and $v$ be the columns of long differences and control weights of that cell. The weighted mean of the controls' long differences is then obtained for all cases with two matrix-vector products:

```math
a = \frac{C\,(v \circ d)}{C\,v}
```

where $\circ$ multiplies element by element and the division is element by element as well, so that $a_b$ is the control mean of case $b$. The ATT of state $i$ in case $b$, if the case treats that state at $g$, is

```math
\widehat{\mathrm{ATT}}_{b}(i,t) = d_i - a_b
```

The sums that make up the subaggregate ATTs are accumulated over the cells in the same way: for each case, the sum of $u_i\,\widehat{\mathrm{ATT}}_{b}(i,t)$ over the states treated at $g$, the sum of their weights $u_i$, and the sum of their $n_{i,t}$. After the last cell, dividing the first sum by the second gives every subaggregate ATT of every case, and the aggregate ATT of each case is the weighted mean of its subaggregate ATTs.

In Julia, the computation for one cell is essentially the following, where `L` has one row per case, and `d`, `v`, `u` hold one value per state:

```julia
# Which states are controls and which are treated, in every case
# (with notyet = true the controls are (L .== 0) .| (L .> t))
controls = (L .== 0)
treated  = (L .== g)

# Control mean of every case, and the ATT of every state in every case
a   = (controls * (v .* d)) ./ (controls * v)
att = d' .- a

# Per case, add to the sums of the group that the cell belongs to
numerator   .+= sum(treated .* u' .* att, dims = 2)
denominator .+= sum(treated .* u', dims = 2)
```

A few products of a table with a vector per cell thus replace one regression per $\widehat{\mathrm{ATT}}(s,t)$ per case. States without a long difference in the cell are left out of both the controls and the treated states, and a case without any control for a cell has no $\widehat{\mathrm{ATT}}(s,t)$ in that cell.

### Two-Way DID-INT with Covariates

Under `ccc = "int"` with covariates, the [Integrated Long Differences](@ref) of the controls are evaluated at the covariate means of the treated state, so the control mean differs from one treated state to the next. It is still linear in those covariate means. With $\Delta\hat{\gamma}_i = \hat{\gamma}_{i,t} - \hat{\gamma}_{i,r}$ the change in the covariate coefficients of state $i$, the control mean that applies to treated state $i$ in case $b$ is

```math
a_b + \beta_b'\,\overline{X}_{i,t}, \qquad \beta_b = \frac{\sum_{c} C_{b,c}\, v_c\, \Delta\hat{\gamma}_c}{\sum_{c} C_{b,c}\, v_c}
```

where $\beta_b$ is the weighted mean of the controls' changes in covariate coefficients. The long difference of the treated state includes its own covariate term, so that

```math
\widehat{\mathrm{ATT}}_{b}(i,t) = \left(d_i + \Delta\hat{\gamma}_i'\,\overline{X}_{i,t}\right) - \left(a_b + \beta_b'\,\overline{X}_{i,t}\right)
```

The term $\Delta\hat{\gamma}_i'\,\overline{X}_{i,t}$ does not depend on the case and is computed once per state and cell. Computing $\beta_b$ for all cases takes one more matrix product, and evaluating it at the covariate means of every state another. For every other `ccc` option the covariate terms cancel and are skipped.

### State-Level Aggregation Options

For the `"state"` and `"sgt"` aggregation options the randomization inference procedure uses, in each alternative assignment, one randomly drawn state in place of each actually treated state (see [Aggregation by State](@ref)). The $\widehat{\mathrm{ATT}}(s,t)$ values of all states are computed as above, and the values of the drawn states are then picked out.

### Memory and the True Jackknife

The alternative assignments are evaluated in batches of 500, so that the size of the tables does not grow with `nperm`.

The true jackknife (`truejack = true`) repeats Step 1 without each state in turn, since under `"hom"`, `"time"` and `"add"` the $\hat{\lambda}_{s,t}$ values of the remaining states change. The first of the three facts above then does not apply, and each removed state requires its own first step, including its own check for collinear covariates (see [Rank Deficiency and Lambda Identification](@ref)). The $\widehat{\mathrm{ATT}}(s,t)$ values and the aggregation are computed in the same way as above, one case at a time.
