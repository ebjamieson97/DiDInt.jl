# Under the Hood

Most pertinent details of **DiDInt.jl** can be found in the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447), however, some details that might be of particular interest to users of the **DiDInt.jl** package are elaborated upon here.

## Four Step Estimation Procedure

The **DID-INT** estimator is computed in four steps:

1. The means of the outcome of interest (conditional on covariates) are computed for each period for each state (denoted as $\lambda$ in the DID-INT paper).
2. Within state differences between these $\lambda$ values are computed and are labelled as coming from treated or untreated states.
3. A (weighted) regression of differences on an intercept term and a dummy variable indicating treated/untreated status is then used to compute subaggregate ATTs for a given subset of the data (more details on the choice of subset in the aggregation section).
4. The (weighted) mean of subaggregate ATTs is taken in order to compute the aggregate ATT.

More explicitly:

**Step 1.** In principle, the following regression is estimated to obtain conditional outcome means $\hat{\lambda}_{s,t}$:

```math
Y_{i,s,t} = \sum_s \sum_t \lambda_{s,t} I(s,t) + f\!\left(X^k_{i,s,t}\right) + \epsilon_{i,s,t}
```

where $I(s,t)$ is a dummy variable equal to 1 if observation $i$ belongs to state $s$ in period $t$, and $f\!\left(X^k_{i,s,t}\right)$ is a function of covariates whose form depends on the common causal covariates (CCC) violation that is being accounted for (see [Common Causual Covariates and Model Specification](@ref)). The estimated $\hat{\lambda}_{s,t}$ values are the covariate-adjusted means of the outcome for state $s$ in period $t$. For common adoption scenarios, the data is flattened into two periods: pre-treatment and post-treatment.

!!! note "Implementation"
    In practice, this regression is not estimated directly on the full dataset as this would be computationally prohibitive for datasets with many states, periods, and covariates. Rather, the Frisch-Waugh-Lovell theorem is applied in the case of homogenous DID-INT, state-varying DID-INT, and time-varying DID-INT to recover covariate coefficients which can then be used in conjunction with the outcome and covariate means per $(s,t)$ group to retrieve the $\hat{\lambda}_{s,t}$ values. For two-way intersection DID-INT, the regression is run separately for each $(s,t)$ group. Finally, for two one-way DID-INT, a sparse matrix procedure is used. The resulting $\hat{\lambda}_{s,t}$ values are numerically equivalent to those from the combined regression shown above. More details on estimating $\hat{\lambda}_{s,t}$ can be found in the [Estimating Lambda](@ref) section.

**Step 2.** For each treated state $s$ belonging to treatment cohort $g$, long differences are computed between each post-treatment period $t$ and the period immediately prior to treatment $t^{-g}$:

```math
\widehat{\mathrm{diff}}(s,g,t) = \hat{\lambda}_{s,t} - \hat{\lambda}_{s,t^{-g}}
```

Analogous differences are computed for each control state for each treatment cohort $g$.

**Step 3.** A regression of long differences on an intercept and a treatment indicator is estimated on different subsets of the long differences depending on the aggregation method. In general - for more details see [Aggregation](@ref) - the regression is of the general form:

```math
\widehat{\mathrm{diff}}_{s,g,t} = \alpha + \beta\, d_{s,g,t} + \varepsilon_{s,g,t}
```

where $d_{s,g,t} = 1$ if state $s$ belongs to treatment cohort $g$ and 0 otherwise. The coefficient $\beta$ estimates the subaggregate ATTs.

In this step, weights are used (by default) according to the number of observations associated with each long difference.

**Step 4.** The aggregate ATT is computed as a weighted mean of the subaggregate ATT estimates.

```math
\widehat{\mathrm{ATT}} = \sum^{J}_{j = 1} w_{j}\, \hat{\beta}_{j}
```

where $J$ denotes the total number of subaggregate ATTs and the weights, $w_j$, are scaled such that $\sum^{J}_{j = 1} w_{j}\, = 1$. By default the weights used here reflect the number of treated observations associated with each subaggregate ATT.

## Common Causual Covariates and Model Specification

A more robust discussion of the common causal covariates assumption can be found in the the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447). This section serves only as a quick reference for the different forms of DID-INT. The following table shows the mapping between the DID-INT variation, the functional form of $f\!\left(X^k_{i,s,t}\right)$, and the string value passed to the `ccc` argument in the **DiDInt.jl** package in order to use that form of DID-INT.

| DID-INT Variation            | Functional Form                                                                                                                          | `ccc`     |
| :--------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------- | :---------- |
| Homogeneous DID-INT          | $f(X_{i,s,t}) = \sum_{k=1}^{K} \gamma^k X^k_{i,s,t}$                                                                                   | `"hom"`   |
| State-Varying DID-INT        | $f(X_{i,s,t}) = \sum_{s=1}^{S} \sum_{k=1}^{K} \gamma^k_s I(s) X^k_{i,s,t}$                                                             | `"state"` |
| Time-Varying DID-INT         | $f(X_{i,s,t}) = \sum_{t=1}^{T} \sum_{k=1}^{K} \gamma^k_t I(t) X^k_{i,s,t}$                                                             | `"time"`  |
| Two-Way Intersection DID-INT | $f(X_{i,s,t}) = \sum_{s=1}^{S} \sum_{t=1}^{T} \sum_{k=1}^{K} \gamma^k_{s,t} I(s)I(t) X^k_{i,s,t}$                                      | `"int"`   |
| Two One-Way DID-INT          | $f(X_{i,s,t}) = \sum_{s=1}^{S} \sum_{k=1}^{K} \gamma^k_s I(s) X^k_{i,s,t} + \sum_{t=1}^{T} \sum_{k=1}^{K} \gamma^k_t I(t) X^k_{i,s,t}$ | `"add"`   |

## Aggregation

In Step 3 of the [Four Step Estimation Procedure](@ref), a choice needs to be made in terms of which long differences should be included in the regression. For example, we could include only the long differences associated with treatment time $g$ which would yield the subaggregate ATT for $g$, $\hat{\beta}_{g}$. Repeating this for each treatment time would give us the set of subaggregate ATTs based on `"cohort"` aggregation. Conversely, we could restrict the long differences to a particular treated state, $s^*$, and the corresponding long differences from the control states in order to retrieve the subaggregate ATTs per treated state. The following table gives an overview of the aggregation groupings, corresponding restrictions, the string value passed to the `agg` argument in the **DiDInt.jl** package, and whether the aggregation method is applicable to staggered adoption or common adoption scenarios.

| Grouping                                                                                        | Restriction                                                                                                                                          | `agg`                | Staggered Adoption | Common Adoption |
| :---------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------- | :----------------- | :-------------- |
| By treatment time,$g$                                                                         | Only long differences from treatment time$g$                                                                                                       | `"cohort"` (default) | Yes                | No              |
| By treated state,$s^*$                                                                        | Only long differences from treated state$s^*$ and the corresponding long differences from control states                                           | `"state"`            | Yes                | Yes             |
| By treatment time and the post-treatment time used in the long difference calculation,$(g,t)$ | Only long differences from that particular$(g,t)$ group                                                                                            | `"simple"`           | Yes                | No              |
| By treated state$s^*$ and $(g,t)$ group                                                     | Only long differences from treated state$s^*$ in that particular $(g,t)$ group as well as the corresponding long differences from control states | `"sgt"`              | Yes                | No              |
| No grouping                                                                                     | None, computes the aggregate ATT directly from the long differences - only an option for common adoption                                             | `"none"`             | No                 | Yes             |
| By periods since treatment,$p$                                                                | Only long differences in that particular periods-post-treatment group,$p$                                                                          | `"time"`             | Yes                | No              |

!!! note "Aggregation by Periods Since Treatment"
    Note that when using the `"time"` aggregation option, the functional form of the long difference regression is modified by adding dummy variables for each treatment cohort (minus a reference cohort) in order to account for cohort-specific intercepts. The modified regression for `"time"` aggregation is then given by: $\widehat{\mathrm{diff}}_{s,g,t,p} = \alpha + \beta\, d_{s,g,t,p} + \sum_{g = 2}^{G} \phi_{g}I(g) + \varepsilon_{s,g,t,p}$.

## Weighting

Weighting can occur at Step 3 (when computing subaggregate ATTs) and at Step 4 (when computing the aggregate ATT) in **DiDInt.jl**, and by default occurs at both steps. The weighting that occurs when computing the subaggregate ATTs weights each long difference by the number of observations used in the calculation of the long difference. The weighting that occurs when computing the aggregate ATT uses weights according to the number of *treated* observations - that is, observations from treated states in post-treatment periods - used in the calculation of each subaggregate ATT. The following table shows the available weighting options in **DiDInt.jl**

| `weighting`        | Applies weights when...                            |
| :------------------- | :------------------------------------------------- |
| `"both"` (default) | Computing subaggregate ATTs and the aggregate ATT. |
| `"att"`            | Computing the aggregate ATT.                       |
| `"diff"`           | Computing the subaggregate ATTs.                   |
| `"none"`           | No weighting is applied.                           |

## Not Yet Treated Cells

Optionally (not by default), long differences from treated states, prior to the state's treatment time, can be used as control long differences in the Step 3 computation of subaggregate ATTs. For example, in a scenario with two treated states and 10 periods, State A being treated at period 2 and State B at period 8, then control long differences can be calculated for State B for periods 2-7 (before B's treatment) to match the long differences calculated for State A.

Using not-yet-treated long differences as controls is available by setting the `notyet` argument to `true`.

## Ensuring Proper Comparisons

In order to best handle data with missing values or inconsistent time period coverage across states, **DiDInt.jl** only includes long differences in the Step 3 regression when there is at least one treated-control pair available for a given $(g,t)$ group. Long differences that do not meet this condition are dropped from any Step 3 regression. This is especially important for the `"cohort"` and `"state"` aggregation methods, which pool long differences across multiple $(g,t)$ groups for a given $g$.

## Period Grid Construction & Date Matching Procedure

**DiDInt.jl** can handle data at various temporal frequencies including daily, weekly, monthly, and yearly observations. For staggered adoption scenarios where treatment timing matters, the package constructs a period grid and matches observations to discrete periods on that grid.

### Automatic Period Detection

By default, if `freq` is not specified, **DiDInt.jl** automatically detects the appropriate period length from the data by examining the maximum distance between adjacent time observations and identifies whether the data is best represented at a yearly, monthly, weekly, or daily frequency (or some combination thereof) and constructs period lengths accordingly. This automatic detection accounts for leap years.

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

## Randomization Inference

The randomization procedure used in **DiDInt.jl** follows the procedure described in [MacKinnon and Webb (2020)](https://doi.org/10.1016%2Fj.jeconom.2020.04.024), although there are a few small differences related to some of the nuances introduced due to some of the aggregation methods.

The **DiDInt.jl** implementation of the randomization inference procedure is as follows:

1. The total number of unique treatment assignment permutations is calculated as: $\frac{N!}{(N-n)! \prod_{m \in M} n_m!}$, where $N$ is the total number of states, $n$ is the number of treated states, $M$ is the set of unique treatment times, and $n_{m}$ is the number of states treated at time $m$. The term $\frac{N!}{(N-n)!}$ counts the number of ways to assign $n$ states from $N$ total states to the initially specified treatment times, while the division by $\prod_{m \in M}n_m!$ removes the overcounting that arises when there are treatment times that were initially assigned to more than one treated state.
2. 

## Jackknife

## Estimating Lambda

## Computation of Edge Case Standard Errors
