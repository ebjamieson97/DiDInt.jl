# Under the Hood

Most pertinent details of **DiDInt.jl** can be found in the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447), however, some details that might be of particular interest to users of the **DiDInt.jl** package are elaborated upon here.

## Four Step Estimation Procedure

The **DID-INT** estimator is computed in four steps:

1. The means of the outcome of interest (conditional on covariates) are computed for each period for each state (denoted as $\lambda$ in the DID-INT paper).
2. Within state differences between these $\lambda$ values are computed and are labelled as coming from treated or untreated states.
3. A (weighted) regression of within state differences on an intercept term and a dummy variable indicating treated/untreated status is then used to compute subaggregate ATTs for a given subset of the data (more details on the choice of subset in the [Aggregation](@ref) section).
4. The (weighted) mean of subaggregate ATTs is taken in order to compute the aggregate ATT.

More explicitly:

**Step 1.** In principle, the following regression is estimated to obtain conditional outcome means $\hat{\lambda}_{s,t}$:

```math
Y_{i,s,t} = \sum_s \sum_t \lambda_{s,t} I(s,t) + f\!\left(X^k_{i,s,t}\right) + \epsilon_{i,s,t}
```

where $I(s,t)$ is a dummy variable equal to 1 if observation $i$ belongs to state $s$ in period $t$, and $f\!\left(X^k_{i,s,t}\right)$ is a function of covariates whose form depends on the common causal covariates (CCC) violation that is being accounted for (see [Common Causal Covariates and Model Specification](@ref)). The estimated $\hat{\lambda}_{s,t}$ values are the covariate-adjusted means of the outcome for state $s$ in period $t$. For common adoption scenarios, the data is flattened into two periods: pre-treatment and post-treatment.

!!! note "Implementation"
    In practice, this regression is not estimated directly on the full dataset as this would be computationally prohibitive for datasets with many states, periods, and covariates. Rather, the Frisch-Waugh-Lovell theorem is applied in the case of homogeneous DID-INT, state-varying DID-INT, time-varying DID-INT, and two one-way DID-INT to recover covariate coefficients which can then be used in conjunction with the outcome and covariate means per $(s,t)$ group to retrieve the $\hat{\lambda}_{s,t}$ values. For two-way intersection DID-INT, the regression is run separately for each $(s,t)$ group. The resulting $\hat{\lambda}_{s,t}$ values are numerically equivalent to those from the combined regression shown above. More details on estimating $\hat{\lambda}_{s,t}$ can be found in the [Estimating Lambda](@ref) section.

**Step 2.** For each treated state $s$ belonging to treatment cohort $g$, long differences are computed between each post-treatment period $t$ and the period immediately prior to treatment $t^{-g}$:

```math
\widehat{\mathrm{diff}}(s,g,t) = \hat{\lambda}_{s,t} - \hat{\lambda}_{s,t^{-g}}
```

Analogous differences are computed for each control state for each treatment cohort $g$.

**Step 3.** A regression of long differences on an intercept and a treatment indicator is estimated on different subsets of the set of long differences depending on the aggregation method. In general - for more details see [Aggregation](@ref) - the regression is of form:

```math
\widehat{\mathrm{diff}}_{s,g,t} = \alpha + \beta\, d_{s,g,t} + \varepsilon_{s,g,t}
```

where $d_{s,g,t} = 1$ if state $s$ was treated at treatment time $g$ and 0 otherwise. The coefficient $\beta$ estimates the subaggregate ATTs.

In this step, weights are used (by default) according to the number of observations associated with each long difference.

**Step 4.** The aggregate ATT is computed as a weighted mean of the subaggregate ATT estimates.

```math
\widehat{\mathrm{ATT}} = \sum^{J}_{j = 1} w_{j}\, \hat{\beta}_{j}
```

where $J$ denotes the total number of subaggregate ATTs and the weights, $w_j$, are scaled such that $\sum^{J}_{j = 1} w_{j}\, = 1$. By default the weights used here reflect the number of treated observations associated with each subaggregate ATT.

## Common Causal Covariates and Model Specification

A more robust discussion of the common causal covariates assumption can be found in the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447). This section serves only as a quick reference for the different forms of DID-INT. The following table shows the mapping between the DID-INT variation, the functional form of $f\!\left(X^k_{i,s,t}\right)$, and the string value passed to the `ccc` argument in the **DiDInt.jl** package in order to use that form of DID-INT.

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
| By treatment time, $g$                                                                         | Only long differences from treatment time $g$                                                                                                       | `"cohort"` (default) | Yes                | No              |
| By treated state, $s^*$                                                                        | Only long differences from treated state $s^*$ and the corresponding long differences from control states                                           | `"state"`            | Yes                | Yes             |
| By treatment time, $g$, and the post-treatment time, $t$, used in the long difference calculation | Only long differences from that particular $(g,t)$ group                                                                                            | `"simple"`           | Yes                | No              |
| By treated state, $s^*$, and $(g,t)$ group                                                     | Only long differences from treated state $s^*$ in that particular $(g,t)$ group as well as the corresponding long differences from control states | `"sgt"`              | Yes                | No              |
| No grouping                                                                                     | None, computes the aggregate ATT directly from the long differences - only an option for common adoption                                             | `"none"`             | No                 | Yes             |
| By periods since treatment, $p$                                                                | Only long differences in that particular periods-post-treatment group, $p$                                                                          | `"time"`             | Yes                | No              |

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

Optionally (not by default), long differences from treated states, prior to the state's treatment time, can be used as control long differences in the Step 3 computation of subaggregate ATTs. For example, in a scenario with two treated states and 10 periods, State A being treated at period 2 and State B at period 8, then control long differences can be calculated for State B for periods 2-7 (before B's treatment) to match a subset of the long differences calculated for State A.

Using not-yet-treated long differences as controls is available by setting the `notyet` argument to `true`.

## Ensuring Proper Comparisons

In order to best handle data with missing values or inconsistent time period coverage across states, **DiDInt.jl** only includes long differences in the Step 3 regression when there is at least one treated-control pair available for a given $(g,t)$ group. Long differences that do not meet this condition are dropped from any Step 3 regression. This is especially important for the `"cohort"` and `"state"` aggregation methods, which pool long differences across multiple $(g,t)$ groups for a given $g$.

## Period Grid Construction & Date Matching Procedure

**DiDInt.jl** can handle data at various temporal frequencies including daily, weekly, monthly, and yearly observations. For staggered adoption scenarios, the package constructs a period grid and matches observations to discrete periods on that grid.

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


## Standard Errors and P-Values

**DiDInt.jl** reports three types of p-values and two sets of standard errors.

- Standard errors using heteroskedasticity-consistent covariance estimators (HCCME): HC0, HC1 (default), HC2, HC3, or HC4 and p-values via a two-tailed test using a $t$-distribution with $n - k$ degrees of freedom, where $n$ is the number of observations (long differences in Step 3 or subaggregate ATTs in Step 4) and $k$ is the number of estimated regression parameters. The desired HCCME can be selected via the `hc` argument.
- Jackknife standard errors and p-values, described in more detail in the [Jackknife](@ref) section.
- Randomization inference (RI) p-values, via a randomization procedure described in the [Randomization Inference](@ref) section.

There are some circumstances in which the jackknife standard errors and p-values or the randomization inference p-values cannot be calculated. These situations are noted in the [Jackknife](@ref) and [Randomization Inference](@ref) sections, respectively.

In general, the standard errors using HCCMEs for the ATTs calculated in Step 3 are retrieved directly from the Step 3 regression. However, in situations in which there is only one treated state, only one control state, and only one $(g,t)$ group in the Step 3 regression, the model becomes saturated and the regression-based standard error is degenerate. In such cases, it is still possible to recover the heteroskedastic robust standard error (for all aggregation methods besides `"time"`) using the variances and covariances of the $\hat{\lambda}_{s,t}$ terms. Note that given one treated state, one control state, and one $(g,t)$ group, the variance of the ATT computed in Step 3 can be decomposed into variance and covariance terms of the $\hat{\lambda}_{s,t}$ values:

```math
\begin{aligned}
\widehat{\mathrm{ATT}}
&=
(\hat{\lambda}_{s_1,t_1} - \hat{\lambda}_{s_1,t_0})
-
(\hat{\lambda}_{s_2,t_1} - \hat{\lambda}_{s_2,t_0})
\\[6pt]
\mathrm{Var}(\widehat{\mathrm{ATT}})
&=
\mathrm{Var}(\hat{\lambda}_{s_1,t_1} - \hat{\lambda}_{s_1,t_0})
+
\mathrm{Var}(\hat{\lambda}_{s_2,t_1} - \hat{\lambda}_{s_2,t_0})
\\
&\quad
-2\,\mathrm{Cov}\!\left(
\hat{\lambda}_{s_1,t_1} - \hat{\lambda}_{s_1,t_0},
\hat{\lambda}_{s_2,t_1} - \hat{\lambda}_{s_2,t_0}
\right)
\\[6pt]
\mathrm{Var}(\widehat{\mathrm{ATT}})
&=
\mathrm{Var}(\hat{\lambda}_{s_1,t_0})
+
\mathrm{Var}(\hat{\lambda}_{s_1,t_1})
+
\mathrm{Var}(\hat{\lambda}_{s_2,t_0})
+
\mathrm{Var}(\hat{\lambda}_{s_2,t_1})
\\
&\quad
-2\,\mathrm{Cov}(\hat{\lambda}_{s_1,t_0}, \hat{\lambda}_{s_1,t_1})
-2\,\mathrm{Cov}(\hat{\lambda}_{s_2,t_0}, \hat{\lambda}_{s_2,t_1})
\\
&\quad
-2\Big(
\mathrm{Cov}(\hat{\lambda}_{s_1,t_1}, \hat{\lambda}_{s_2,t_1})
-
\mathrm{Cov}(\hat{\lambda}_{s_1,t_1}, \hat{\lambda}_{s_2,t_0})
-
\mathrm{Cov}(\hat{\lambda}_{s_1,t_0}, \hat{\lambda}_{s_2,t_1})
+
\mathrm{Cov}(\hat{\lambda}_{s_1,t_0}, \hat{\lambda}_{s_2,t_0})
\Big)
\end{aligned}
```

Further details on recovering the standard errors in such an edge case can be found in the section [Computation of Edge Case Standard Errors](@ref).

## Jackknife

There are two variations to the **DiDInt.jl** jackknife procedure. The default option is to calculate jackknife estimates of the ATTs while dropping one state at a time directly from Step 3 in the [Four Step Estimation Procedure](@ref). This is referred to as the "fast jackknife". The alternative method, the "true jackknife", calculates jackknife estimates of the ATTs while dropping one state at a time starting from Step 1. The true jackknife can be used in lieu of the fast jackknife by setting the `truejack` argument to `true`.

Note that for the `ccc` options of `"state"` and `"int"` that the true jackknife is always equal to the fast jackknife as the covariate effects are local to the state or state-time, respectively, so that removing a state from the dataset does not affect the estimated covariate effects at other states and thus does not affect the estimated $\hat{\lambda}_{s,t}$ values at other states.

The jackknife standard error is computed as:

```math
\widehat{\mathrm{se}}_{\mathrm{jack}} = \sqrt{\frac{N-1}{N} \sum_{i=1}^{N}\left(\hat{\beta}_{(-i)} - \hat{\beta}\right)^2}
```

where $N$ is the number of states used in the computation of the original subaggregate or aggregate ATT estimate, $\hat{\beta}_{(-i)}$ is the ATT estimate with state $i$ removed, and $\hat{\beta}$ is the original ATT estimate. Jackknife p-values are computed via a two-tailed test using a $t$-distribution with $N - 1$ degrees of freedom. A minimum of two treated states and two control states are required to compute jackknife standard errors and p-values; if this condition is not met, these quantities are returned as `missing`.

Further, note that for the aggregation methods of `"state"` and `"sgt"`, where only one treated state is considered at a time, that the jackknife standard errors and jackknife p-values are never estimated at the subaggregate level. In order for the Step 3 regression to be estimable (without the use of pseudo-inverted matrices) there must be at least one treated state and one control state present in the regression. The structure of `"state"` and `"sgt"` aggregation prevents the jackknife standard error from being fully estimable at the subaggregate level as dropping the one treated state would leave at least one of the remove-one-state subsamples inestimable (the $d_{s,g,t}$ values would all be 0s and the matrix would be rank deficient for the subsample where the treated state is removed). This limitation does **not** prevent the jackknife standard error and jackknife p-value from being calculated for the *aggregate* ATT under `"state"` or `"sgt"` aggregation, so long as there are at least two treated and two control states.

## Randomization Inference

The randomization procedure used in **DiDInt.jl** follows the procedure described in [MacKinnon and Webb (2020)](https://doi.org/10.1016%2Fj.jeconom.2020.04.024), although for the `"state"` and `"sgt"` aggregation options, there is an additional randomization step.

The **DiDInt.jl** implementation of the randomization inference procedure is as follows:

1. The total number of unique treatment assignment randomizations is calculated by $\frac{N!}{(N-n)! \prod_{m \in M} n_m!} - 1$, where $N$ is the total number of states, $n$ is the number of treated states, $M$ is the set of unique treatment times, and $n_m$ is the number of states treated at time $m$. The number of permutations considered during the procedure, $nperm^{*}$, is the lesser of this total and the user-supplied `nperm` value.
2. Next, the order of the set of $N$ total states is randomized, as are the $n$ treatment times. The first $n$ states from the randomized $N$ states are selected and paired with the randomized treatment times, checked to ensure that the new randomized assignment is unique, and then a new vector of treatment indicator values, $d_{s,g,t,r}$ (where $r$ is the $r$-th permutation), is created based on the randomization.
3. Subaggregate and aggregate ATTs are then calculated using the $d_{s,g,t,r}$ values and the distribution of coefficients resulting from the randomized assignments is compared against the coefficients computed when using the actual assignment of treatment times:

$p_{RI} = \frac{\left(\sum_{r=1}^{nperm^{*}} \mathbf{1}\left(|ATT_{RI,r}| > |ATT_{Actual}|\right)\right)}{nperm^{*}}$.

When the randomization inference procedure is performed while the aggregation is set to `"state"` or `"sgt"`, and the initial treatment assignment involved more than one state receiving the same treatment time, there is an additional level of randomization. The set of states in any permutation now associated with the treatment time from the actually treated state are shuffled and the first state is drawn to be used as the treated state in that permutation. This is in order to preserve the structure of the `"state"` and `"sgt"` aggregation methods whereby only one treated state is considered at a time.

!!! note "Computing Subaggregate ATTs During Randomization"
    For all aggregation methods besides `"time"`, in order to avoid having to invert any matrices (which would be more computationally expensive), a difference in weighted means (or a simple difference in means if observation-count weights are not being used) is used in lieu of the Step 3 regression. Letting $\tilde{w}_i = \frac{n_i}{\sum_{h=1}^{H} n_h}$ denote observation-count weights scaled over the full subsample, the subaggregate ATT estimator is:

    For `weighting ∈ {"diff", "both"}`:
    $$\widehat{ATT}_{sub} = \frac{\sum_{i:\, d_i = 1} \tilde{w}_i \,\widehat{\mathrm{diff}}_i}{\sum_{i:\, d_i = 1} \tilde{w}_i} - \frac{\sum_{i:\, d_i = 0} \tilde{w}_i \,\widehat{\mathrm{diff}}_i}{\sum_{i:\, d_i = 0} \tilde{w}_i}$$

    For `weighting ∈ {"att", "none"}`:
    $$\widehat{ATT}_{sub} = \frac{1}{n_1} \sum_{i:\, d_i = 1} \widehat{\mathrm{diff}}_i \;-\; \frac{1}{n_0} \sum_{i:\, d_i = 0} \widehat{\mathrm{diff}}_i$$

    where $n_1$ and $n_0$ are the number of treated and control observations in the subsample, respectively.


It should also be noted that the randomization inference procedure will not run unless it is possible to calculate the long differences for each $(g,t)$ used in the computation of the actual ATT values for each state. This may occur in datasets with missing values or inconsistent time period coverage across states. This safeguard is to ensure that the randomization procedure only runs if it is possible to make a proper comparison between the actual $ATT$ values and the $ATT_{RI}$ values.

## Estimating Lambda

As mentioned previously, estimating the regression model shown in Step 1 would be computationally infeasible, especially for a dataset with many states, time periods, and covariates. For example, the number of columns in the design matrix for a two-way intersection DID-INT model would be $n_{s} \cdot n_{t} + n_{s} \cdot n_{t} \cdot n_{k}$ where $n_{s}$ is the number of states, $n_{t}$ is the number of time periods, and $n_{k}$ is the number of covariates. Consequently, **DiDInt.jl** takes three different approaches for estimating the $\hat{\lambda}_{s,t}$ values depending on the specific DID-INT variation. For each approach, the resulting $\hat{\lambda}_{s,t}$ values are numerically equivalent to those from the combined regression formulation shown in Step 1.

### Two-Way Intersection DID-INT

Two-way intersection DID-INT is the simplest case. As the effects of the covariates are local to each $(s,t)$ group, separate regressions of the form:

```math
Y_{i,s,t} = \lambda_{s,t} + \sum_{k=1}^{K} \gamma^k_{s,t} X^k_{i,s,t} + \epsilon_{i,s,t}
```

are estimated for each $(s,t)$ cell, where $\hat{\lambda}_{s,t}$ is simply the estimated intercept for the local regression. 

If, for any $(s,t)$ group, a covariate causes the design matrix to be rank-deficient, then that covariate is dropped from the regression for that $(s,t)$ group.

### State-Varying, Time-Varying, and Homogeneous DID-INT

For state-varying, time-varying, and homogeneous DID-INT, the Frisch-Waugh-Lovell (FWL) theorem is applied. The procedure is the same across all three cases, differing only in the grouping over which the covariate coefficients $\hat{\beta}$ are estimated:

1. Both the outcome and covariates are demeaned within each $(s,t)$ cell - which is equivalent to residualizing the outcome and covariates by the $I(s,t)$ cell dummies.
2. The covariate coefficients $\hat{\beta}$ are estimated by regressing the demeaned outcome on the demeaned covariates (without an intercept), where the regression is run:
   - Globally across all observations for homogeneous DID-INT
   - Within each time period $t$ for time-varying DID-INT
   - Within each state $s$ for state-varying DID-INT
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

## Computation of Edge Case Standard Errors

In general, the heteroskedasticity-robust standard errors for the subaggregate ATTs computed in Step 3 are retrieved directly from the Step 3 regression. However, as noted in the [Standard Errors and P-Values](@ref) section, when there is only one treated state, one control state, and one $(g,t)$ group in the Step 3 regression the model is saturated and the Step 3 regression-based standard error is degenerate. In such cases, the variance of the ATT can instead be recovered from the variances and covariances of the $\hat{\lambda}_{s,t}$ values, as shown in the decomposition in the aforementioned section. Recovering the variance of the ATT from the variances and covariances of the $\hat{\lambda}_{s,t}$ values necessitates computing the full variance-covariance matrix of the $\hat{\lambda}_{s,t}$ values, $\widehat{\mathrm{Var}}(\hat{\lambda})$.

In all DID-INT variations, the computation of the edge case standard errors is a more computationally expensive procedure than the standard lambda estimation procedure, so by default **DiDInt.jl** does not compute them. Computation of edge case standard errors can be toggled on by setting `edgecase = true`. When edge case computation is toggled on, the set of required $\hat{\lambda}_{s,t}$ values needed to compute ATTs is determined prior to any computations, and those involved in the saturated Step 3 regression(s) are identified. Only the necessary blocks of data are sent to the edge case computation procedure, where both the $\hat{\lambda}_{s,t}$ values and their variances/covariances are recovered. The separation of the data into blocks sent to the edge case procedure, and the rest of the data which is fed into the standard [Estimating Lambda](@ref) procedure is determined by the $(s,t)$ groups to which the implicated $\hat{\lambda}_{s,t}$ belong to as well as the specific DID-INT variation. For example, time-varying DID-INT requires that all observations in group $t$ are grouped together, so any observation belonging to a $t$ period implicated in a saturated Step 3 regression is sent to the edge case procedure. Parallel logic applies to the other DID-INT variations, however it should be noted that the edge case procedure is not an option for two one-way DID-INT as this model is always rank deficient and $\widehat{\mathrm{Var}}(\hat{\lambda})$ is not recoverable. 

### Computing $\widehat{\mathrm{Var}}(\hat{\lambda})$

The variance-covariance matrix $\widehat{\mathrm{Var}}(\hat{\lambda})$ is an $n_\lambda \times n_\lambda$ matrix, where $n_\lambda = n_s \cdot n_t$ is the total number of $\hat{\lambda}_{s,t}$ values. For homogeneous DID-INT it is computed by following a single regression, and for all other DID-INT variations it is computed in blocks following regressions on subsets of the data. The raw outcome $Y_{i,s,t}$ is regressed on $I(s,t)$ cell dummies and the original (non-demeaned) covariates $X$, after which the covariance matrix of $\hat{\lambda}_{s,t}$ values from that regression is recovered using the chosen HCCME. 

The specific regression used to compute $\widehat{\mathrm{Var}}(\hat{\lambda})$ mirrors the CCC structure of the chosen DID-INT variant:

- **Two-way intersection** (`"int"`): a separate regression is run for each $(s,t)$ cell, with the cell dummy and covariates as regressors. The HCCME for each cell's dummy coefficient is placed into the corresponding diagonal entry of $\widehat{\mathrm{Var}}(\hat{\lambda})$. Off-diagonal entries between different cells are zero since the regressions are run separately.
- **State-varying** (`"state"`): a separate regression is run for each state $s$, pooling all time periods within that state. The HCCME for the cell dummies within each state block is placed into the corresponding entries of $\widehat{\mathrm{Var}}(\hat{\lambda})$ (including off-diagonal covariance terms).
- **Time-varying** (`"time"`): analogous to state-varying, with separate regressions run for each time period $t$.
- **Homogeneous** (`"hom"`): a single global regression is run across all observations. The HCCME for all cell dummies from the regression is exactly $\widehat{\mathrm{Var}}(\hat{\lambda})$.

!!! note "Consistency with Lambda Estimation"
    For the FWL-based cases (`"hom"`, `"state"`, `"time"`), covariates that are constant within every $(s,t)$ cell in the relevant block are dropped from the block regressions used in the edge case procedure. These same covariates are demeaned to exactly zero during the regular lambda estimation procedure and thus contribute nothing to the estimated $\hat{\beta}$ there either. For the `"int"` case, both procedures are identical: any covariate that does not vary within a given $(s,t)$ cell is dropped from that cell's regression in both the normal lambda estimation procedure and the edge case procedure.

!!! note "Rank Deficiency and Lambda Identification"
    In all DID-INT variations besides two-way intersection DID-INT, the design matrices used to recover $\hat{\lambda}_{s,t}$ values are allowed to be rank deficient, meaning $\hat{\lambda}_{s,t}$ is not always uniquely identified. However, the ATT estimates remain uniquely identified due to cancellation that occurs in the long difference contrasts. **DiDInt.jl** makes use of the `\` operator in order to facilitate this, which returns the minimum-norm least squares solution, ensuring consistent ATT estimates between procedures, even when $\hat{\lambda}_{s,t}$ is not uniquely identified. In the case of two-way intersection DID-INT, in both the normal lambda estimation procedure and the edge case procedure, a rank revealing QR procedure is used to drop any remaining rank deficiency inducing covariates within each $(s,t)$ block after having already removed covariates that did not vary within the given $(s,t)$ block.