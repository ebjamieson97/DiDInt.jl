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
Y_{i,s,t} = \sum_s \sum_t \lambda_{s,t} I(s,t) + f\!\left(X^k_{i,s,t}\right) + \epsilon_{i,s,t} \hspace{0.5cm} \tag{1}
```

where $I(s,t)$ is a dummy variable equal to 1 if observation $i$ belongs to state $s$ in period $t$, and $f\!\left(X^k_{i,s,t}\right)$ is a function of covariates whose form depends on the common causal covariates (CCC) violation that is being accounted for. The estimated $\hat{\lambda}_{s,t}$ values are the covariate-adjusted means of the outcome for state $s$ in period $t$. For common adoption scenarios, the data is flattened into two periods: pre-treatment and post-treatment.

!!! note "Implementation"
    In practice, regression (1) is not estimated directly on the full dataset as this would be computationally prohibitive for datasets with many states, periods, and covariates. Rather, the Frisch-Waugh-Lovell theorem is applied in the case of homogenous DID-INT, state-varying DID-INT, and time-varying DID-INT to recover covariate coefficients which can then be used in conjunction with the outcome means per $(s,t)$ group to retrieve the $\hat{\lambda}_{s,t}$ values. For two-way intersection DID-INT, regression (1) is run separately for each $(s,t)$ group. Finally, for two one-way DID-INT, a sparse matrix procedure is used. The resulting $\hat{\lambda}_{s,t}$ values are numerically equivalent to those from regression (1). More details on estimating $\hat{\lambda}_{s,t}$ can be found in the [Estimating Lambda](@ref) section.


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

In this step weights are used (by default) according to the number of observations associated with each long difference.

**Step 4.** The aggregate ATT is computed as a weighted mean of the subaggregate ATT estimates.

```math
\widehat{\mathrm{ATT}} = \sum^{K}_{k = 1} w_{k}\, \hat{\beta}_{k}
```

where $K$ denotes the total number of subaggregate ATTs. By default the weights used here reflect the number of treated observations associated with each subaggregate ATT.

## Aggregation

## Weighting

## Estimating Lambda

## Randomization Inference

## Jackknife

## Period Grid Construction & Date Matching Procedure