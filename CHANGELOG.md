# Changelog

## 0.11.0

### Weighting options

- The `treated_weight` argument is removed. The ATT of a treated state in a period is
  weighted within its sub-aggregate ATT by that state's number of observations in that
  period under `weighting = "both"` and `"att"`, and equally under `"diff"` and `"none"`.
  This is what the default `treated_weight = "treated"` did, so results obtained with the
  default are unchanged; `treated_weight = "population"` is no longer available.
- Each value of `weighting` now has one meaning: `"diff"` weights the long differences of
  the control states, `"att"` weights the ATTs (within each sub-aggregate ATT and within
  the aggregate ATT), `"both"` does both and `"none"` neither.

### Edge case standard errors under `ccc = "int"`

- With covariates, these standard errors are now conditional on the observed covariates:
  the variance of the ATT is the sum, over the four state-period cells involved, of the
  variance of each cell's fitted value at the covariate means of the treated state, with
  those means treated as fixed. 0.10.1 also included the sampling variability of the
  covariate means. Only standard errors computed with `edgecase = true` under
  `ccc = "int"` with covariates change; no ATT changes.

### Warnings about collinear covariates and rank deficiency

- Under `ccc = "int"`, a warning is now given only for the ATTs that depend on a dropped
  covariate, naming the ATT, the cell and the covariate. A cell's coefficients are
  evaluated at the covariate means of the treated state, and a covariate dropped from the
  cell's regression (because it does not vary in the cell, or is collinear with
  covariates listed before it) matters only when the cell's data cannot determine its
  fitted value at those means. In 0.10.1 every covariate dropped as collinear was
  reported, whether or not it affected an ATT.
- Under `ccc = "add"`, the warning about rank deficiency beyond the split between states
  and periods is given only when an ATT depends on which least-squares solution of the
  first stage is used. Every direction in which the solution is undetermined is checked.
  In 0.10.1 the warning was given whenever the rank was below the expected rank, also
  when no ATT was affected.
- Under `"hom"`, `"state"` and `"time"` the warning is unchanged.

### Bug fixes

- Under `ccc = "add"`, which columns of the first-stage regression count as linearly
  dependent is now decided by the same relative tolerance as the collinearity rule of
  the other options, after scaling every column to unit length. With a covariate that has
  a very large mean and is collinear with another, the rounding error left by demeaning
  could count as a covariate of its own, and the estimates then changed with the order
  of the rows of the data (and could differ between platforms). Estimates are otherwise
  unchanged.

### Documentation

- The estimation procedure is described in five steps: the average of the ATTs within
  each sub-aggregate ATT and the average of the sub-aggregate ATTs are now separate steps.
- Shorter, rewritten sections on weighting, rank deficiency, the edge case standard
  errors of two-way DID-INT with covariates, how alternative assignments are drawn, and
  how randomization inference and the jackknife are computed quickly. The example of the
  aggregate ATT across aggregation options now uses the default weighting.

## 0.10.1

### Naming

- The model with state-specific covariate coefficients is now called **region-varying
  DID-INT** (previously state-varying), and the model with coefficients specific to each
  state and period **two-way DID-INT**. The printed summary uses these names. The `ccc`
  values (`"state"`, `"int"`) and every other argument and column name are unchanged.
- DID-INT is spelled out as "integrated difference-in-differences" throughout.

### Collinear covariates

- **Collinear covariates are dropped in the order of `covariates`** under `ccc = "hom"`,
  `"state"`, `"time"` and `"int"`: a covariate is dropped where it is a linear combination
  of the covariates listed before it, in the part of the data over which its coefficient
  is estimated (all data, a state, a period, or a state-period cell). The rule uses a
  relative tolerance and no column pivoting, so the estimates no longer depend on the
  units in which covariates are measured, and they can be changed only by reordering
  `covariates`.
- **A warning names the covariates that were dropped** as collinear and where. Covariates
  that do not vary within any state-period cell of that part of the data are still dropped
  silently, as they cannot affect the estimates.
- **Changes to estimates.** Estimates change only where covariates are collinear:
  - Under `"hom"`, `"state"` and `"time"`, no covariate was dropped in 0.10.0. Collinear
    covariates have infinitely many coefficient vectors that fit equally well, and
    0.10.0 used the one with the smallest sum of squared coefficients (the minimum-norm
    least-squares solution). For example, with `x2 = 2 * x1` and a combined effect of 5
    per unit of `x1`, it gave `x1` a coefficient of 1 and `x2` a coefficient of 2. That
    choice changes when a covariate is rescaled (dollars to thousands of dollars), and so
    did the ATTs. 0.10.1 drops `x2` and gives `x1` a coefficient of 5, whatever the units.
  - Under `"int"`, 0.10.0 already dropped collinear covariates within a cell, but chose
    which one by a column-pivoted QR decomposition, which favours the covariates with the
    most variation and therefore also depended on units. Estimates change in cells where
    it dropped a covariate other than the one listed last.

  Estimates without collinear covariates, and all estimates under `"add"`, are unchanged.
- **`edgecase = true` no longer changes any ATT.** With `edgecase = true`, the cells
  involved in an edge case are estimated by a separate regression that also yields the
  covariance matrix needed for the standard error. In 0.10.0 that regression handled
  collinear covariates differently from the regression used for all other cells, so with
  collinear covariates the ATTs could depend on whether `edgecase` was on. Both now use
  the same covariates, determined once, and return identical ATTs.
- **Fewer missing edge case standard errors.** An edge case standard error was `missing`
  in 0.10.0 when collinear covariates made that separate regression rank deficient. It is
  now computed.
- **`ccc = "add"` warns about unusual rank deficiency.** Its first-stage regression is
  always rank deficient in a way that does not affect the estimates. A warning is now
  given when it is rank deficient beyond that (for example, a covariate that does not vary
  within any cell of some state or period), as the ATTs may then depend on the
  least-squares solution used.

### Edge case standard errors under `ccc = "int"` with covariates

- **Bug fix: these standard errors were slightly too large.** The variance of the
  ATT of a treated state in a period with a single control state has three terms, the third
  being a trace term for the product of two estimation errors. The first two terms are
  computed with estimates in place of unknown true values, which already makes each of
  them too large by that trace term on average. 0.10.0 added the trace term on top; it is
  now subtracted, which gives an estimate that is right on average. The difference is of
  the order of a few percent of the variance in cells with around ten observations and
  negligible in large cells. Only standard errors computed with `edgecase = true` under
  `ccc = "int"` with covariates change; no ATT changes. The Details page has the full
  derivation.

### Jackknife warnings

- When removing a state leaves a sub-aggregate ATT (or a policy) without any treated
  state, its jackknife standard error is `missing` and no warning is given. This reverses
  the warning introduced in 0.10.0. The jackknife standard error of the aggregate ATT (or
  the pooled ATT) is still reported.
- The warning that the jackknife of the aggregate ATT averages over different sets of
  sub-aggregate ATTs is removed.

### Documentation

- This changelog is now a page of the documentation site.
- The Details page has rewritten sections on weighting (as three successive averages),
  on the edge case standard errors of two-way DID-INT with covariates, and on how
  alternative assignments are drawn for randomization inference, and new sections on rank
  deficiency and on how randomization inference and the jackknife are computed quickly.
- The Examples page has new examples on cohort ATTs with and without policies and on the
  aggregate ATT across aggregation options, and explains when the parallel trends plots
  are valid.
- The `edgecase` documentation states that edge case standard errors are generally biased
  with panel data.

## 0.10.0

This version changes how the estimates are computed and aggregated, so results differ
from those of 0.9.6 in the cases listed under "Changes to estimates". The documentation
site describes the current procedure in detail.

### `didint()` returns a `DiDIntResult`

- `didint()` now returns a `DiDIntResult` instead of a `DataFrame`. Displaying it prints a
  summary of the specification and the estimates, and its fields hold the estimates as
  tables: `agg` (the aggregate ATT, or the ATT of each policy and the pooled ATT), `sub`
  (the sub-aggregate ATTs), `att_st` (the ATT of each treated state in each period) and
  `specs` (the specification).
- `coef`, `stderror` and `coefnames` return the estimates, their standard errors and
  their labels as vectors.
- **To update existing code:** `DataFrame(didint(...))` returns the table that earlier
  versions returned, with the same column names.
- The Stata and R wrappers (`wrapper = "stata"` or `"r"`) still receive that table.

### Changes to estimates

- **ATT for each treated state and period first.** An ATT is now estimated for every
  treated state in every period from its treatment time on, and every sub-aggregate ATT is
  a weighted mean of those. The aggregate ATT is a weighted mean of the sub-aggregate ATTs.
  Previously each sub-aggregate ATT came from a regression that pooled long differences.
  - `agg = "sgt"` estimates and standard errors are unchanged.
  - `agg = "simple"` estimates are unchanged with `treated_weight = "population"`; its
    standard errors change for groups with more than one treated state.
  - `agg = "cohort"`, `"state"` and `"time"` estimates change in general.
- **New `treated_weight` argument**, with default `"treated"`: an ATT is weighted within
  its sub-aggregate by the treated state's observations in that period. `"population"`
  gives the previous weighting (the observations of that period and of the period before
  treatment).
- **Integrated long differences** (Karim & Webb 2025, eqs. 46-47). Under `ccc = "int"`
  with covariates, long differences are evaluated at the treated state's covariate means,
  which changes the estimates. Estimates under the other `ccc` options, and without
  covariates, are unaffected.
- **Common adoption is handled like staggered adoption**, as the case of a single
  treatment cohort. The data is no longer collapsed into one pre-treatment and one
  post-treatment period: long differences are taken from the period immediately before
  treatment, an ATT is estimated for each post-treatment period, and every `agg` option is
  available. All common adoption results change.
- **Bug fix:** with staggered adoption, the long differences that treated states
  contribute to other treatment cohorts were computed for the wrong periods. This affected
  estimates with `notyet = true` and all randomization inference p-values.

### New features

- `agg = "calendar"`: sub-aggregate ATTs by period.
- **Policies.** New arguments `pvar` and `treated_policies` give the policy of each
  treated state, and new aggregations `"policy"`, `"policy_gt"`, `"policy_cohort"`,
  `"policy_calendar"` and `"policy_time"` report an ATT for each policy (`att_p`).
  `policy_aggregation = true` adds a pooled ATT across policies (`att_pool`), weighted as
  set by `policy_weighting`.
- **Exact randomization inference.** When no more than `nperm` alternative treatment
  assignments exist, all of them are used.
- **Kernel-smoothed randomization inference p-values** (Racine & MacKinnon 2007) when
  `ri_kernel_cutoff` or fewer permutations are used (default 99). Set `ri_kernel_cutoff`
  to `nothing`, `false` or `0` to turn this off.
- **`edgecase` standard errors** now apply to any ATT of a treated state and period that
  has a single control state, whatever the `agg` option (including `"time"`), and account
  for the covariate terms under `ccc = "int"`.

### Changes to inference

- The randomization inference p-value is now (1 + number of placebo estimates at least as
  large in absolute value as the estimate) / (number of permutations + 1). Previously it
  was the share of placebo estimates strictly larger than the estimate. Ties count against
  rejecting, with a small tolerance for rounding.
- Permutations are drawn by a faster method and from a random number
  generator that only depends on `seed`. Randomization inference p-values for a given
  `seed` therefore differ from earlier versions, and when a `seed` is given `didint()`
  no longer changes the global random number generator.
- Under `agg = "state"` and `"sgt"`, each treated state is represented by the same
  randomly drawn state in all of its periods within a permutation.
- A jackknife standard error is reported for an estimate only if it involves at least two
  treated states and two control states.
- HC2, HC3 and HC4 standard errors of the ATT of a single treated state and period are
  `NaN` (the treated state has leverage one).

### Removed

- `agg = "none"`: every `agg` option reports the aggregate ATT in `agg_att`.
- The `verbose` argument: the randomization inference no longer prints its progress.

### Other

- The jackknife and the randomization inference are computed for all deletions and
  permutations at once, and the long differences are built from lookups, which makes
  `didint()` several times faster.
- Bug fixes: `truejack = true` with `agg = "time"` no longer errors; string covariates
  with exactly two categories are handled correctly; `edgecase = true` with `ccc = "add"`
  is switched off with a warning.

### Bug fixes that can change estimates

These affect particular kinds of data; results for other data are unchanged by them.

- **Covariates measured at the state-period level.** A covariate that is constant within
  every state-period cell should have no effect, but under `ccc = "hom"`, `"state"`,
  `"time"` and `"add"` it could change the estimates substantially (rounding noise left
  by the within-cell demeaning was being fitted). Such covariates now have no effect.
- **Categorical covariates under `ccc = "int"`.** When the reference category of a
  categorical covariate does not occur in a state-period cell, that cell's dummy variables
  add up to one. They were all kept, which made the estimates depend on the level of the
  outcome. One of them is now dropped in that cell.
- **Treated states with a treatment time outside the data.** A state treated at or before
  the first date of the data (or of `start_date`) was silently used as a control. Its
  observations are now dropped, with a warning. A state treated after the last date is
  used as a control, as before, now with a warning; previously a treatment time shortly
  after the last date could be counted as occurring in the last period.
- **Automatic period length.** For data spaced several years apart (for example every
  five or ten years) the detected period could be one year too long. Specifying `freq` and
  `freq_multiplier` was not affected.
- **Randomization inference with incomplete data.** The randomization inference could run
  although a treated state was not observed in a period that it needs under alternative
  assignments. Its p-values are now `missing` in that case, with a warning.

### Other bug fixes

- `hc = "hc4"` no longer errors under `weighting = "both"` or `"diff"`.
- Listing a treated state more than once, or a `gvar` column that is not constant within a
  state, is now an error instead of silently changing the estimates.
- `gvar` works with a `Date` time column, and with rows whose state is missing.
- Numbers given as floats are accepted for `hc`, `freq_multiplier`, `start_date` and
  `end_date`; `date_format` is not case sensitive and an unrecognized format is an error.
- `didint_plot(event = true)` returns its rows sorted by `ccc` and by period relative to
  treatment.
- `didint_plot()` accepts a single treatment time for several treated states, a single
  treated state, and covariates given as symbols.
- A jackknife deletion that removes every treated state of a sub-aggregate ATT now gives a
  warning.
- The kernel-smoothed p-value falls back to the unsmoothed one when the placebo estimates
  are equal up to rounding.
