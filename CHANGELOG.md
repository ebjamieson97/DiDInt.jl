# Changelog

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
