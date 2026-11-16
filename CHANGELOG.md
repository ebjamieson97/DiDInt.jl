# Changelog

## 0.10.0

This version changes how the estimates are computed and aggregated, so results differ
from those of 0.9.6 in the cases listed under "Changes to estimates". The documentation
site describes the current procedure in detail.

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
  `seed` therefore differ from earlier versions, and `didint()` no longer changes the
  global random number generator.
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
