###########################################################################################
# estimate_didint.jl: the DID-INT estimator                                               #
#                                                                                         #
#   Section 1: didint(), the main function                                                #
#   Section 2: Long differences                                                           #
#   Section 3: ATT(s,t), sub-aggregates and the aggregate ATT                             #
#   Section 4: Regressions and edge-case standard errors                                  #
#   Section 5: Jackknife                                                                  #
#   Section 6: Randomization inference                                                    #
#                                                                                         #
# The first stage is in common_computation.jl, the argument and data checks are in        #
# data_checks.jl and date_helpers.jl, and the engine that recomputes the estimates for    #
# the jackknife and the randomization inference is in batched_engine.jl.                  #
#                                                                                         #
# Papers referred to in the comments:                                                     #
#   Karim & Webb (2025), "Good Controls Gone Bad: Difference-in-Differences with          #
#       Covariates" (the DID-INT paper; equation numbers are from this paper)             #
#   Austin, Strumpf, Karim & Webb (2026), "Which Policy Works, and Where?"                #
#   Callaway & Sant'Anna (2021), "Difference-in-Differences with multiple time periods"   #
#   MacKinnon & Webb (2020), "Randomization inference for difference-in-differences       #
#       with few treated clusters"                                                        #
#   Racine & MacKinnon (2007), "Inference via kernel smoothing of bootstrap P values"     #
#   MacKinnon, Nielsen & Webb (2023), "Cluster-robust inference: A guide to empirical     #
#       practice" and "Fast and reliable jackknife and bootstrap methods for              #
#       cluster-robust inference"                                                         #
###########################################################################################


###########################################################################################
# Section 1: didint(), the main function                                                  #
###########################################################################################

"""
    didint(
           outcome::Union{AbstractString, Symbol},
           state::Union{AbstractString, Symbol},
           time::Union{AbstractString, Symbol},
           data::Union{DataFrame, RConnector.RDataFrame};
           gvar::Union{AbstractString, Symbol, Nothing} = nothing,
           treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
           treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
           pvar::Union{AbstractString, Symbol, Nothing} = nothing,
           treated_policies::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
           date_format::Union{AbstractString, Nothing} = nothing,
           covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing,
           ccc::AbstractString = "int",
           agg::AbstractString = "cohort",
           weighting::AbstractString = "both",
           policy_aggregation::Bool = false,
           policy_weighting::AbstractString = "equal",
           ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
           freq::Union{AbstractString, Nothing} = nothing,
           freq_multiplier::Number = 1,
           start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
           end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
           nperm::Number = 999,
           ri_kernel_cutoff::Union{Number, Nothing} = 99,
           seed::Number = rand(1:1000000),
           use_pre_controls::Bool = false,
           notyet::Union{Nothing, Bool} = nothing,
           hc::Union{AbstractString, Number} = "hc1",
           truejack::Bool = false,
           edgecase::Bool = false
          )

The `didint()` function estimates the average effect of treatment on the treated (ATT)
while accounting for covariates that may vary by state, time, or by both state and time simultaneously.

# Details
The estimator of Karim & Webb (2025) is computed in four steps:

1. **First stage.** A covariate-adjusted mean of the outcome, `λ(s,t)`, is estimated for every state `s`
   and period `t`, with the covariate effects restricted as the `ccc` argument specifies.
2. **Long differences.** For every state, the change in `λ` between a period `t` and the period just
   before a cohort's treatment time.
3. **ATT(s,t).** For each treated state `s` and each period `t` from its treatment time on, its long
   difference minus the weighted mean of the long differences of its controls: the never-treated states
   and, with `notyet = true`, the states not yet treated at `t`.
4. **Aggregation.** The sub-aggregate ATTs of the chosen `agg` are weighted means of the ATT(s,t)'s, and the
   aggregate ATT is a weighted mean of the sub-aggregate ATTs.

Common adoption (every treated state treated at the same time) is the one-cohort case of the same
estimator, so every `agg` option is available whatever the adoption pattern.

The arguments `treated_states` and `treatment_times` should be entered such that the first element in
`treated_states` refers to the state treated at the date entered as the first element in `treatment_times`,
the second element in `treated_states` refers to the state treated at the date entered as the second element
in `treatment_times`, and so on. A single treatment time applies to every treated state.

# Parameters

## Required Parameters
- `outcome::Union{AbstractString, Symbol}`
    Input the name of the column which identifies the outcome of interest.
- `state::Union{AbstractString, Symbol}`
    Input the name of the column which identifies the state membership of the observation.
- `time::Union{AbstractString, Symbol}`
    Input the name of the column which identifies the date of the observation.
- `data::Union{DataFrame, RConnector.RDataFrame}`
    The DataFrame to be used for the analysis.

## Treatment Specification
- `gvar::Union{AbstractString, Symbol, Nothing} = nothing`
    Name of the column which indicates time of first treatment for each state (missing or `0` for
    never-treated states). Use either `gvar`, or `treated_states` and `treatment_times`.
- `treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing`
    A vector of strings or numbers (or a single one) noting the treated state(s), matching the type of the `state` column.
- `treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing`
    A vector (or single entry) denoting the associated treatment times of the 'treated_states'.
- `pvar::Union{AbstractString, Symbol, Nothing} = nothing`
    Name of the column which indicates the policy each state adopts (constant within a state;
    missing, `0`, `"0"` or `""` for never-treated states). Use either `pvar` or `treated_policies`.
- `treated_policies::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing`
    A vector (or single entry) denoting the policies of the 'treated_states', in the same order. It cannot
    be combined with `gvar` (use `pvar` instead).
    If neither `pvar` nor `treated_policies` is given, every treated state has the same policy.
    Policies only matter for the policy aggregations: controls are never-treated states (and, with
    `notyet`, not-yet-treated states) whatever their policy.

## Model Specification
- `ccc::AbstractString = "int"`
    Specify which version of DID-INT should be used.
    Options are: `"hom"`, `"time"`, `"state"`, "`add`", and `"int"`.
    The long differences are the integrated ones of Karim & Webb (2025, eqs. 46-47): the change in a state's
    covariate slopes is evaluated at the treated state's covariate means. That covariate term cancels
    in every ATT(s,t) except under `"int"`, so only `"int"` estimates with covariates depend on it.
- `agg::AbstractString = "cohort"`
    Enter the aggregation method as a string. Each option groups the ATT(s,t)'s differently:
    `"cohort"` (by treatment time), `"simple"` (by treatment time and period), `"state"` (by treated state),
    `"sgt"` (no grouping: one ATT per treated state and period), `"time"` (by periods since treatment),
    `"calendar"` (by period), and the policy aggregations `"policy"`, `"policy_gt"`, `"policy_cohort"`,
    `"policy_calendar"` and `"policy_time"`.
    Every sub-aggregate is a weighted mean of ATT(s,t) estimates, and `agg_att` is a weighted mean of the
    sub-aggregates, so `agg_att` can differ across `agg` options (it does not under `weighting = "att"`
    or `"both"`).
    The policy aggregations form the same sub-aggregates within each policy (`"policy"` has one per
    policy) and report each policy's ATT in `att_p`, on the first row of that policy, instead of `agg_att`.
    Unlike the group effects of Callaway & Sant'Anna (2021), which weight a cohort's post-treatment periods
    equally, `"cohort"` weights them as `weighting` specifies.
- `weighting::AbstractString = "both"`
    Specify which weighting method should be used. Observation counts can be used as weights when
    computing the ATT(s,t)'s and sub-aggregates from long differences (`"diff"`), when computing the
    aggregate ATT from the sub-aggregates (`"att"`), at both levels (`"both"`), or not at all (`"none"`).
    Options are: `"both"`, `"att"`, `"diff"`, or `"none"`.
- `policy_aggregation::Bool = false`
    With a policy aggregation, also report the pooled estimate across policies, `att_pool`, the
    `policy_weighting`-weighted mean of the policies' ATTs, with its own inference columns. Its HC standard
    error has (number of policies - 1) degrees of freedom, so two policies give one degree of freedom.
    `att_pool` is the overall ATT across policies for the policy aggregations (they have no `agg_att`).
    With a non-policy aggregation there is nothing to pool and the argument is ignored.
- `policy_weighting::AbstractString = "equal"`
    Weight of each policy in `att_pool`: `"equal"` (the typical policy), `"treated_obs"` (the policy's
    treated observations), `"treated_cells"` (its number of treated state-periods) or `"treated_states"`
    (its number of treated states). Only used with `policy_aggregation = true`.
- `covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing`
    A vector of covariates entered as either strings or symbols (or a single covariate string or symbol),
    or, `nothing` (default). The order of the covariates matters when covariates are collinear: under
    `ccc = "hom"`, `"state"`, `"time"` and `"int"`, a covariate that is a linear combination of covariates
    listed before it is dropped, with a warning, wherever the covariate slopes are estimated (in all the data,
    in each state, in each period, or in each state-period cell).
- `notyet::Union{Nothing, Bool} = nothing`
    Determine if pre-treatment periods from treated states should be used as controls
    (`use_pre_controls` is an alias).
- `ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing`
    A dictionary specifying which category in a categorical variable should be used
    as the reference (baseline) category.

## Date Processing & Period Grid Construction
- `date_format::Union{AbstractString, Nothing} = nothing`
    Date format (e.g. "yyyy" or "yyyy-mm-dd") to be used when parsing string dates from
    the time column, or `start_date`, `end_date`, and `treatment_times` arguments.
- `freq::Union{AbstractString, Nothing} = nothing`
    A string indicating the desired timeframe of a period for the analysis.
    Options are: `"year"`, `"month"`, `"week"`, `"day"`.
- `freq_multiplier::Number = 1`
    An integer by which the 'freq' argument should be multiplied, e.g. if a two-year
    period is desired, set `freq = "year"` and `freq_multiplier = 2`.
- `start_date::Union{AbstractString, Number, Date, Nothing} = nothing`
    Any data prior this date is dropped, and serves as the starting date for the period
    grid construction if activated.
- `end_date::Union{AbstractString, Number, Date, Nothing} = nothing`
    Any data after this date is dropped, and serves as the end date for the period
    grid construction if activated.

## Inference
- `nperm::Number = 999`
    The number of unique permutations (not including the initial assignment of treatment times)
    to be considered when performing the randomization inference. If no more than `nperm` unique
    permutations exist, all of them are used.
- `ri_kernel_cutoff::Union{Number, Nothing} = 99`
    When the randomization inference uses `ri_kernel_cutoff` or fewer permutations (because
    `nperm` is that small, or because no more unique permutations exist), its p-values are kernel-smoothed
    (Racine & MacKinnon 2007) instead of counting the placebo estimates at least as large as the observed one.
    Set it to `nothing`, `false` or `0` to never smooth.
- `seed::Number = rand(1:1000000)`
    An integer to set the random seed for the randomization inference procedure. The same seed gives the
    same permutations. When a seed is given, the global random number generator is left untouched (the default seed is drawn from it).
- `hc::Union{AbstractString, Number} = "hc1"`
    Specify which heteroskedasticity-consistent covariance matrix estimator (HCCME) should be used.
    Options are `0`, `1`, `2`, `3`, and `4` (or `"hc0"`, `"hc1"`, `"hc2"`, `"hc3"`, `"hc4"`).
    An ATT(s,t) is estimated from one treated state and its controls, and the treated state's leverage is 1,
    so its HC2, HC3 and HC4 standard errors are `NaN`. The HC standard error of a sub-aggregate or aggregate
    comes from the spread of the estimates it averages and treats them as independent; estimates that
    share control states are not, so it can be too small. The jackknife and the randomization inference
    do not rely on that independence.
- `truejack::Bool = false`
    When `ccc` is set to either `"add"`, `"time"`, or `"hom"`, then in order to get valid jackknife
    standard errors, we need to re-estimate the DID-INT model from square one (re-running the first stage).
    This is because the covariate effects in those cases depend on values from across states, so dropping a state will change
    the lambda values, this is not true for the `ccc` options of `"int"` or `"state"`.
- `edgecase::Bool = false`
    When an ATT(s,t) has at least two control states, its standard error can be calculated
    directly from the regression of differences on the intercept term and treatment indicator. However, when there is only one control
    state for a particular ATT(s,t), the standard error can still be constructed from the variance and covariance terms of the relevant
    means. Calculating these `edgecase` standard errors can be computationally expensive, therefore by default this argument is set to `false`, but
    may be toggled on by setting `edgecase = true`.
    Under `ccc = "int"` with covariates the ATT also depends on the treated state's covariate means, and the `edgecase` standard error
    is conditional on the observed covariates: it treats those means as fixed. Every `edgecase` standard error assumes that observations in different
    state-period cells are independent, as in repeated cross-sections. With panel data, where the same individuals are observed in several
    periods, it omits the covariances between a state's periods, and the `edgecase` standard errors are then generally biased.
    Not available for `ccc = "add"`.

# Returns
A [`DiDIntResult`](@ref). Displaying it prints a summary of the specification, the aggregate ATT and the sub-aggregate ATTs.
Its fields hold the results as tables:

- `agg`: the aggregate ATT (for the policy aggregations, the ATT of each policy and, with `policy_aggregation = true`, the pooled ATT).
- `sub`: the sub-aggregate ATTs of the chosen `agg`, with the weight of each in the aggregate ATT (`missing` when they are weighted equally).
- `att_st`: the ATT of each treated state in each period from its treatment time on.
- `specs`: the specification that was estimated, including the period length, the date range and the number of permutations used.
- `table`: every result in one wide table (also returned by `DataFrame(result)`).

In `agg` and `sub` the estimates are in the column `att`, with `se` and `pval` (HC standard error and p-value), `jknife_se` and
`jknife_pval` (jackknife standard error and p-value) and `ri_pval` (randomization inference p-value).
`coef(result)`, `stderror(result)` and `coefnames(result)` return the estimates, standard errors and labels as vectors.

The wide table has one row per sub-aggregate ATT:

- The columns identifying the sub-aggregate: `treatment_time` (`"cohort"`), `r1`, `time` and `gvar` (`"simple"`:
  the period before treatment, the period, and the treatment time), `state` (`"state"`), `state`, `gvar` and `t` (`"sgt"`),
  `periods_post_treat` (`"time"`), `t` (`"calendar"`). The policy aggregations add a `policy` column.
- The sub-aggregate ATT, in `att_cohort`, `att_gt`, `att_s`, `att_sgt`, `att_t`, `att_calendar` or `att_p`, and its inference in the
  columns with the same name prefixed by `se_` and `pval_` (HC standard error and p-value), `jknifese_` and `jknifepval_`
  (jackknife standard error and p-value) and `ri_pval_` (randomization inference p-value).
- The aggregate ATT, `agg_att`, with `se_agg_att`, `pval_agg_att`, `jknifese_agg_att`, `jknifepval_agg_att` and `ri_pval_agg_att`,
  on the first row. For the policy aggregations each policy's ATT is in `att_p` (with the same prefixes) on the first row of
  that policy and, with `policy_aggregation = true`, the pooled ATT is in `att_pool` (with the same prefixes) on the first row and
  each policy's weight in it is in `policy_weights`.
- `nperm`: the number of permutations used in the randomization inference (first row).
- `weights`: the weight of each sub-aggregate in the aggregate ATT (for the policy aggregations, in its policy's ATT). It is empty
  (`missing`; `nothing` in the table passed to the wrappers) when the sub-aggregates are weighted equally, under `weighting = "none"` or `"diff"`.
- `period`, `start_date` and `end_date`: the period length and the date range of the data used.

A jackknife needs at least two treated and two control states, overall and for each estimate, and the sub-aggregates of
`"state"` and `"sgt"`, which each have a single treated state, never have jackknife standard errors. Estimates that cannot be
computed are `missing`.

# Citations
- Karim & Webb (2025). "Good Controls Gone Bad: Difference-in-Differences with Covariates". https://arxiv.org/abs/2412.14447
- MacKinnon & Webb (2020). "Randomization inference for difference-in-differences with few treated clusters". https://doi.org/10.1016/j.jeconom.2020.04.024
- Racine & MacKinnon (2007). "Inference via kernel smoothing of bootstrap P values". https://doi.org/10.1016/j.csda.2006.11.013
- Callaway & Sant'Anna (2021). "Difference-in-Differences with multiple time periods". https://doi.org/10.1016/j.jeconom.2020.12.001

"""
function didint(outcome::Union{AbstractString, Symbol},
                state::Union{AbstractString, Symbol},
                time::Union{AbstractString, Symbol},
                data;
                gvar::Union{AbstractString, Symbol, Nothing} = nothing,
                treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
                treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
                pvar::Union{AbstractString, Symbol, Nothing} = nothing,
                treated_policies::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
                date_format::Union{AbstractString, Nothing} = nothing,
                covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing,
                ccc::AbstractString = "int",
                agg::AbstractString = "cohort",
                weighting::AbstractString = "both",
                policy_aggregation::Bool = false,
                policy_weighting::AbstractString = "equal",
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                nperm::Number = 999,
                ri_kernel_cutoff::Union{Number, Nothing} = 99,
                seed::Number = rand(1:1000000),
                use_pre_controls::Bool = false,
                notyet::Union{Nothing, Bool} = nothing,
                hc::Union{AbstractString, Number} = "hc1",
                wrapper::Union{AbstractString, Nothing} = nothing,
                truejack::Bool = false,
                edgecase::Bool = false
                )

    #--------------------------------------------------------------------------------------
    # Step 1: Check the arguments
    #--------------------------------------------------------------------------------------

    hc = hc_checks(hc)
    wrapper = init_wrapper_check(wrapper)
    data_type_check(data)
    date_format = check_date_format(date_format)
    freq_multiplier = check_positive_integer(freq_multiplier, "freq_multiplier")

    # Let notyet override use_pre_controls
    use_pre_controls = isnothing(notyet) ? use_pre_controls : notyet

    # Column names as strings, whether they were given as strings or symbols
    outcome, state, time, gvar = symbol_to_string(outcome, state, time, gvar)
    covariates = isnothing(covariates) ? covariates : string.(covariates)

    # Treated states and their treatment times, from gvar or from treated_states and treatment_times
    treated_states, treatment_times = gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                                              event = false, check = "estimate_didint")

    # Policy of each treated state, from pvar or treated_policies
    state_policies = policy_assignments(pvar, treated_policies, treated_states, gvar, data, state)

    start_date = start_end_date_checks(start_date, "start_date", date_format)
    end_date = start_end_date_checks(end_date, "end_date", date_format)

    ccc = check_option(ccc, "ccc", ["hom", "time", "state", "add", "int"])
    seed = check_positive_integer(seed, "seed")
    nperm = check_positive_integer(nperm, "nperm")
    ri_kernel_cutoff = check_ri_kernel_cutoff(ri_kernel_cutoff)

    agg = lowercase(agg)
    if agg == "none"
        error("agg = \"none\" has been removed; every agg reports the aggregate ATT in agg_att. Use e.g. agg = \"simple\".")
    end
    agg = check_option(agg, "agg", ["cohort", "state", "simple", "sgt", "time", "calendar",
                                    "policy", "policy_gt", "policy_cohort", "policy_calendar", "policy_time"])
    weighting = check_option(weighting, "weighting", ["both", "none", "att", "diff"])

    # Pooling across policies only exists for the policy aggregations; with any other agg,
    # policy_aggregation has nothing to pool and is ignored
    policy_weighting = check_option(policy_weighting, "policy_weighting", ["equal", "treated_obs", "treated_cells", "treated_states"])
    if !policy_aggregation && policy_weighting != "equal"
        error("'policy_weighting' only applies when 'policy_aggregation = true'.")
    end
    if !is_policy_agg(agg)
        policy_aggregation = false
    end

    #--------------------------------------------------------------------------------------
    # Step 2: Clean the data and match its dates to periods
    #--------------------------------------------------------------------------------------

    # Work on a DataFrame that shares the columns of `data`; columns added here are not added to `data`
    data_copy = DataFrame(data; copycols = false)

    # Check that the specified columns exist and can be processed
    data_copy, treatment_times, covariates = validate_data(data_copy, outcome, state, time, covariates, treatment_times, date_format;
                                               warn_missing_covariates = false)

    # Treated states and treatment times as vectors of the same length, each state listed once.
    # A single treatment time applies to every treated state: common adoption is the one-cohort
    # case of the same estimator.
    treated_states, treatment_times = check_treated_states(treated_states, treatment_times)

    # Make sure control states exist. Not-yet-treated states can only serve as controls when
    # there are at least two treatment times.
    control_states = setdiff(unique(data_copy[!, state]), treated_states)
    if isempty(control_states) && !use_pre_controls
        error("No control states were found.")
    elseif isempty(control_states) && length(unique(treatment_times)) == 1
        error("No valid control states were found.")
    end

    # Drop rows with missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # The state column and treated_states must both be numbers or both be strings
    validate_state_types(treated_states, data_copy, state)
    data_copy.state_71X9yTx = data_copy[!, state]
    if (length(unique(data_copy.state_71X9yTx)) == 1)
        error("Not enough non-missing observations: found only 1 state after dropping missing observations.")
    end

    # Convert dates to Date objects and keep the data between start_date and end_date
    data_copy, treatment_times, start_date, end_date, all_times, freq = validate_and_convert_dates(data_copy, time, treatment_times, date_format,
                                                                                                   freq, start_date, end_date)

    # Build the period grid and match every date and treatment time to a period on it
    data_copy, treatment_times, treated_states, match_to_these_dates, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier,
                                                                                                                    start_date, end_date, treatment_times,
                                                                                                                    treated_states; drop_already_treated = true)
    if isempty(treated_states)
        error("None of the treated states has a treatment time within the data used (after its first date and no later than its last date).")
    end

    # Check that the treated states exist, and that each has data before and from its treatment time
    treated_states, treatment_times = validate_treated_states(treated_states, treatment_times, data_copy)
    check_treatment_windows(treated_states, treatment_times, data_copy)

    # States as strings from here on (as are treated_states)
    data_copy, treated_states = validate_string_treated_states(data_copy; treated_states = treated_states)

    # Keep the period of each observation as a Date (time_dmG5fpM) and as a string (time_71X9yTx)
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)

    # Convert factor covariates into numeric dummy variable columns
    data_copy, covariates_to_include = process_covariates(covariates, data_copy, ref)

    # A Float64 outcome column makes the regressions faster than a Union{Number, Missing} one
    data_copy.outcome_71X9yTx = convert(Vector{Float64}, data_copy.outcome_71X9yTx)

    # The edgecase standard errors are not available for ccc = "add" (its first stage is rank deficient)
    if edgecase == true && ccc == "add"
        @warn "'edgecase' standard errors are not available for the 'ccc = add' option."
        edgecase = false
    end

    # Under "int" and "state" deleting a state does not change the other states' estimates, so
    # the jackknife never needs to re-estimate the first stage
    if ccc in ["int", "state"]
        truejack = false
    end

    #--------------------------------------------------------------------------------------
    # Step 3: First stage. lambda_{s,t}, the covariate-adjusted mean outcome of each state
    # and period (Karim & Webb 2025, eq. 45)
    #--------------------------------------------------------------------------------------

    lambda_df, vcov_lambda, covariate_fit = first_stage(data_copy, ccc, covariates_to_include, hc, edgecase,
                                                        time_to_index = time_to_index, treatment_times = treatment_times,
                                                        match_to_these_dates = match_to_these_dates,
                                                        treated_states = treated_states, use_pre_controls = use_pre_controls)

    # Integrated long differences (Karim & Webb 2025, eqs. 46-47) add a covariate term to each
    # long difference. It cancels in ATT(s,t) for every ccc except "int": it is zero for "hom"
    # and "state", and for "time" and "add" the treated state and its controls carry the same
    # time slopes at the same covariate means. So it is only applied under "int" with covariates.
    integration = (ccc == "int" && !isempty(covariates_to_include)) ? covariate_fit : nothing

    #--------------------------------------------------------------------------------------
    # Step 4: Long differences of the lambdas
    #--------------------------------------------------------------------------------------

    # diff_df holds the long differences of the treated states and of their controls;
    # ri_diff_df the extra ones that randomization inference needs
    diff_df, ri_diff_df = construct_diff_df(lambda_df, treated_states, treatment_times, time_to_index, match_to_these_dates,
                                             use_pre_controls, data_copy)

    # Randomization inference is only possible if every state has every long difference
    randomize = check_prior_to_ri(diff_df, ri_diff_df)

    #--------------------------------------------------------------------------------------
    # Step 5: ATT(s,t) for every treated state and post-treatment period (Karim & Webb 2025,
    # eq. 48), then the sub-aggregates of the chosen agg and the aggregate ATT, all weighted
    # means of ATT(s,t)'s (eq. 49)
    #--------------------------------------------------------------------------------------

    # Policies of the treated states. They only enter the policy aggregations: controls do not
    # depend on policy, so every other agg treats all treated states as having one policy.
    policy = (labels = Dict(s => (is_policy_agg(agg) ? state_policies[s] : DEFAULT_POLICY) for s in treated_states),
              pooled = policy_aggregation, weighting = policy_weighting)

    # Warn where the data leave an ATT(s,t) undetermined: under "int", the ATT(s,t)'s that
    # depend on covariates dropped from a cell's regression; under "add", rank deficiency
    # beyond the split of the slopes between states and periods
    isnothing(integration) || warn_dropped_covariates_int(integration, diff_df)
    ccc == "add" && warn_additive_identification(covariate_fit, diff_df, weighting)

    att_st = build_att_st(diff_df, weighting, hc, vcov_lambda; integration = integration, policies = policy.labels)
    results = sub_aggregate_results(att_st, agg, weighting, hc; pooled = policy.pooled)
    if is_policy_agg(agg)
        # Each policy's ATT and, with policy_aggregation, the pooled ATT across policies
        results = policy_att_results(results, att_st, agg, weighting, hc, policy)
    else
        results = aggregate_att(results, agg, weighting, hc)
    end

    #--------------------------------------------------------------------------------------
    # Step 6: Jackknife and randomization inference. Both recompute the same weighted means,
    # on each deletion of a state and on each permutation of the treatment assignment.
    #--------------------------------------------------------------------------------------

    results = jackknife_procedure(results, att_st, diff_df, agg, weighting, use_pre_controls,
                                  truejack, data_copy, ccc, covariates_to_include, treated_states, treatment_times,
                                  time_to_index, match_to_these_dates; integration = integration, policy = policy)
    if randomize
        results = randomization_inference(results, att_st, vcat(diff_df, ri_diff_df), agg, weighting,
                                          use_pre_controls, treated_states, treatment_times, nperm, ri_kernel_cutoff, seed;
                                          integration = integration, policy = policy)
    end

    #--------------------------------------------------------------------------------------
    # Step 7: Add the period length and the date range to the results
    #--------------------------------------------------------------------------------------

    if isnothing(date_format) || date_format == "yyyy"
        date_format = assume_date_format(period)
    end
    results.period .= string(period)
    results.start_date .= parse_date_to_string_didint(start_date, date_format)
    results.end_date .= parse_date_to_string_didint(maximum(data_copy.time_dmG5fpM), date_format)

    #--------------------------------------------------------------------------------------
    # Step 8: Return the results. The Stata and R wrappers receive the table of results as
    # it is; otherwise it is arranged into a DiDIntResult (see results.jl).
    #--------------------------------------------------------------------------------------

    if wrapper in ("stata", "r")
        return wrapper_check(results, wrapper)
    end
    specs = (outcome = outcome, covariates = isnothing(covariates) ? String[] : String.(covariates),
             ccc = ccc, agg = agg, weighting = weighting,
             notyet = use_pre_controls, hc = hc, truejack = truejack, edgecase = edgecase,
             policy_aggregation = policy_aggregation, policy_weighting = policy_weighting,
             nperm = randomize ? Int(results.nperm[1]) : missing, ri_kernel_cutoff = ri_kernel_cutoff, seed = seed,
             period = results.period[1], start_date = results.start_date[1], end_date = results.end_date[1])
    return DiDIntResult(results, att_st, specs)
end

###########################################################################################
# Section 2: Long differences                                                             #
#                                                                                         #
# A long difference is the change in a state's lambda between a period t and r1, the      #
# period just before a cohort's treatment time g (Karim & Webb 2025, eq. 46). Each        #
# (t, r1) pair is a "cell" of cohort g. The treated states of cohort g and their          #
# controls each have one long difference per cell.                                        #
###########################################################################################

function long_difference_rows(state::AbstractString, t::AbstractVector, r1::AbstractVector, treated_time::AbstractVector,
                              treat::Int, lambdas::AbstractDict, n_obs::AbstractDict)

    # Long differences of one state in the cells (t[j], r1[j]) of the cohorts treated_time[j]:
    #   diff        lambda_{s,t} - lambda_{s,r1} (missing if the state is not observed in t or r1)
    #   n           observations behind the long difference, n_{s,t} + n_{s,r1}
    #   n_t         observations in period t, n_{s,t}
    #   lambda_idx  positions of (lambda_{s,r1}, lambda_{s,t}) in the covariance matrix of the
    #               lambdas, always in that order (edge_case_se relies on it)
    # `treat` labels the rows: 1 treated, 0 control, -1 only used by randomization inference.
    diffs = Vector{Union{Float64, Missing}}(missing, length(t))
    lambda_indices = Vector{Union{Tuple{Int, Int}, Missing}}(missing, length(t))
    n = Vector{Int}(undef, length(t))
    n_t = Vector{Int}(undef, length(t))
    for j in eachindex(t)
        lambda_t = get(lambdas, (state, t[j]), nothing)
        lambda_r1 = get(lambdas, (state, r1[j]), nothing)
        if !isnothing(lambda_t) && !isnothing(lambda_r1)
            diffs[j] = lambda_t.lambda - lambda_r1.lambda
            lambda_indices[j] = (lambda_r1.index, lambda_t.index)
        end
        n_t[j] = get(n_obs, (state, t[j]), 0)
        n[j] = n_t[j] + get(n_obs, (state, r1[j]), 0)
    end
    return DataFrame(state = state, treated_time = treated_time, t = t, r1 = r1, diff = diffs, treat = treat,
                     n = n, n_t = n_t, lambda_idx = lambda_indices)
end

function construct_diff_df(lambda_df, treated_states, treatment_times, time_to_index, match_to_these_dates, use_pre_controls,
                           data_copy; jackknife = false, without = nothing)

    # Build the long differences behind every ATT(s,t). Returns
    #   diff_df     the long differences of the treated states in their own cohort's cells
    #               (treat = 1) and those of their controls in the same cells (treat = 0)
    #   ri_diff_df  the long differences of the treated states in the cells of the other cohorts
    #               (treat = -1), which randomization inference needs to give a state another
    #               state's treatment time
    # With `jackknife = true` the state `without` has been deleted: only diff_df is returned.

    if jackknife
        mask = treated_states .!= without
        treated_states = treated_states[mask]
        treatment_times = treatment_times[mask]
    end

    # Lookups: each state-period's lambda with its position in the covariance matrix of the
    # lambdas, and its number of observations
    lambdas = Dict((lambda_df.state[i], lambda_df.time[i]) => (lambda = Float64(lambda_df.lambda[i]), index = Int(lambda_df.lambda_index[i]))
                   for i in 1:nrow(lambda_df))
    cell_counts = combine(groupby(data_copy, [:state_71X9yTx, :time_dmG5fpM]), nrow => :n)
    n_obs = Dict(zip(zip(cell_counts.state_71X9yTx, cell_counts.time_dmG5fpM), cell_counts.n))

    empty_diff_df() = DataFrame(state = String[], treated_time = Date[], t = Date[], r1 = Date[], diff = Float64[],
                                treat = Int[], n = Int[], n_t = Int[], lambda_idx = Union{Tuple{Int, Int}, Missing}[])
    periods_since_treatment(df) = [time_to_index[t] - time_to_index[g] for (t, g) in zip(df.t, df.treated_time)]

    # Treated states: one long difference for every period after r1 in which the state is observed
    pieces = DataFrame[]
    for (s, g) in zip(treated_states, treatment_times)
        r1 = match_to_these_dates[time_to_index[g] - 1]
        t = sort(lambda_df.time[(lambda_df.state .== s) .& (lambda_df.time .> r1)])
        push!(pieces, long_difference_rows(s, t, fill(r1, length(t)), fill(g, length(t)), 1, lambdas, n_obs))
    end
    diff_df = vcat(empty_diff_df(), pieces...)

    # Control states: the same cells as the treated states
    cells = unique(select(diff_df, :t, :r1, :treated_time))
    unique_states = unique(lambda_df.state)
    control_states = setdiff(unique_states, treated_states)
    pieces = [long_difference_rows(c, cells.t, cells.r1, cells.treated_time, 0, lambdas, n_obs) for c in control_states]
    diff_df = vcat(diff_df, pieces...)
    diff_df.time_since_treatment = periods_since_treatment(diff_df)

    # Treated states in the cells of the other cohorts, for randomization inference (and, with
    # use_pre_controls, as not-yet-treated controls)
    if !jackknife || use_pre_controls
        pieces = DataFrame[]
        for (s, g) in zip(treated_states, treatment_times)
            other_cells = unique(select(diff_df[(diff_df.state .!= s) .& (diff_df.treated_time .!= g), :], :t, :r1, :treated_time))
            push!(pieces, long_difference_rows(s, other_cells.t, other_cells.r1, other_cells.treated_time, -1, lambdas, n_obs))
        end
        ri_diff_df = vcat(empty_diff_df(), pieces...)
        ri_diff_df.time_since_treatment = periods_since_treatment(ri_diff_df)
    end

    # With use_pre_controls, a treated state is a control in the cells of other cohorts whose
    # period t comes before its own treatment time: move those rows to diff_df
    if use_pre_controls
        for (t,s) in zip(treatment_times, treated_states)
            ri_diff_df[ri_diff_df.state .== s .&& ri_diff_df.t .< t, "treat"] .= 0
        end
        diff_df = vcat(diff_df, ri_diff_df[ri_diff_df.treat .== 0,:])
        ri_diff_df = ri_diff_df[ri_diff_df.treat .== -1,:]
        unique_states = unique(diff_df.state)
        control_states = setdiff(unique_states, treated_states)
        if isempty(control_states)
            # Without never-treated states, drop the cells that have no treated-control pair
            diff_df = subset(
                        groupby(diff_df, [:t, :r1]),
                        :treat => x -> any(x .== 1) && any(x .== 0)
                      )
            ri_diff_df = semijoin(ri_diff_df, diff_df, on=[:t, :r1])
        end
    end

    # Drop missing long differences, and then the cells left without a treated-control pair
    if any(ismissing.(diff_df.diff))
        diff_df = dropmissing(diff_df, :diff)
        diff_df = subset(
                        groupby(diff_df, [:t, :r1]),
                        :treat => x -> any(x .== 1) && any(x .== 0)
                      )
        if nrow(diff_df) == 0 && !jackknife
            error("Could not find a single pair of non-missing valued treated & untreated differences for the same (g,t) group.\nUnable to compute any ATTs.")
        end
        if !jackknife
            ri_diff_df = semijoin(ri_diff_df, diff_df, on=[:t, :r1])
            ri_diff_df = dropmissing(ri_diff_df, :diff)
        end
    end

    if !jackknife
        return diff_df, ri_diff_df
    end
    return diff_df

end

function check_prior_to_ri(diff_df, ri_diff_df)

    # Randomization inference gives every state every cohort's treatment time in turn, so it
    # needs each state's long difference in every cell. Returns false, with a warning, if a
    # state lacks one.
    combined_df = vcat(diff_df, ri_diff_df)
    required_cells = Set(zip(combined_df.t, combined_df.r1))
    cells_of_state = Dict{String, Set{Tuple{Date, Date}}}(s => Set{Tuple{Date, Date}}() for s in unique(combined_df.state))
    for (s, t, r1, diff) in zip(combined_df.state, combined_df.t, combined_df.r1, combined_df.diff)
        ismissing(diff) || push!(cells_of_state[s], (t, r1))
    end
    randomize = all(length(cells) == length(required_cells) for cells in values(cells_of_state))
    if !randomize
        @warn "Randomization inference cannot be performed.\nDid not find all the non-missing values for differences required to implement randomization inference procedure."
    end
    return randomize
end

function long_difference_lookup(long_diffs::AbstractDataFrame)

    # (state, t, r1) => (diff, n, n_t) for every non-missing long difference
    lookup = Dict{Tuple{String, Date, Date}, NamedTuple{(:diff, :n, :n_t), Tuple{Float64, Float64, Float64}}}()
    for row in eachrow(long_diffs)
        ismissing(row.diff) && continue
        lookup[(row.state, row.t, row.r1)] = (diff = Float64(row.diff), n = Float64(row.n), n_t = Float64(row.n_t))
    end
    return lookup
end

###########################################################################################
# Section 3: ATT(s,t), sub-aggregates and the aggregate ATT                               #
#                                                                                         #
# ATT(s,t) is estimated for every treated state s and post-treatment period t (Karim &    #
# Webb 2025, eq. 48). Every sub-aggregate is a weighted mean of ATT(s,t)'s and the        #
# aggregate ATT a weighted mean of sub-aggregates (eq. 49; Callaway & Sant'Anna 2021),    #
# so no regression ever pools long differences from different cells.                      #
###########################################################################################

# For each aggregation: the suffix of its sub-aggregate columns (att_cohort, se_att_cohort, ...)
# and its key columns as (results column, ATT(s,t) table column, element type). The key columns
# define the groups of ATT(s,t)'s that are averaged. The results columns are the public output
# names that the Stata and R wrappers rely on.
const SUB_AGG_SPECS = Dict(
    "cohort"   => (suffix = "att_cohort", keys = [(:treatment_time, :gvar, Date)]),
    "simple"   => (suffix = "att_gt", keys = [(:r1, :r1, Date), (:time, :t, Date), (:gvar, :gvar, Date)]),
    "state"    => (suffix = "att_s", keys = [(:state, :state, String)]),
    "sgt"      => (suffix = "att_sgt", keys = [(:state, :state, String), (:gvar, :gvar, Date), (:t, :t, Date)]),
    "time"     => (suffix = "att_t", keys = [(:periods_post_treat, :periods_post_treat, Float64)]),
    "calendar" => (suffix = "att_calendar", keys = [(:t, :t, Date)]),
    # Policy aggregations (Austin, Strumpf, Karim & Webb 2026): the same groups within each
    # policy. "policy" has one sub-aggregate per policy, the policy's ATT itself.
    "policy"          => (suffix = "att_p", keys = [(:policy, :policy, String)]),
    "policy_gt"       => (suffix = "att_gt", keys = [(:policy, :policy, String), (:r1, :r1, Date), (:time, :t, Date), (:gvar, :gvar, Date)]),
    "policy_cohort"   => (suffix = "att_cohort", keys = [(:policy, :policy, String), (:treatment_time, :gvar, Date)]),
    "policy_calendar" => (suffix = "att_calendar", keys = [(:policy, :policy, String), (:t, :t, Date)]),
    "policy_time"     => (suffix = "att_t", keys = [(:policy, :policy, String), (:periods_post_treat, :periods_post_treat, Float64)]),
)

is_policy_agg(agg::AbstractString) = startswith(agg, "policy")

# Columns of the ATT(s,t) table that define the groups of an aggregation
group_columns(agg::AbstractString) = [att_col for (_, att_col, _) in SUB_AGG_SPECS[agg].keys]

# Inference columns reported for an estimate column, e.g. "att_p" -> "se_att_p", ...
inference_columns(estimate::AbstractString) = [prefix * estimate for prefix in ["se_", "pval_", "jknifese_", "jknifepval_", "ri_pval_"]]

# Policy options when none are given: one common policy, no pooling
default_policy_options(treated_states) = (labels = Dict(s => DEFAULT_POLICY for s in treated_states),
                                          pooled = false, weighting = "equal")

# The (cohort, period) cells of the observed ATT(s,t)'s, with their base period and event time
att_st_cells(att_st::DataFrame) = sort(unique(att_st[:, [:gvar, :t, :r1, :periods_post_treat]]), [:gvar, :t])

function result_key(results::DataFrame, i::Integer, agg::AbstractString)

    # Key tuple of results row i, with the event time as an Int as in the ATT(s,t) table
    return Tuple(coltype === Float64 ? Int(results[i, col]) : results[i, col]
                 for (col, _, coltype) in SUB_AGG_SPECS[agg].keys)
end

function init_results(agg::AbstractString, nrows::Integer; pooled::Bool = false)

    # Results DataFrame with the public column names in their order: key columns, sub-aggregate
    # ATT, top-level ATT and its inference, sub-aggregate inference, weights.
    # The top-level ATT is the aggregate ATT (agg_att) for the non-policy aggregations. For the
    # policy aggregations it is each policy's ATT (att_p, on the first row of the policy; for
    # agg = "policy" it is the sub-aggregate itself) and, with `pooled`, the pooled ATT across
    # policies (att_pool) with each policy's weight in it (policy_weights).
    float_column() = Vector{Union{Missing, Float64}}(missing, nrows)
    spec = SUB_AGG_SPECS[agg]
    top_level = String[]
    if !is_policy_agg(agg)
        top_level = ["agg_att", "se_agg_att", "pval_agg_att", "jknifese_agg_att", "jknifepval_agg_att", "ri_pval_agg_att"]
    else
        agg != "policy" && append!(top_level, vcat("att_p", inference_columns("att_p")))
        pooled && append!(top_level, vcat("att_pool", inference_columns("att_pool")))
    end

    results = DataFrame()
    for (col, _, coltype) in spec.keys
        results[!, col] = coltype === Float64 ? float_column() : Vector{coltype}(undef, nrows)
    end
    results[!, spec.suffix] = float_column()
    for col in vcat(top_level, "nperm", inference_columns(spec.suffix))
        results[!, col] = float_column()
    end
    results.weights = Vector{Union{Nothing, Float64}}(nothing, nrows)
    if pooled
        results.policy_weights = float_column()
    end
    return results
end

function control_weights(long_diffs::DataFrame, weighting::AbstractString)

    # Weights of the long differences in an ATT(s,t) regression: their observation counts n
    # under weighting = "diff" or "both", otherwise unweighted. The treated state's own weight
    # cancels, so these only weight the controls against each other.
    if in(weighting, ["diff", "both"])
        W = convert(Vector{Float64}, long_diffs.n)
        W ./= sum(W)
        return W
    else
        return fill(nothing, nrow(long_diffs))
    end
end

function slope_change_term(fit::CovariateFit, source_state, t, r1, target_state)

    # Covariate term of an integrated long difference (Karim & Webb 2025, eqs. 46-47):
    # (gamma_{source,t} - gamma_{source,r1})' xbar_{target,t}, the change in the source state's
    # slopes evaluated at the target (treated) state's period-t covariate means. Because f(X) is
    # linear, integrating over the target's covariate distribution reduces to its mean.
    return dot(fit.slopes[(source_state, t)] .- fit.slopes[(source_state, r1)], fit.xbar[(target_state, t)])
end

function dropped_covariate_effects(fit::CovariateFit, diff_df::DataFrame)

    # ccc = "int": the ATT(s,t)'s that depend on which covariates a cell's regression drops.
    # A cell k that drops covariates (they do not vary there, or are collinear with covariates
    # listed before them) has slopes that its data do not determine: any multiple of a vector
    # n in the null space of its centered covariates can be added to them without changing its
    # fit. Its fitted value at covariate values x, ybar_k + gamma_k' (x - xbar_k), is the same
    # for all of those slopes if and only if n' (x - xbar_k) = 0 for every such n.
    # In ATT(s,t) the cells (s,r1), (c,t) and (c,r1) of the treated state s and of each of its
    # controls c are evaluated at x = xbar_{s,t} (Karim & Webb 2025, eqs. 46-47). The treated
    # cell (s,t) is evaluated at its own means and is never affected.
    # The condition is tested with a relative tolerance: |n' (xbar_{s,t} - xbar_k)| must not
    # exceed IDENTIFICATION_RTOL times sum_i |n_i| scale_i, where scale_i is the root mean
    # square of covariate i over the data.
    # Returns one entry per (ATT(s,t), cell) pair that fails it, with the covariates dropped in
    # that cell that are responsible.
    effects = NamedTuple{(:state, :t, :cell, :covariates), Tuple{String, Date, Tuple{String, Date}, Vector{String}}}[]
    isempty(fit.null_space) && return effects

    # For every cell that drops covariates and every vector n of its basis: n' xbar_k, the
    # tolerance, and the dropped covariate that the vector belongs to (its last non-zero entry)
    checks = Dict(cell => (basis = basis, own = basis' * fit.xbar[cell],
                           tolerance = IDENTIFICATION_RTOL .* (abs.(basis)' * fit.covariate_scale),
                           dropped = [fit.covariates[findlast(!iszero, n)] for n in eachcol(basis)])
                  for (cell, basis) in fit.null_space)

    # Add an entry if the fitted value of `cell` at the covariate means `target` of ATT(s,t)
    # depends on the cell's dropped covariates
    function check!(s, t, target, cell)
        check = get(checks, cell, nothing)
        isnothing(check) && return nothing
        fails = [abs(dot(view(check.basis, :, j), target) - check.own[j]) > check.tolerance[j] for j in eachindex(check.own)]
        any(fails) && push!(effects, (state = s, t = t, cell = cell, covariates = check.dropped[fails]))
        return nothing
    end

    control_states = Dict{Tuple{Date, Date}, Vector{String}}()
    for row in findall(diff_df.treat .== 0)
        push!(get!(control_states, (diff_df.t[row], diff_df.r1[row]), String[]), diff_df.state[row])
    end
    for row in findall(diff_df.treat .== 1)
        s, t, r1 = diff_df.state[row], diff_df.t[row], diff_df.r1[row]
        target = fit.xbar[(s, t)]
        check!(s, t, target, (s, r1))
        for c in get(control_states, (t, r1), String[])
            check!(s, t, target, (c, t))
            check!(s, t, target, (c, r1))
        end
    end
    return effects
end

function warn_dropped_covariates_int(fit::CovariateFit, diff_df::DataFrame)

    # One warning that lists the ATT(s,t)'s found by dropped_covariate_effects, with the cell
    # and the covariates responsible. Nothing is said about covariates that are dropped
    # without affecting any ATT(s,t).
    effects = dropped_covariate_effects(fit, diff_df)
    isempty(effects) && return nothing
    describe(e) = "ATT($(e.state), $(e.t)): $(join(e.covariates, ", ")) in state $(e.cell[1]), period $(e.cell[2])"
    shown = join(describe.(first(effects, 5)), "; ")
    more = length(effects) > 5 ? "; and $(length(effects) - 5) more" : ""
    n_att = length(unique((e.state, e.t) for e in effects))
    @warn "Under ccc = \"int\", $n_att ATT(s,t) estimate" * (n_att == 1 ? " depends" : "s depend") *
          " on covariates that were dropped from the regression of a state-period cell (they do not vary there, or " *
          "are collinear with covariates listed before them). A dropped covariate's coefficient is set to zero in " *
          "that cell, but the cell is evaluated at the covariate means of the treated state, where other " *
          "coefficients that fit the cell's data equally well give a different estimate: " * shown * more *
          ". Which of several collinear covariates is dropped is decided by the order of 'covariates'."
    return nothing
end

function additive_att_shifts(shift::AbstractDict, diff_df::DataFrame, weighting::AbstractString)

    # The change in every ATT(s,t) when each cell's lambda changes by shift[cell] (ccc = "add",
    # where ATT(s,t) only depends on the lambdas): the change in the treated state's long
    # difference minus the weighted mean of the changes in its controls' long differences,
    # with the weights that build_att_st uses. Returns one value per treated long difference
    # that has a control.
    long_difference_shift(row) = shift[(diff_df.state[row], diff_df.t[row])] - shift[(diff_df.state[row], diff_df.r1[row])]
    control_rows = Dict{Tuple{Date, Date}, Vector{Int}}()
    for row in findall(diff_df.treat .== 0)
        push!(get!(control_rows, (diff_df.t[row], diff_df.r1[row]), Int[]), row)
    end
    changes = Float64[]
    for row in findall(diff_df.treat .== 1)
        controls = get(control_rows, (diff_df.t[row], diff_df.r1[row]), Int[])
        isempty(controls) && continue
        v = weighting in ["diff", "both"] ? Float64.(diff_df.n[controls]) : ones(length(controls))
        push!(changes, long_difference_shift(row) - sum(v .* long_difference_shift.(controls)) / sum(v))
    end
    return changes
end

function warn_additive_identification(fit::CovariateFit, diff_df::Union{DataFrame, Nothing} = nothing,
                                      weighting::AbstractString = "both")

    # ccc = "add": one warning if the first-stage regression leaves the lambdas undetermined
    # beyond the split of the slopes between states and periods (see additive_lambda_shifts!),
    # saying whether the ATTs are affected.
    # For every direction in which the lambdas can move, the ATT(s,t)'s are recomputed with the
    # lambdas moved along it. Under "add" the ATT(s,t)'s are linear in the lambdas, so the
    # change is the ATT(s,t) of the shifts themselves. An ATT(s,t) counts as changed when it
    # moves by more than IDENTIFICATION_RTOL times the largest shift of a lambda. If none
    # changes, every least-squares solution gives the same ATTs (for example a covariate that
    # is constant within the cells of a state and over time: it shifts every lambda of that
    # state by the same amount, which cancels in its long differences).
    # Without `diff_df` (didint_plot, which shows the lambdas themselves) the warning only says
    # that the lambdas are not uniquely determined.
    isempty(fit.lambda_shifts) && return nothing
    opening = "Under ccc = \"add\" the covariate slopes are less well determined than usual: the first-stage regression " *
              "has rank $(fit.additive_rank[1]) where $(fit.additive_rank[2]) is expected. The covariate-adjusted means " *
              "(the lambdas) are therefore not uniquely determined: they depend on which of the equally good " *
              "least-squares solutions is used. "
    involved(shifts) = join(unique(vcat([s.covariates for s in shifts]...)), ", ")
    causes = " Likely causes: a covariate that does not vary within any state-period cell of some state or of some " *
             "period, or covariates that are collinear within cells."
    if isnothing(diff_df)
        @warn opening * "Covariates involved: " * involved(fit.lambda_shifts) * "." * causes
        return nothing
    end

    affecting = eltype(fit.lambda_shifts)[]
    changed = falses(count(diff_df.treat .== 1))
    for s in fit.lambda_shifts
        changes = additive_att_shifts(s.shift, diff_df, weighting)
        moved = abs.(changes) .> IDENTIFICATION_RTOL * maximum(abs, values(s.shift))
        if any(moved)
            push!(affecting, s)
            changed[1:length(moved)] .|= moved
        end
    end
    n_changed = count(changed)
    if n_changed == 0
        @warn opening * "The ATT estimates are not affected: every solution gives the same ATTs. Covariates involved: " *
              involved(fit.lambda_shifts) * "."
    else
        @warn opening * "The ATT estimates are affected: $n_changed ATT(s,t) estimate" * (n_changed == 1 ? " depends" : "s depend") *
              " on the solution used. Covariates involved: " * involved(affecting) * "." * causes
    end
    return nothing
end

function build_att_st(diff_df::DataFrame, weighting::AbstractString, hc::AbstractString, vcov_lambda;
                      integration::Union{CovariateFit, Nothing} = nothing,
                      policies::Union{AbstractDict, Nothing} = nothing)

    # The ATT(s,t) table: one row per treated state s and post-treatment period t, with
    #   ATT(s,t) = long difference of s - weighted mean of the long differences of its controls
    # (Karim & Webb 2025, eq. 48). It is computed as the treat coefficient of the weighted
    # regression of the long differences on [1, treat], using s's long difference and those of
    # the controls in the same (t, r1) cell, which also gives its HC standard error.
    # With `integration`, every long difference in the regression gains its covariate term
    # evaluated at s's covariate means (eqs. 46-47), so a control's long difference depends
    # on the treated state it is compared with.
    # The single treated row has leverage 1, so HC2-HC4 standard errors are NaN. With exactly
    # two rows (one control) the regression SE is degenerate and the edge-case SE is used
    # instead. An (s,t) without a control or whose regression cannot be estimated is dropped.
    # `policies` gives each treated state's policy (one common policy if not given); the
    # controls do not depend on it.
    treated = sort(diff_df[diff_df.treat .== 1, :], [:state, :treated_time, :t])
    nrows = nrow(treated)
    float_column() = Vector{Union{Missing, Float64}}(missing, nrows)
    att_st = DataFrame(state = treated.state, gvar = treated.treated_time, t = treated.t, r1 = treated.r1,
                       periods_post_treat = treated.time_since_treatment,
                       policy = String[isnothing(policies) ? DEFAULT_POLICY : policies[s] for s in treated.state],
                       att = float_column(), se = float_column(), pval = float_column(),
                       n_t = treated.n_t, n = treated.n, n_controls = zeros(Int, nrows))

    # Rows of the control long differences of each (t, r1) cell
    control_rows = Dict{Tuple{Date, Date}, Vector{Int}}()
    for row in findall(diff_df.treat .== 0)
        push!(get!(control_rows, (diff_df.t[row], diff_df.r1[row]), Int[]), row)
    end

    for i in 1:nrows
        controls = diff_df[get(control_rows, (treated.t[i], treated.r1[i]), Int[]), :]
        att_st.n_controls[i] = nrow(controls)
        nrow(controls) == 0 && continue
        temp = vcat(treated[i:i, :], controls)
        X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
        Y = convert(Vector{Float64}, temp.diff)
        if !isnothing(integration)
            Y .+= [slope_change_term(integration, source, treated.t[i], treated.r1[i], treated.state[i]) for source in temp.state]
        end
        result_dict = final_regression_results(X, Y, W = control_weights(temp, weighting), hc = hc)
        att_st.att[i] = result_dict["beta_hat"]
        if length(Y) != 2
            att_st.se[i] = result_dict["beta_hat_se"]
        elseif isnothing(integration)
            att_st.se[i] = edge_case_se(temp.lambda_idx, vcov_lambda)
        else
            att_st.se[i] = edge_case_se_integrated(integration, treated.state[i], controls.state[1], treated.t[i], treated.r1[i])
        end
        att_st.pval[i] = result_dict["pval_att"]
    end
    att_st = att_st[(att_st.n_controls .> 0) .& .!ismissing.(att_st.att), :]
    if nrow(att_st) == 0
        error("Could not estimate a single ATT(s,t): no treated state has a long difference with a valid control in the same period.")
    end
    return att_st
end

function treated_unit_weights(att_st::AbstractDataFrame, weighting::AbstractString)

    # Weight u_{s,t} of each ATT(s,t) within its sub-aggregate: n_{s,t}, the treated state's
    # observations in period t, under weighting = "att" or "both", and 1 otherwise. With these
    # weights a sub-aggregate is the effect on the average treated observation of its group
    # (the weights pi_{s|g,t} of Austin, Strumpf, Karim & Webb 2026, Section 3.1).
    # The three averages and what `weighting` weights in each:
    #   weighting   control c in an ATT(s,t)   ATT(s,t) in a sub-aggregate   sub-aggregate j in the aggregate
    #   "both"      n_{c,t} + n_{c,r1}         n_{s,t}                       N_j = sum of its n_{s,t}
    #   "att"       1                          n_{s,t}                       N_j
    #   "diff"      n_{c,t} + n_{c,r1}         1                             1
    #   "none"      1                          1                             1
    return weighting in ["att", "both"] ? Float64.(att_st.n_t) : ones(nrow(att_st))
end

function sub_aggregate_results(att_st::DataFrame, agg::AbstractString, weighting::AbstractString,
                               hc::AbstractString; pooled::Bool = false)

    # Each sub-aggregate is the u-weighted mean of the ATT(s,t)'s in its group (for the policy
    # aggregations, the same groups within each policy).
    # Its HC SE comes from the weighted regression of those ATT(s,t)'s on a constant. This treats
    # them as independent, which can understate the SE when they share controls; the jackknife
    # and the randomization inference do not have that problem.
    # A group with a single ATT(s,t) reports that ATT(s,t)'s own SE and p-value. The weight of a
    # group in the aggregate ATT (the second weighting level) is the sum of its n_t under
    # weighting = "att"/"both".
    spec = SUB_AGG_SPECS[agg]
    group_cols = group_columns(agg)
    u = treated_unit_weights(att_st, weighting)

    # Groups in output order: by key, except states, which use custom_sort_order
    groups = unique(att_st[:, group_cols])
    if agg == "state"
        groups = groups[sortperm(custom_sort_order.(groups.state)), :]
    else
        sort!(groups, group_cols)
    end

    results = init_results(agg, nrow(groups); pooled = pooled)
    for i in 1:nrow(groups)
        in_group = trues(nrow(att_st))
        for col in group_cols
            in_group .&= att_st[!, col] .== groups[i, col]
        end
        for (col, att_col, _) in spec.keys
            results[i, col] = groups[i, att_col]
        end
        rows = att_st[in_group, :]
        if nrow(rows) == 1
            results[i, spec.suffix] = rows.att[1]
            results[i, "se_" * spec.suffix] = rows.se[1]
            results[i, "pval_" * spec.suffix] = rows.pval[1]
        else
            result_dict = final_regression_results(ones(nrow(rows), 1), convert(Vector{Float64}, rows.att),
                                                   W = u[in_group], hc = hc)
            results[i, spec.suffix] = result_dict["beta_hat"]
            results[i, "se_" * spec.suffix] = result_dict["beta_hat_se"]
            results[i, "pval_" * spec.suffix] = result_dict["pval_att"]
        end
        if weighting in ["att", "both"]
            results.weights[i] = sum(rows.n_t)
        end
    end
    return results
end

function custom_sort_order(s)

    # Sort key for state names: names that are not integers first, alphabetically, then
    # integer names in numeric order
    parsed = tryparse(Int, s)
    if isnothing(parsed)
        (0, lowercase(s))
    else
        (1, parsed)
    end
end

function scale_weights_final(results::DataFrame, weighting::AbstractString)

    # Weights of the sub-aggregates in the aggregate ATT: their sums of n_t scaled to sum to
    # one under weighting = "att"/"both", none (equal weights) otherwise
    if weighting in ["att", "both"]
        results.weights ./= sum(results.weights)
        results.weights = convert(Vector{Float64}, results.weights)
    elseif weighting in ["diff", "none"]
        results.weights = convert(Vector{Nothing}, results.weights)
    end
    return results
end

function aggregate_att(results::DataFrame, agg::AbstractString, weighting::AbstractString, hc::AbstractString)

    # Aggregate ATT: the weighted mean of the sub-aggregate ATTs, computed as the weighted
    # regression of the sub-aggregates on a constant, which also gives its HC SE. A single
    # sub-aggregate is the aggregate, so it reports that sub-aggregate's own SE and p-value.
    suffix = SUB_AGG_SPECS[agg].suffix
    results = scale_weights_final(results, weighting)
    if nrow(results) == 1
        results.agg_att[1] = results[1, suffix]
        results.se_agg_att[1] = results[1, "se_" * suffix]
        results.pval_agg_att[1] = results[1, "pval_" * suffix]
        return results
    end
    X = ones(nrow(results), 1)
    Y = convert(Vector{Float64}, results[!, suffix])
    W = results.weights
    result_dict = final_regression_results(X, Y, W = W, hc = hc)
    results.agg_att[1] = result_dict["beta_hat"]
    results.se_agg_att[1] = result_dict["beta_hat_se"]
    results.pval_agg_att[1] = result_dict["pval_att"]
    return results
end

function policy_pool_weights(table::AbstractDataFrame, policies::AbstractVector, policy_weighting::AbstractString)

    # Weight omega_p of each policy in the pooled ATT, from a table of ATT(s,t)'s: "equal" 1,
    # "treated_obs" the policy's treated observations (sum of n_t), "treated_cells" its number
    # of ATT(s,t)'s, "treated_states" its number of treated states. Normalized to sum to 1.
    weights = map(policies) do p
        rows = table.policy .== p
        if policy_weighting == "equal"
            1.0
        elseif policy_weighting == "treated_obs"
            Float64(sum(table.n_t[rows]))
        elseif policy_weighting == "treated_cells"
            Float64(count(rows))
        else
            Float64(length(unique(table.state[rows])))
        end
    end
    return weights ./ sum(weights)
end

function policy_att_results(results::DataFrame, att_st::DataFrame, agg::AbstractString, weighting::AbstractString,
                            hc::AbstractString, policy)

    # Top-level estimates of the policy aggregations (Austin, Strumpf, Karim & Webb 2026). Each
    # policy's ATT_p is the weighted mean of that policy's sub-aggregates, with its HC SE from
    # the weighted regression of those sub-aggregates on a constant, as aggregate_att does across
    # all sub-aggregates (a single sub-aggregate reports its own SE). It is written on the first
    # row of the policy. For agg = "policy" the sub-aggregate is ATT_p itself.
    suffix = SUB_AGG_SPECS[agg].suffix
    policies = unique(results.policy)

    # Weights of the sub-aggregates, normalized within each policy. agg = "policy" has no
    # sub-aggregates below the policy's ATT, so it has none.
    if weighting in ["att", "both"] && agg != "policy"
        for p in policies
            rows = results.policy .== p
            results.weights[rows] = results.weights[rows] ./ sum(results.weights[rows])
        end
        results.weights = convert(Vector{Float64}, results.weights)
    else
        results.weights = Vector{Nothing}(nothing, nrow(results))
    end

    if agg != "policy"
        for p in policies
            rows = findall(results.policy .== p)
            if length(rows) == 1
                results[rows[1], "att_p"] = results[rows[1], suffix]
                results[rows[1], "se_att_p"] = results[rows[1], "se_" * suffix]
                results[rows[1], "pval_att_p"] = results[rows[1], "pval_" * suffix]
            else
                result_dict = final_regression_results(ones(length(rows), 1), convert(Vector{Float64}, results[rows, suffix]),
                                                       W = results.weights[rows], hc = hc)
                results[rows[1], "att_p"] = result_dict["beta_hat"]
                results[rows[1], "se_att_p"] = result_dict["beta_hat_se"]
                results[rows[1], "pval_att_p"] = result_dict["pval_att"]
            end
        end
    end

    # Pooled ATT across policies: ATT_pool = sum_p omega_p ATT_p, with its HC SE from the
    # weighted regression of the ATT_p's on a constant (number of policies - 1 degrees of
    # freedom). With a single policy it is that policy's ATT.
    if policy.pooled
        first_rows = [findfirst(results.policy .== p) for p in policies]
        omega = policy_pool_weights(att_st, policies, policy.weighting)
        results.policy_weights[first_rows] = omega
        if length(policies) == 1
            results.att_pool[1] = results.att_p[1]
            results.se_att_pool[1] = results.se_att_p[1]
            results.pval_att_pool[1] = results.pval_att_p[1]
        else
            result_dict = final_regression_results(ones(length(policies), 1), convert(Vector{Float64}, results.att_p[first_rows]),
                                                   W = omega, hc = hc)
            results.att_pool[1] = result_dict["beta_hat"]
            results.se_att_pool[1] = result_dict["beta_hat_se"]
            results.pval_att_pool[1] = result_dict["pval_att"]
        end
    end
    return results
end

###########################################################################################
# Section 4: Regressions and edge-case standard errors                                    #
###########################################################################################

function final_regression_results(X::Matrix, Y::Vector; W::Vector = [nothing], hc::AbstractString = "hc1")

    # Least squares of Y on X, weighted by W if weights are given, for the small second-stage
    # regressions (an ATT(s,t) from long differences, a sub-aggregate from ATT(s,t)'s, the
    # aggregate from sub-aggregates). Rows with a missing value are dropped. Returns the last
    # coefficient ("beta_hat"), its HC standard error ("beta_hat_se") and its two-sided p-value
    # from a t distribution with n - k degrees of freedom ("pval_att"). The standard error and
    # p-value are missing when there are no residual degrees of freedom.
    failed = Dict("beta_hat" => missing, "beta_hat_se" => missing, "pval_att" => missing)
    weighted = !(eltype(W) <: Nothing)

    # Drop rows with missing values
    valid = [!any(ismissing, X[i, :]) for i in 1:length(Y)] .& .!ismissing.(Y)
    if weighted
        valid .&= .!ismissing.(W)
        W = W[valid]
    end
    X = X[valid, :]
    Y = Y[valid]

    n, k = length(Y), size(X, 2)
    if n < k
        @warn "Insufficient observations after dropping missing values (n = $n, k = $k)"
        return failed
    end

    # Weighted least squares is ordinary least squares on rows scaled by the square root of
    # their (normalized) weights
    if weighted
        sw = sqrt.(W ./ sum(W))
        X = X .* sw
        Y = Y .* sw
    end
    beta_hat = safe_solve(X, Y)
    if ismissing(beta_hat)
        return failed
    end

    if n > k
        resid = Y - X * beta_hat
        beta_hat_se = sqrt(compute_hc_covariance(X, resid, hc)[k, k])
        pval_att = 2 * (1 - cdf(TDist(n - k), abs(beta_hat[k] / beta_hat_se)))
    else
        beta_hat_se = missing
        pval_att = missing
    end
    return Dict("beta_hat" => beta_hat[k], "beta_hat_se" => beta_hat_se, "pval_att" => pval_att)
end

# Edge case: an ATT(s,t) with exactly one control c. Its regression has two rows and two
# parameters, fits both rows exactly, and has no residuals to build a standard error from.
# With edgecase = true the variance is instead recovered from the first stage:
#   ATT(s,t) = (lambda_{s,t} - lambda_{s,r1}) - (lambda_{c,t} - lambda_{c,r1}),
# a combination of four first-stage estimates whose covariances the first stage can provide.
# Both functions below assume that observations in different state-period cells are
# independent, as in repeated cross-sections. With panel data the same individuals appear in a
# state's r1 and t cells, their estimates are correlated, and that covariance is omitted.

function edge_case_se(lambda_indices, vcov_lambda)

    # Edge-case SE without integrated covariate terms (every ccc without covariates; "hom",
    # "state" and "time" with covariates). With a = (-1, 1, 1, -1) on (lambda_{s,r1},
    # lambda_{s,t}, lambda_{c,r1}, lambda_{c,t}),
    #   Var(ATT) = a' V a
    #            = the four variances - 2 Cov(lambda_{s,r1}, lambda_{s,t}) - 2 Cov(lambda_{c,r1}, lambda_{c,t})
    #              - 2 Cov(diff_s, diff_c),
    # where V is the covariance matrix of the lambdas (vcov_lambda) and diff_s, diff_c are the
    # two long differences. Which covariances are nonzero depends on the ccc (under "int" all
    # of them are zero, as every cell is estimated separately), so one formula serves them all.
    # Each element of lambda_indices is (index of lambda_r1, index of lambda_t) for one state;
    # Cov(diff_s, diff_c) is only computed correctly if both states use this same order.
    # Without the covariance matrix (edgecase = false) the SE is missing.
    isnothing(vcov_lambda) && return missing

    #s1_l1 short for state1_lambda1, etc
    s1_l1_var = vcov_lambda[lambda_indices[1][1], lambda_indices[1][1]]
    s1_l2_var = vcov_lambda[lambda_indices[1][2], lambda_indices[1][2]]
    s2_l1_var = vcov_lambda[lambda_indices[2][1], lambda_indices[2][1]]
    s2_l2_var = vcov_lambda[lambda_indices[2][2], lambda_indices[2][2]]

    # If any of the variance terms for lambda are 0, NaN, or Inf then the var(ATT) is not computable
    if any(x -> x == 0 || isnan(x) || isinf(x), [s1_l1_var, s1_l2_var, s2_l1_var, s2_l2_var])
        return missing
    end

    se_att = sqrt(s1_l1_var + s1_l2_var + s2_l1_var + s2_l2_var - 2*vcov_lambda[lambda_indices[1][1], lambda_indices[1][2]] -
             2*vcov_lambda[lambda_indices[2][1], lambda_indices[2][2]] + -2*(vcov_lambda[lambda_indices[1][2], lambda_indices[2][2]] -
             vcov_lambda[lambda_indices[1][2], lambda_indices[2][1]] - vcov_lambda[lambda_indices[1][1], lambda_indices[2][2]] +
             vcov_lambda[lambda_indices[1][1], lambda_indices[2][1]]))
    return se_att

end

function edge_case_se_integrated(fit::CovariateFit, treated, control, t, r1)

    # Edge-case SE of an integrated ATT(s,t) with a single control under ccc = "int" with
    # covariates, where the covariate terms of the long differences (Karim & Webb 2025,
    # eqs. 46-47) do not cancel. With xbar = xbar_{s,t}, the covariate means of the treated
    # cell, ATT(s,t) combines the fitted values of four cells' regressions at xbar:
    #   ATT(s,t) = (lambda_{s,t} + gamma_{s,t}' xbar) - (lambda_{s,r1} + gamma_{s,r1}' xbar)
    #            - (lambda_{c,t} + gamma_{c,t}' xbar) + (lambda_{c,r1} + gamma_{c,r1}' xbar).
    # The variance is taken conditional on the observed covariates, so xbar is fixed. The four
    # regressions use separate observations and are independent, so the variances of the four
    # fitted values add up:
    #   Var(ATT(s,t) | X) = sum over the 4 cells k of
    #                       Var(lambda_k) + 2 xbar' Cov(gamma_k, lambda_k) + xbar' Var(gamma_k) xbar
    #                     = sum over k of a' V_k a,   a = (1, xbar),
    # with V_k the HC covariance of cell k's (lambda, gamma). The sampling variability of
    # xbar_{s,t} itself is not included. Without covariates this is the sum of the four
    # Var(lambda_k), which is what edge_case_se computes.
    # Like every edge-case SE it assumes that the cells are independent, which holds for
    # repeated cross-sections but not for panel data.
    cells = [(treated, t), (treated, r1), (control, t), (control, r1)]
    all(haskey(fit.cell_vcov, cell) for cell in cells) || return missing

    a = vcat(1.0, fit.xbar[(treated, t)])
    variance = sum(a' * fit.cell_vcov[cell] * a for cell in cells)
    return isfinite(variance) && variance > 0 ? sqrt(variance) : missing
end

###########################################################################################
# Section 5: Jackknife                                                                    #
###########################################################################################

function jackknife_procedure(results::DataFrame, att_st::DataFrame, diff_df::DataFrame, agg::AbstractString,
                             weighting::AbstractString, use_pre_controls::Bool,
                             truejack::Bool, data_copy, ccc, covariates_to_include, treated_states, treatment_times,
                             time_to_index, match_to_these_dates;
                             integration::Union{CovariateFit, Nothing} = nothing,
                             policy = default_policy_options(treated_states))

    # Delete-one-state (cluster) jackknife: the CV3 (CRV3) variance estimator, centered at the
    # full-sample estimate (MacKinnon, Nielsen & Webb 2023). Each state is deleted as a treated
    # state and as a control, and the same weighted means are recomputed, at every level of the
    # output: the sub-aggregates, the aggregate ATT or each policy's ATT, and the pooled ATT.
    # Deleting a state leaves the other states' lambdas and slopes unchanged under ccc = "int"
    # and "state", or is taken to under truejack = false, so the batched engine evaluates all
    # deletions in one pass. With truejack = true (ccc not "int" or "state") the first stage is
    # re-estimated without the state, since its lambdas then depend on the other states.

    # A jackknife needs at least 2 treated and 2 control states; otherwise every jackknife
    # SE and p-value stays missing
    n_treated = length(unique(diff_df[diff_df.treat .== 1, :state]))
    n_control = length(unique(diff_df[diff_df.treat .== 0, :state]))
    if n_control < 2 || n_treated < 2
        @warn "Need at least 2 treated states and 2 control states to do jackknife procedure."
        return results
    end

    # PART ONE: THE DELETE-ONE ESTIMATES, ONE ROW PER DELETED STATE
    # top_jacks holds the aggregate ATT, or for the policy aggregations the pooled ATT;
    # policy_jacks holds each policy's ATT_p; sub_jacks the sub-aggregates
    cells = att_st_cells(att_st)
    observed_labels = Dict(zip(treated_states, treatment_times))
    all_states = unique(diff_df.state)
    n_states = length(all_states)
    policy_names, policy_code = policy_coding(results, agg, policy.labels)
    if !truejack
        # Fast jackknife: every row of the assignment matrix is the observed assignment, and
        # row h switches state h off as a treated state and as a control. The weights
        # renormalize automatically.
        engine = build_engine(diff_df, cells, all_states, weighting; integration = integration)
        L_observed, P_observed = assignment_row(all_states, observed_labels, policy.labels, policy_code)
        est = engine_estimates(engine, results, agg, weighting, use_pre_controls,
                               repeat(L_observed', n_states), repeat(P_observed', n_states);
                               active = .!Matrix{Bool}(I, n_states, n_states),
                               policy_names = policy_names, policy_weighting = policy.weighting)
        sub, att_p, top = est.sub, est.att_p, est.top
    else
        # True jackknife: re-estimate the first stage without state h, rebuild the long
        # differences, and evaluate the observed assignment of the remaining states. The
        # re-estimation decides anew which collinear covariates to drop, as removing a state
        # can create or remove collinearity; it does not repeat the full sample's warnings.
        sub = fill(NaN, n_states, nrow(results))
        att_p = fill(NaN, n_states, length(policy_names))
        top = fill(NaN, n_states)
        for (h, without) in enumerate(all_states)
            lambda_h, _, fit_h = first_stage(data_copy[data_copy.state_71X9yTx .!= without, :], ccc, covariates_to_include, "hc0", false; warn = false)
            diff_h = construct_diff_df(lambda_h, treated_states, treatment_times, time_to_index, match_to_these_dates,
                                       use_pre_controls, data_copy; jackknife = true, without = without)
            remaining = filter(!=(without), all_states)
            engine_h = build_engine(diff_h, cells, remaining, weighting;
                                    integration = isnothing(integration) ? nothing : fit_h)
            L_h, P_h = assignment_row(remaining, observed_labels, policy.labels, policy_code)
            est = engine_estimates(engine_h, results, agg, weighting, use_pre_controls, reshape(L_h, 1, :), reshape(P_h, 1, :);
                                   policy_names = policy_names, policy_weighting = policy.weighting)
            sub[h, :] = est.sub
            att_p[h, :] = est.att_p
            top[h] = est.top[1]
        end
    end
    # An estimate that a deletion leaves undefined is missing
    defined_or_missing(x) = isfinite(x) ? x : missing
    sub_jacks = defined_or_missing.(sub)
    policy_jacks = defined_or_missing.(att_p)
    top_jacks = defined_or_missing.(top)
    result_keys = [result_key(results, i, agg) for i in 1:nrow(results)]

    # PART TWO: JACKKNIFE STANDARD ERRORS AND P-VALUES
    # Treated states behind a set of ATT(s,t)'s, and control states in their cells
    control_rows = diff_df[diff_df.treat .== 0, :]
    function state_counts(in_group)
        group_cells = Set(zip(att_st.t[in_group], att_st.r1[in_group]))
        n_group_control = length(unique(control_rows.state[[cell in group_cells for cell in zip(control_rows.t, control_rows.r1)]]))
        return length(unique(att_st.state[in_group])), n_group_control
    end

    # Sub-aggregates. A sub-aggregate gets no jackknife SE unless its ATT(s,t)'s involve at
    # least 2 treated states and their cells at least 2 control states, and every deletion
    # leaves it defined. The single-state sub-aggregates of "state" and "sgt" never do.
    suffix = SUB_AGG_SPECS[agg].suffix
    group_cols = group_columns(agg)
    att_st_keys = [Tuple(att_st[r, col] for col in group_cols) for r in 1:nrow(att_st)]
    for i in 1:nrow(results)
        n_group_treated, n_group_control = state_counts(att_st_keys .== Ref(result_keys[i]))
        if n_group_treated < 2 || n_group_control < 2 || any(ismissing, sub_jacks[:, i])
            results[i, "jknifese_" * suffix] = missing
            results[i, "jknifepval_" * suffix] = missing
        else
            se, pval = jackknife_se_pval(collect(skipmissing(sub_jacks[:, i])), results[i, suffix])
            results[i, "jknifese_" * suffix] = se
            results[i, "jknifepval_" * suffix] = pval
        end
    end

    # Aggregate ATT, over the deletions that leave it defined. A deletion that removes every
    # treated state of a sub-aggregate leaves that sub-aggregate undefined (so it has no
    # jackknife SE, above), and the aggregate for that deletion is the mean over the remaining
    # sub-aggregates.
    if !is_policy_agg(agg)
        results.jknifese_agg_att[1], results.jknifepval_agg_att[1] = jackknife_se_pval(collect(skipmissing(top_jacks)), results.agg_att[1])
        return results
    end

    # Each policy's ATT_p, if the policy has at least 2 treated states and its cells at least 2
    # control states (for agg = "policy" it is the sub-aggregate handled above)
    if agg != "policy"
        for (j, p) in enumerate(policy_names)
            n_policy_treated, n_policy_control = state_counts(att_st.policy .== p)
            if n_policy_treated >= 2 && n_policy_control >= 2
                row = findfirst(results.policy .== p)
                results[row, "jknifese_att_p"], results[row, "jknifepval_att_p"] =
                    jackknife_se_pval(collect(skipmissing(policy_jacks[:, j])), results[row, "att_p"])
            end
        end
    end

    # Pooled ATT across policies. A deletion that removes a policy's only treated state leaves
    # the pool defined over the other policies.
    if policy.pooled
        results.jknifese_att_pool[1], results.jknifepval_att_pool[1] = jackknife_se_pval(collect(skipmissing(top_jacks)), results.att_pool[1])
    end
    return results
end

function jackknife_se_pval(estimates::AbstractVector, full_sample::Real)

    # Jackknife SE sqrt((n-1)/n * sum_h (theta_h - theta)^2) over the n delete-one estimates
    # theta_h, centered at the full-sample estimate theta (the CV3, or CRV3, estimator of MacKinnon,
    # Nielsen & Webb 2023), and its two-sided p-value from a t(n-1) distribution
    n = length(estimates)
    se = n >= 2 ? sqrt(sum((estimates .- full_sample).^2) * ((n - 1) / n)) : missing
    pval = n - 1 > 0 && !ismissing(se) ? 2 * (1 - cdf(TDist(n - 1), abs(full_sample / se))) : missing
    return se, pval
end

###########################################################################################
# Section 6: Randomization inference                                                      #
###########################################################################################

# Number of permutations that randomization inference evaluates at a time
const RI_BATCH_SIZE = 500

function randomization_inference(results::DataFrame, att_st::DataFrame, long_diffs::DataFrame, agg::AbstractString,
                                 weighting::AbstractString, use_pre_controls::Bool,
                                 treated_states, treatment_times, nperm::Int, ri_kernel_cutoff::Int, seed::Number;
                                 integration::Union{CovariateFit, Nothing} = nothing,
                                 policy = default_policy_options(treated_states))

    # Randomization inference (MacKinnon & Webb 2020): give the observed (treatment time, policy)
    # labels to other states, recompute the same estimates for each alternative assignment, and
    # compare them with the observed estimates. A p-value is reported at every level of the
    # output: every sub-aggregate, the aggregate ATT or each policy's ATT, and the pooled ATT.
    # long_diffs holds every state's long difference in every cell. The first stage does not
    # depend on the assignment, so the batched engine evaluates all permutations from the same
    # lambdas, slopes and covariate means.

    # PART ONE: THE ALTERNATIVE ASSIGNMENTS OF THE (TREATMENT TIME, POLICY) LABELS
    # All of them if there are at most `nperm`, otherwise `nperm` distinct random ones
    engine = build_engine(long_diffs, att_st_cells(att_st), unique(long_diffs.state), weighting;
                          integration = integration)
    original_treated = unique(long_diffs[long_diffs.treat .== 1, [:state, :treated_time]])
    policy_names, policy_code = policy_coding(results, agg, policy.labels)
    L_observed, P_observed = assignment_row(engine.states, Dict(zip(original_treated.state, original_treated.treated_time)),
                                            policy.labels, policy_code)
    rng = Xoshiro(seed)
    L, P, _ = ri_assignment_draws(L_observed, P_observed, nperm, rng)
    if size(L, 1) < nperm
        @warn "'nperm' was set to $nperm but only $(size(L, 1)) unique permutations exist. \n
                Setting 'nperm' to $(size(L, 1))."
        nperm = size(L, 1)
    end
    if nperm < 399
        @warn "'nperm' is less than 399!"
    end

    # PART TWO: COMPUTE THE PLACEBO SUB-AGGREGATE AND TOP-LEVEL ATTS FOR EVERY PERMUTATION
    # att_ri holds the placebo aggregate ATT, or for the policy aggregations the pooled ATT;
    # att_ri_policy holds each policy's placebo ATT_p. The statistic is the estimate itself (not
    # a t-statistic), and a placebo sub-aggregate is compared with the observed sub-aggregate
    # that has the same key. "state" and "sgt" have one treated state per sub-aggregate and use
    # stand-in states instead (see state_level_ri_estimates).
    att_ri = Vector{Float64}(undef, nperm)
    att_ri_policy = Matrix{Float64}(undef, nperm, length(policy_names))
    att_ri_sub = Matrix{Float64}(undef, nperm, nrow(results))
    observed_cohort = Dict(zip(treated_states, treatment_times))
    for batch in Iterators.partition(1:nperm, RI_BATCH_SIZE)
        if agg in ["state", "sgt"]
            est = state_level_ri_estimates(engine, results, agg, weighting, use_pre_controls, L[batch, :], observed_cohort, rng)
        else
            est = engine_estimates(engine, results, agg, weighting, use_pre_controls, L[batch, :], P[batch, :];
                                   policy_names = policy_names, policy_weighting = policy.weighting)
            att_ri_policy[batch, :] = est.att_p
        end
        att_ri_sub[batch, :] = est.sub
        att_ri[batch] = est.top
    end

    # PART THREE: COMPUTE P-VALS BASED ON THE PLACEBO ESTIMATES
    # The kernel-smoothed p-value when at most ri_kernel_cutoff permutations were used, whether
    # they are all the permutations that exist or a sample of them; otherwise the ordinary
    # p-value. With a single permutation there is nothing to smooth.
    use_kernel = 2 <= nperm <= ri_kernel_cutoff
    pvalue(placebos, observed) = use_kernel ? ri_pvalue_kernel(placebos, observed) : ri_pvalue(placebos, observed)
    suffix = SUB_AGG_SPECS[agg].suffix
    for i in 1:nrow(results)
        results[i, "ri_pval_" * suffix] = pvalue(att_ri_sub[:, i], results[i, suffix])
    end
    if !is_policy_agg(agg)
        results.ri_pval_agg_att[1] = pvalue(att_ri, results.agg_att[1])
    else
        if agg != "policy"
            for (j, p) in enumerate(policy_names)
                row = findfirst(results.policy .== p)
                results[row, "ri_pval_att_p"] = pvalue(att_ri_policy[:, j], results[row, "att_p"])
            end
        end
        if policy.pooled
            results.ri_pval_att_pool[1] = pvalue(att_ri, results.att_pool[1])
        end
    end
    # nperm is the number of permutations used
    results.nperm[1] = nperm
    return results
end

# Relative tolerance for ties between a placebo estimate and the observed estimate
const RI_TIE_RTOL = sqrt(eps(Float64))

function ri_pvalue(att_ri::AbstractVector, observed::Real)

    # RI p-value p = (1 + #{|T*| >= |T|}) / (S + 1) over the S finite placebo estimates T*,
    # the more conservative of the two versions in MacKinnon & Webb (2020).
    # Ties count as hits. A placebo assignment that leaves a group's treated states and
    # its controls unchanged (common with use_pre_controls = true) reproduces the observed
    # estimate exactly in exact arithmetic, but the two are computed along different paths
    # and can differ by rounding (~1e-16), which also varies with compiler flags such as
    # --check-bounds. Treating |T*| within a relative RI_TIE_RTOL of |T| as a tie keeps
    # the count deterministic. Returns missing if there is no finite placebo estimate or if
    # the observed estimate is not finite.
    finite_att_ri = filter(isfinite, att_ri)
    n_placebo = length(finite_att_ri)
    (n_placebo == 0 || !isfinite(observed)) && return missing
    threshold = abs(observed) * (1 - RI_TIE_RTOL)
    n_hits = count(x -> abs(x) >= threshold, finite_att_ri)
    return (1 + n_hits) / (n_placebo + 1)
end

function ri_pvalue_kernel(att_ri::AbstractVector, observed::Real)

    # Kernel-smoothed RI p-value (Racine & MacKinnon 2007):
    #   p = 1 - (1/S) sum_j Phi((|T| - |T*_j|) / h),   h = 1.575 sd(T*) S^(-4/9),
    # over the S finite placebo estimates T*, with Phi the standard normal CDF and sd the sample
    # standard deviation of the signed T*. With few permutations the ordinary p-value can only
    # take a few values; smoothing the indicator |T*| >= |T| lets it move continuously with T.
    # Falls back to the ordinary p-value if S < 2 or if the T* do not vary (their sd is zero up
    # to rounding, or not finite). Returns missing if the observed estimate is not finite.
    finite_att_ri = filter(isfinite, att_ri)
    n_placebo = length(finite_att_ri)
    (n_placebo == 0 || !isfinite(observed)) && return missing
    sd = n_placebo >= 2 ? std(finite_att_ri) : NaN
    if !(isfinite(sd) && sd > RI_TIE_RTOL * maximum(abs, finite_att_ri))
        return ri_pvalue(att_ri, observed)
    end
    h = 1.575 * sd * n_placebo^(-4 / 9)
    return 1 - mean(cdf.(Normal(), (abs(observed) .- abs.(finite_att_ri)) ./ h))
end
