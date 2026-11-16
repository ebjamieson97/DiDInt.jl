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
           treated_weight::AbstractString = "treated",
           policy_aggregation::Bool = false,
           policy_weighting::AbstractString = "equal",
           ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
           freq::Union{AbstractString, Nothing} = nothing,
           freq_multiplier::Number = 1,
           start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
           end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
           nperm::Number = 999,
           ri_kernel_cutoff::Number = 99,
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
The arguments `treated_states` and `treated_times` should be entered such that the first element in
`treated_states` refers to the state treated at the date entered as the first element in `treated_times`,
the second element in `treated_states` refers to the state treated at the date entered as the second element
in `treated_times`, and so on. 

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
    Name of the column which indicates time of first treatment for each state.
- `treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing`
    A vector of strings (or a single string) noting the treated state(s).
- `treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing`
    A vector (or single entry) denoting the associated treatment times of the 'treated_states'.
- `pvar::Union{AbstractString, Symbol, Nothing} = nothing`
    Name of the column which indicates the policy each state adopts (constant within a state;
    missing, `0`, `"0"` or `""` for never-treated states). Use either `pvar` or `treated_policies`.
- `treated_policies::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing`
    A vector (or single entry) denoting the policies of the 'treated_states', in the same order.
    If neither `pvar` nor `treated_policies` is given, every treated state has the same policy.
    Policies only matter for the policy aggregations: controls are never-treated states (and, with
    `notyet`, not-yet-treated states) whatever their policy.

## Model Specification
- `ccc::AbstractString = "int"`
    Specify which version of DID-INT should be used.
    Options are: `"hom"`, `"time"`, `"state"`, "`add`", and `"int"`.
    The long differences are the integrated ones of Karim & Webb (eqs. 46-47): the change in a state's
    covariate slopes is evaluated at the treated state's covariate means. That covariate term cancels
    in every ATT(s,t) except under `"int"`, so only `"int"` estimates with covariates depend on it.
- `agg::AbstractString = "cohort"` 
    Enter the aggregation method as a string.
    Options are: `"cohort"`, `"simple"`, `"state"`, `"sgt"`, `"time"`, `"calendar"`, and the policy
    aggregations `"policy"`, `"policy_gt"`, `"policy_cohort"`, `"policy_calendar"`, `"policy_time"`.
    Every sub-aggregate is a weighted mean of ATT(s,t) estimates (one per treated state and
    post-treatment period), and `agg_att` is a weighted mean of the sub-aggregates.
    The policy aggregations form the same sub-aggregates within each policy (`"policy"` has one per
    policy) and report each policy's ATT in `att_p`, on the first row of that policy, instead of `agg_att`.
- `weighting::AbstractString = "both"`
    Specify which weighting method should be used.
    Options are: `"both"`, `"att"`, `"diff"`, or `"none"`.
- `treated_weight::AbstractString = "treated"`
    How an ATT(s,t) is weighted within its sub-aggregate when `weighting` uses observation counts.
    `"treated"` weights it by the treated state's post-period observations n_{s,t} (under `"att"` and
    `"both"`); `"population"` weights it by the observations behind its long difference,
    n_{s,t} + n_{s,r1} (under `"diff"` and `"both"`). Options are `"treated"` and `"population"`.
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
    or, `nothing` (default).
- `notyet::Bool = false`
    Determine if pre-treatment periods from treated states should be used as controls.
- `ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing`
    A dictionary specifying which category in a categorical variable should be used
    as the reference (baseline) category.

## Date Processing & Period Grid Construction
- `date_format::Union{AbstractString, Nothing} = nothing`
    Date format (e.g. "yyyy" or "yyyy-mm-dd") to be used when parsing string dates from
    the time column, or `start_date`, `end_date`, and `treatment_times` arguments.
- `freq::Union{AbstractString, Nothing} = nothing`
    A string indicating the desired timeframe of a period for the analysis for staggered adoption scenarios.
    Options are: `"year"`, `"month"`, `"week"`, `"day"`.
- `freq_multiplier::Number = 1`
    An integer by which the 'freq' argument should be multiplied in a staggered adoption scenario, e.g. if a two-year
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
    to be considered when performing the randomization inference.
- `ri_kernel_cutoff::Number = 99`
    When the randomization inference uses `ri_kernel_cutoff` or fewer alternative assignments (because
    `nperm` is that small, or because no more unique assignments exist), its p-values are kernel-smoothed
    (Racine & MacKinnon 2007) instead of counting the placebo estimates at least as large as the observed one.
    Set it to `0` to never smooth. The column `ri_procedure` reports which p-value was used.
- `seed::Number = rand(1:1000000)`
    An integer to set the random seed for the randomization inference procedure.
- `hc::Union{AbstractString, Number} = "hc1"`
    Specify which heteroskedasticity-consistent covariance matrix estimator (HCCME) should be used.
    Options are `0`, `1`, `2`, `3`, and `4` (or `"hc0"`, `"hc1"`, `"hc2"`, `"hc3"`, `"hc4"`).
- `truejack::Bool = false`
    When aggregation is set to either `"add"`, `"time"`, or `"hom"`, then in order to get valid jackknife
    standard errors, we need to re-estimate the DID-INT model from square one (running the large FixedEffectsModels).
    This is because the covariate effects in those cases depend on values from across states, so dropping a state will change
    the lambda values, this is not true for the aggregation options of `"int"` or `"state"`
- `edgecase::Bool = false`
    When there are at least three states for which a long difference can be calculated, then standard errors for an ATT can be calculated
    directly from the regression of differences on the intercept term and treatment indicator. However, in cases when there are only two states
    for a particular long difference regression, the standard error can still be constructed from the variance and covariance terms of the relevant
    means. Calculating these `edgecase` standard errors can be computationally expensive, therefore by default this argument is set to `false`, but
    may be toggled on by setting `edgecase = true`.
    Under `ccc = "int"` with covariates the ATT also depends on the treated state's covariate means, and the `edgecase` standard error
    includes their sampling variance (the unconditional variance). Every `edgecase` standard error assumes that observations in different
    state-period cells are independent, as in repeated cross-sections. With panel data, where the same individuals are observed in several
    periods, it omits the covariances between a state's periods and should be read with caution. Not available for `ccc = "add"`.

# Returns
A DataFrame of results including the estimate of the ATT as well as standard errors and p-values.

# Citations
- Karim & Webb (2025). "Good Controls Gone Bad: Difference-in-Differences with Covariates". https://arxiv.org/abs/2412.14447
- MacKinnon & Webb (2020). "Randomization inference for difference-in-differences with few treated clusters". https://doi.org/10.1016/j.jeconom.2020.04.024 

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
                treated_weight::AbstractString = "treated",
                policy_aggregation::Bool = false,
                policy_weighting::AbstractString = "equal",
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                nperm::Number = 999,
                ri_kernel_cutoff::Number = 99,
                seed::Number = rand(1:1000000),
                use_pre_controls::Bool = false,
                notyet::Union{Nothing, Bool} = nothing,
                hc::Union{AbstractString, Number} = "hc1",
                wrapper::Union{AbstractString, Nothing} = nothing,
                truejack::Bool = false,
                edgecase::Bool = false
                )

    # Check hc args
    hc = hc_checks(hc)

    # Check wrapper input
    wrapper = init_wrapper_check(wrapper)

    # Check typeof(data)
    data_type_check(data)

    # Let notyet override use_pre_controls
    use_pre_controls = isnothing(notyet) ? use_pre_controls : notyet

    # Turn outcome, state, time, and gvar to strings if they were inputted as symbols
    outcome, state, time, gvar = symbol_to_string(outcome, state, time, gvar)

    # Ensure column for covariates are passed as strings
    covariates = isnothing(covariates) ? covariates : string.(covariates)

    # Determine if the gvar method or the treatment_times & treated_states method is being used
    treated_states, treatment_times = gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                                              event = false, check = "estimate_didint")

    # Policy of each treated state, from pvar or treated_policies (update spec §6.1)
    state_policies = policy_assignments(pvar, treated_policies, treated_states, gvar, data, state)

    # Do checks for start_date and end_date
    start_date = start_end_date_checks(start_date, "start_date", date_format)
    end_date = start_end_date_checks(end_date, "end_date", date_format)           

    # Check that seed is set correctly
    seed = Int(round(seed))
    if seed <= 0 
        error("'seed' must be a positive integer.")
    end 

    # Check that agg args are passed correctly
    agg = lowercase(agg)
    if agg == "none"
        error("agg = \"none\" has been removed; every agg reports the aggregate ATT in agg_att. Use e.g. agg = \"simple\".")
    end
    agg_options = ["cohort", "state", "simple", "sgt", "time", "calendar",
                   "policy", "policy_gt", "policy_cohort", "policy_calendar", "policy_time"]
    if !(agg in agg_options)
        error("'agg' must be one of $(agg_options)")
    end

    # Check the policy pooling args (update spec §6.4)
    policy_weighting = lowercase(policy_weighting)
    policy_weighting_options = ["equal", "treated_obs", "treated_cells", "treated_states"]
    if !(policy_weighting in policy_weighting_options)
        error("'policy_weighting' must be one of $(policy_weighting_options)")
    end
    if !policy_aggregation && policy_weighting != "equal"
        error("'policy_weighting' only applies when 'policy_aggregation = true'.")
    end
    # Pooling across policies only exists for the policy aggregations; with any other agg,
    # policy_aggregation (and policy_weighting) have nothing to pool and are ignored
    if !is_policy_agg(agg)
        policy_aggregation = false
    end

    # Check that weighting arg is passed correctly
    weighting = lowercase(weighting)
    weighting_options = ["both", "none", "att", "diff"]
    if !(weighting in weighting_options)
        error("'weighting' must be one of $(weighting_options)")
    end

    # Check that treated_weight arg is passed correctly (update spec §4.1)
    treated_weight = lowercase(treated_weight)
    treated_weight_options = ["treated", "population"]
    if !(treated_weight in treated_weight_options)
        error("'treated_weight' must be one of $(treated_weight_options)")
    end

    # Round nperm
    nperm = Int(round(nperm))
    if nperm < 1
        error("'nperm' must be a positive integer > 0.")
    end

    # Round ri_kernel_cutoff
    ri_kernel_cutoff = Int(round(ri_kernel_cutoff))
    if ri_kernel_cutoff < 0
        error("'ri_kernel_cutoff' must be a non-negative integer.")
    end

    # Create data copy
    data_copy = DataFrame(data; copycols = false)

    # Do validation checks to see if specified columns exist, and that they are able to be processed correctly
    data_copy, treatment_times, covariates = validate_data(data_copy, outcome, state, time, covariates, treatment_times, date_format;
                                               warn_missing_covariates = false)
    treatment_times = !(treatment_times isa AbstractVector) ? [treatment_times] : treatment_times
    if isempty(treatment_times)
        error("'treatment_times' must have at least one entry.")
    end

    # Make sure treated_states is a vector. A single treatment time applies to every treated state;
    # common adoption is the one-cohort case of the same estimator (update spec §3.6).
    treated_states = !(treated_states isa AbstractVector) ? [treated_states] : treated_states
    if length(treatment_times) == 1 && length(treated_states) > 1
        treatment_times = fill(treatment_times[1], length(treated_states))
    end
    if length(treatment_times) != length(treated_states)
        error("'treatment_times' should either be the same length as the 'treated_states' vector or of length 1.")
    end

    # Make sure control states exist. Not-yet-treated states can only serve as controls when
    # there are at least two treatment times.
    unique_states = unique(data_copy[!, state])
    control_states = setdiff(unique_states, treated_states)
    if isempty(control_states) && !use_pre_controls
        error("No control states were found.")
    elseif isempty(control_states) && length(unique(treatment_times)) == 1
        error("No valid control states were found.")
    end

    # Check missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # Ensure the state column is a string or number and that the nonmissingtype(treated_states) == nonmissingtype(state column)
    validate_state_types(treated_states, data_copy, state)
    data_copy.state_71X9yTx = data_copy[!, state]

    # Do some checks after (potentially) dropping rows to make sure that the analysis
    # is still valid
    if (length(unique(data_copy.state_71X9yTx)) == 1)
        error("Not enough non-missing observations.\nFound", length(unique(data_copy.state_71X9yTx)), "states after dropping missing observations.")
    end

    # Validate and convert dates to Date objects, filter data by date range
    data_copy, treatment_times, start_date, end_date, all_times, freq = validate_and_convert_dates(data_copy, time, treatment_times, date_format,
                                                                                                   freq, start_date, end_date)

    # Match dates to the period grid
    data_copy, treatment_times, treated_states, match_to_these_dates, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier,
                                                                                                                    start_date, end_date, treatment_times,
                                                                                                                    treated_states)

    # Check that treated states actually exist
    treated_states, treatment_times = validate_treated_states(treated_states, treatment_times, data_copy)

    # Also need to make sure that start_times < treatment_times < end_times is true for each state
    for i in eachindex(treated_states)
        s = treated_states[i]
        treat_time = treatment_times[i]
        state_dates = data_copy[data_copy.state_71X9yTx .== s, :time_71X9yTx]
        earliest = minimum(state_dates)
        latest = maximum(state_dates)
        if !(earliest < treat_time)
            earliest = string(earliest)
            treat_time = string(treat_time)
            error("For state $s, the earliest date for which there is non-missing value ($earliest) is not strictly less than the treatment time ($treat_time).")
        end
        if !(treat_time <= latest)
            treat_time = string(treat_time)
            latest = string(latest)
            error("For state $s, the treatment time ($treat_time) is greater than the last date ($latest) for which there is a non-missing observation.")
        end
    end

    # Ensure state column is a string (as are treated_states)
    data_copy, treated_states = validate_string_treated_states(data_copy; treated_states = treated_states)

    # Once the date matching is done, keep a Date copy of the time column and convert
    # time_71X9yTx to a string for processing in `categorical()`
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)

    # Create dummies for each time and state interaction
    data_copy.state_time = categorical(data_copy.state_71X9yTx .* "0IQR7q6Wei7Ejp4e" .* data_copy.time_71X9yTx)

    # Convert factor covariates into multiple numeric dummy variable columns
    data_copy, covariates_to_include = process_covariates(covariates, data_copy, ref)

    # Force outcome to float64 to speed up regression (runs faster if <:Number rather than <:Union{Number, Missing})
    data_copy.outcome_71X9yTx = convert(Vector{Float64}, data_copy.outcome_71X9yTx)

    # The edgecase standard errors are not available for ccc = "add" (its first stage is rank deficient)
    if edgecase == true && ccc == "add"
        @warn "'edgecase' standard errors are not available for the 'ccc = add' option."
        edgecase = false
    end

    # Overwrite truejack depending on the ccc arg
    if ccc in ["int", "state"]
        truejack = false
    end

    # Run the fixed effects model and get back the dataframe of means (or means residualized by covariates) for each period at each state
    lambda_df, vcov_lambda, covariate_fit = iterative_demean(data_copy, ccc, covariates_to_include, hc, edgecase,
                                                             time_to_index = time_to_index, treatment_times = treatment_times, match_to_these_dates = match_to_these_dates,
                                                             treated_states = treated_states, use_pre_controls = use_pre_controls)

    # Integrated long differences (DID-INT eqs. 46-47; update spec §3.2, D5). The covariate
    # term they add cancels in ATT(s,t) for every ccc except "int": it is zero for "hom" and
    # "state", and for "time" and "add" the treated state and its controls carry the same
    # time slopes at the same covariate means. So it is only applied under "int" with covariates.
    integration = (ccc == "int" && !isempty(covariates_to_include)) ? covariate_fit : nothing

    # Long differences: diff_df holds the rows of the ATT(s,t) regressions, ri_diff_df the extra
    # differences randomization inference needs
    diff_df, ri_diff_df = construct_diff_df(lambda_df, treated_states, treatment_times, time_to_index, match_to_these_dates,
                                             use_pre_controls, data_copy)

    # Do check to make sure that each state has each (g,t) group for RI
    randomize = check_prior_to_ri(diff_df, ri_diff_df)

    # Show period length in results
    if isnothing(date_format) || date_format == "yyyy"
        date_format = assume_date_format(period)
    end
    period = string(period)
    start_date = parse_date_to_string_didint(start_date, date_format)
    end_date = parse_date_to_string_didint(maximum(data_copy.time_dmG5fpM), date_format)

    # ATT(s,t) for every treated state and post period (DID-INT eq. 48), then the sub-aggregates
    # of the chosen agg and the aggregate ATT, all weighted means of ATT(s,t)'s (update spec §3.3-§5)
    # Policies of the treated states (update spec §6). They only enter the policy aggregations:
    # controls do not depend on policy (§6.2), so every other agg ignores them.
    policy = (labels = Dict(s => (is_policy_agg(agg) ? state_policies[s] : DEFAULT_POLICY) for s in treated_states),
              pooled = policy_aggregation, weighting = policy_weighting)

    att_st = build_att_st(diff_df, weighting, hc, vcov_lambda; integration = integration, policies = policy.labels)
    results = sub_aggregate_results(att_st, agg, weighting, treated_weight, hc; pooled = policy.pooled)
    if is_policy_agg(agg)
        # Each policy's ATT, and with policy_aggregation the pooled ATT across policies (§6.3-§6.4)
        results = policy_att_results(results, att_st, agg, weighting, hc, policy)
    else
        results = aggregate_att(results, agg, weighting, hc)
    end

    # Jackknife and randomization inference recompute the same weighted means on each deletion / draw
    results = jackknife_procedure(results, att_st, diff_df, agg, weighting, treated_weight, use_pre_controls,
                                  truejack, data_copy, ccc, covariates_to_include, treated_states, treatment_times,
                                  time_to_index, match_to_these_dates; integration = integration, policy = policy)
    if randomize
        results = randomization_inference(results, att_st, vcat(diff_df, ri_diff_df), agg, weighting, treated_weight,
                                          use_pre_controls, treated_states, treatment_times, nperm, ri_kernel_cutoff, seed;
                                          integration = integration, policy = policy)
    end

    # Show the period length and date range in the results
    results.period .= period
    results.start_date .= start_date
    results.end_date .= end_date
    results = wrapper_check(results, wrapper)
    return results
end

### Second stage: ATT(s,t) table and its aggregations (update spec §3.3-§5, §11.1)

# For each aggregation: the suffix of its sub-aggregate columns (att_cohort, se_att_cohort, ...)
# and its key columns as (results column, ATT(s,t) table column, element type). The results
# columns are the public output names that the Stata/R wrappers rely on.
const SUB_AGG_SPECS = Dict(
    "cohort"   => (suffix = "att_cohort", keys = [(:treatment_time, :gvar, Date)]),
    "simple"   => (suffix = "att_gt", keys = [(:r1, :r1, Date), (:time, :t, Date), (:gvar, :gvar, Date)]),
    "state"    => (suffix = "att_s", keys = [(:state, :state, String)]),
    "sgt"      => (suffix = "att_sgt", keys = [(:state, :state, String), (:gvar, :gvar, Date), (:t, :t, Date)]),
    "time"     => (suffix = "att_t", keys = [(:periods_post_treat, :periods_post_treat, Float64)]),
    "calendar" => (suffix = "att_calendar", keys = [(:t, :t, Date)]),
    # Policy aggregations (update spec §6.3): the same groups within each policy. "policy" has
    # one sub-aggregate per policy, the policy's ATT itself.
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

function init_results(agg::AbstractString, nrows::Integer; pooled::Bool = false)

    # Results DataFrame with the public column names in their existing order: key columns,
    # sub-aggregate ATT, top-level ATT and its inference, sub-aggregate inference, weights.
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
    # How the RI p-values were computed (update spec §8.2); filled by randomization_inference
    results.ri_procedure = Vector{Union{Missing, String}}(missing, nrows)
    return results
end

function step3_weights(temp::DataFrame, weighting::AbstractString)

    # Weights of the long differences in an ATT(s,t) regression: observation counts n under
    # weighting = "diff" or "both", otherwise unweighted. The treated row's weight cancels.
    if in(weighting, ["diff", "both"])
        W = convert(Vector{Float64}, temp.n)
        W ./= sum(W)
        return W
    else
        return fill(nothing, nrow(temp))
    end
end

function slope_change_term(fit::CovariateFit, source_state, t, r1, target_state)

    # Covariate term of an integrated long difference (DID-INT eqs. 46-47; update spec §3.2):
    # (gamma_{source,t} - gamma_{source,r1})' xbar_{target,t}, the change in the source state's
    # slopes evaluated at the target (treated) state's period-t covariate means. Because f(X) is
    # linear, integrating over the target's covariate distribution reduces to its mean.
    return dot(fit.slopes[(source_state, t)] .- fit.slopes[(source_state, r1)], fit.xbar[(target_state, t)])
end

function build_att_st(diff_df::DataFrame, weighting::AbstractString, hc::AbstractString, vcov_lambda;
                      integration::Union{CovariateFit, Nothing} = nothing,
                      policies::Union{AbstractDict, Nothing} = nothing)

    # ATT(s,t) for each treated state s and post period t (DID-INT eq. 48; update spec §3.3):
    # the treat coefficient of the weighted regression of the long differences on [1, treat],
    # using s's long difference and those of the controls in the same (t, r1) cell (§3.4).
    # With `integration`, every long difference in the regression gains its covariate term
    # evaluated at s's covariate means (eqs. 46-47), so a control's long difference depends
    # on the treated state it is compared with.
    # With one treated row its leverage is 1, so HC2-HC4 standard errors are NaN (D17). With
    # exactly two rows the regression SE is degenerate and the edge-case SE is used instead.
    # An (s,t) without a control or whose regression cannot be estimated is dropped.
    # `policies` gives each treated state's policy (one common policy if not given); the
    # controls do not depend on it (§3.4, §6.2).
    treated = sort(diff_df[diff_df.treat .== 1, :], [:state, :treated_time, :t])
    nrows = nrow(treated)
    float_column() = Vector{Union{Missing, Float64}}(missing, nrows)
    att_st = DataFrame(state = treated.state, gvar = treated.treated_time, t = treated.t, r1 = treated.r1,
                       periods_post_treat = treated.time_since_treatment,
                       policy = String[isnothing(policies) ? DEFAULT_POLICY : policies[s] for s in treated.state],
                       att = float_column(), se = float_column(), pval = float_column(),
                       n_t = treated.n_t, n = treated.n, n_controls = zeros(Int, nrows))
    for i in 1:nrows
        controls = diff_df[(diff_df.treat .== 0) .& (diff_df.t .== treated.t[i]) .& (diff_df.r1 .== treated.r1[i]), :]
        att_st.n_controls[i] = nrow(controls)
        nrow(controls) == 0 && continue
        temp = vcat(treated[i:i, :], controls)
        X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
        Y = convert(Vector{Float64}, temp.diff)
        if !isnothing(integration)
            Y .+= [slope_change_term(integration, source, treated.t[i], treated.r1[i], treated.state[i]) for source in temp.state]
        end
        W = step3_weights(temp, weighting)
        result_dict = final_regression_results(X, Y, W = W, hc = hc)
        att_st.att[i] = result_dict["beta_hat"]
        if length(Y) != 2
            att_st.se[i] = result_dict["beta_hat_se"]
        elseif isnothing(integration)
            att_st.se[i] = edge_case_se(temp.lambda_idx, vcov_lambda)
        else
            att_st.se[i] = edge_case_se_integrated(integration, treated.state[i], controls.state[1], treated.t[i], treated.r1[i], hc)
        end
        att_st.pval[i] = result_dict["pval_att"]
    end
    att_st = att_st[(att_st.n_controls .> 0) .& .!ismissing.(att_st.att), :]
    if nrow(att_st) == 0
        error("Could not estimate a single ATT(s,t): no treated state has a long difference with a valid control in the same period.")
    end
    return att_st
end

function treated_unit_weights(att_st::AbstractDataFrame, weighting::AbstractString, treated_weight::AbstractString)

    # Level 1 weight u of each ATT(s,t) within its sub-aggregate (update spec §4.1):
    #   weighting   treated_weight = "population"   treated_weight = "treated"
    #   "none"      1                               1
    #   "diff"      n = n_{s,t} + n_{s,r1}          1
    #   "att"       1                               n_t = n_{s,t}
    #   "both"      n                               n_t
    # "treated" makes a sub-aggregate the effect on the average treated observation
    # (Which Policy, §3.1, pi_{s|g,t}).
    if weighting == "none"
        return ones(nrow(att_st))
    elseif treated_weight == "population"
        return weighting in ["diff", "both"] ? Float64.(att_st.n) : ones(nrow(att_st))
    else
        return weighting in ["att", "both"] ? Float64.(att_st.n_t) : ones(nrow(att_st))
    end
end

function sub_aggregate_results(att_st::DataFrame, agg::AbstractString, weighting::AbstractString,
                               treated_weight::AbstractString, hc::AbstractString; pooled::Bool = false)

    # Each sub-aggregate is the u-weighted mean of the ATT(s,t)'s in its group (update spec §5.1;
    # for the policy aggregations, the same groups within each policy, §6.3).
    # Its HC SE comes from the weighted regression of those ATT(s,t)'s on a constant (D6, §5.2;
    # this treats them as independent, which can understate the SE when they share controls).
    # A group with a single ATT(s,t) reports that ATT(s,t)'s own SE and p-value (A1). The Level 2
    # weight of a group is the sum of its n_t under weighting = "att"/"both" (§4.2).
    spec = SUB_AGG_SPECS[agg]
    group_cols = group_columns(agg)
    u = treated_unit_weights(att_st, weighting, treated_weight)

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

function aggregate_att(results::DataFrame, agg::AbstractString, weighting::AbstractString, hc::AbstractString)

    # Aggregate ATT: weighted regression of the sub-aggregate ATTs on a constant, with the
    # Level 2 weights (update spec §4.2, D4, D6). A single sub-aggregate is the aggregate, so it
    # reports that sub-aggregate's own SE and p-value (as A1 does for a single ATT(s,t)).
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

    # Weight omega_p of each policy in the pooled ATT (update spec §6.4, D15), from a table of
    # ATT(s,t)'s: "equal" 1, "treated_obs" the policy's treated observations (sum of n_t),
    # "treated_cells" its number of ATT(s,t)'s, "treated_states" its number of treated states.
    # Normalized to sum to 1.
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

    # Top-level estimates of the policy aggregations (update spec §6.3-§6.4). Each policy's ATT_p
    # is the Level 2 weighted mean of that policy's sub-aggregates, with its HC SE from the
    # weighted regression of those sub-aggregates on a constant, as aggregate_att does across
    # all sub-aggregates (a single sub-aggregate reports its own SE). It is written on the first
    # row of the policy. For agg = "policy" the sub-aggregate is ATT_p itself.
    suffix = SUB_AGG_SPECS[agg].suffix
    policies = unique(results.policy)

    # Level 2 weights (§4.2), normalized within each policy. agg = "policy" has no sub-aggregates
    # below the policy's ATT, so it has none.
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

    # Pooled ATT across policies: ATT_pool = sum_p omega_p ATT_p (D15, D16), with its HC SE from
    # the weighted regression of the ATT_p's on a constant (number of policies - 1 degrees of
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
                                                   W = copy(omega), hc = hc)
            results.att_pool[1] = result_dict["beta_hat"]
            results.se_att_pool[1] = result_dict["beta_hat_se"]
            results.pval_att_pool[1] = result_dict["pval_att"]
        end
    end
    return results
end

### Helpers shared with the batched engine (src/batched_engine.jl), which recomputes the estimates
### for the jackknife deletions and the randomization inference draws

function long_difference_lookup(long_diffs::AbstractDataFrame)

    # (state, t, r1) => (diff, n, n_t) for every non-missing long difference
    lookup = Dict{Tuple{String, Date, Date}, NamedTuple{(:diff, :n, :n_t), Tuple{Float64, Float64, Float64}}}()
    for row in eachrow(long_diffs)
        ismissing(row.diff) && continue
        lookup[(row.state, row.t, row.r1)] = (diff = Float64(row.diff), n = Float64(row.n), n_t = Float64(row.n_t))
    end
    return lookup
end

# The (cohort, period) cells of the observed ATT(s,t)'s, with their base period and event time
att_st_cells(att_st::DataFrame) = sort(unique(att_st[:, [:gvar, :t, :r1, :periods_post_treat]]), [:gvar, :t])

function result_key(results::DataFrame, i::Integer, agg::AbstractString)

    # Key tuple of results row i, with the event time as an Int as in the ATT(s,t) table
    return Tuple(coltype === Float64 ? Int(results[i, col]) : results[i, col]
                 for (col, _, coltype) in SUB_AGG_SPECS[agg].keys)
end

### The following functions are parent functions within the didint() function:
function scale_weights_final(results::DataFrame, weighting::AbstractString)
    
    # Scale weights and convert to either Float64 or Nothing
    if weighting in ["att", "both"]
        results.weights ./= sum(results.weights)
        results.weights = convert(Vector{Float64}, results.weights)
    elseif weighting in ["diff", "none"]
        results.weights = convert(Vector{Nothing}, results.weights)
    end
    return results
end

function final_regression_results(X::Matrix, Y::Vector; W::Vector = [nothing], hc::AbstractString = "hc1")

    n = length(Y)
    # Check for and remove missing values
    # Build mask for each component
    x_valid = [!any(ismissing, X[i, :]) for i in 1:n]
    y_valid = .!ismissing.(Y)

    if eltype(W) <: Nothing
        valid_mask = x_valid .& y_valid
    else
        w_valid = .!ismissing.(W)
        valid_mask = x_valid .& y_valid .& w_valid
    end

    # Filter out missing observations
    X = X[valid_mask, :]
    Y = Y[valid_mask]
    if eltype(W) <: Number
        W = W[valid_mask]
    end
    
    # Check if we have enough observations left
    ncolx = size(X, 2)
    n = length(Y)
    if n < ncolx
        @warn "Insufficient observations after dropping missing values (n = $n, k = $ncolx)"
        return Dict("beta_hat" => missing, "beta_hat_se" => missing, 
                   "pval_att" => missing, "beta_hat_se_jknife" => missing, 
                   "pval_att_jknife" => missing)
    end

    beta_hat = nothing
    beta_hat_cov = nothing

    # Run OLS (normally if weights aren't provided, and scale (X,Y) -> (Xw,Yw) otherwise)
    if eltype(W) <: Nothing
        beta_hat = safe_solve(X, Y)
        if ismissing(beta_hat)
            return Dict("beta_hat" => missing, "beta_hat_se" => missing, "pval_att" => missing,
                       "beta_hat_se_jknife" => missing, "pval_att_jknife" => missing)
        end
        if n > ncolx
            resid = Y - X * beta_hat
            beta_hat_cov = compute_hc_covariance(X, resid, hc)
        end
    elseif eltype(W) <: Number
        W ./= sum(W)
        sw = sqrt.(W)           
        Xw = X .* sw            
        Yw = Y .* sw
        beta_hat = safe_solve(Xw, Yw)
        if ismissing(beta_hat)
            return Dict("beta_hat" => missing, "beta_hat_se" => missing, "pval_att" => missing,
                       "beta_hat_se_jknife" => missing, "pval_att_jknife" => missing)
        end
        if n > ncolx
            resid_w = Yw - Xw * beta_hat
            beta_hat_cov = compute_hc_covariance(Xw, resid_w, hc)
        end
    end 

    if n > ncolx
        beta_hat_var = diag(beta_hat_cov)
        beta_hat_se = sqrt(beta_hat_var[ncolx]) 
    else 
        beta_hat_se = missing
    end
    dof = n - ncolx
    pval_att = dof > 0 && !ismissing(beta_hat_se) ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se))) : missing 
    result_dict = Dict("beta_hat" => beta_hat[ncolx], "beta_hat_se" => beta_hat_se, "pval_att" => pval_att)
    return result_dict
end

function custom_sort_order(s)
    parsed = tryparse(Int, s)
    if isnothing(parsed)
        (0, lowercase(s))
    else
        (1, parsed)
    end 
end

function check_prior_to_ri(diff_df, ri_diff_df)
    
    randomize = true

    # Get all required (t, r1) combinations and check for each state, each pair exists in the combined data
    combined_df = vcat(diff_df, ri_diff_df)
    all_states = unique(combined_df.state)
    required_pairs = unique(select(combined_df, :t, :r1))
    for state in all_states
        state_data = combined_df[combined_df.state .== state, :]
        state_pairs = unique(select(state_data, :t, :r1))
        missing_pairs = antijoin(required_pairs, state_pairs, on=[:t, :r1])
        if nrow(missing_pairs) > 0
            randomize = false
            break
        end
    end
    if !randomize
        @warn "Randomization inference cannot be performed.\nDid not find all the non-missing values for differences required to implement randomization inference procedure."
    end

    return randomize

end

# Policy options when none are given: one common policy, no pooling
default_policy_options(treated_states) = (labels = Dict(s => DEFAULT_POLICY for s in treated_states),
                                          pooled = false, weighting = "equal")

# Number of draws that randomization inference evaluates at a time
const RI_BATCH_SIZE = 500

function randomization_inference(results::DataFrame, att_st::DataFrame, long_diffs::DataFrame, agg::AbstractString,
                                 weighting::AbstractString, treated_weight::AbstractString, use_pre_controls::Bool,
                                 treated_states, treatment_times, nperm::Int, ri_kernel_cutoff::Int, seed::Number;
                                 integration::Union{CovariateFit, Nothing} = nothing,
                                 policy = default_policy_options(treated_states))

    # Randomization inference (MacKinnon & Webb 2020): give the observed (treatment time, policy)
    # labels to other states, recompute the same estimates for each alternative assignment, and
    # compare them with the observed estimates (update spec §5.2, §7, §8). An RI p-value is
    # reported at every level of the output: every sub-aggregate, the aggregate ATT or each
    # policy's ATT, and the pooled ATT. long_diffs holds every state's long difference in every
    # cell. The first stage does not depend on the assignment, so the batched engine evaluates
    # all draws from the same lambdas, slopes and covariate means (§7.1, §7.3).

    # PART ONE: THE ALTERNATIVE ASSIGNMENTS OF THE (TREATMENT TIME, POLICY) LABELS
    # All of them if there are at most `nperm`, otherwise `nperm` distinct random ones (§7.2)
    engine = build_engine(long_diffs, att_st_cells(att_st), unique(long_diffs.state), weighting, treated_weight;
                          integration = integration)
    original_treated = unique(long_diffs[long_diffs.treat .== 1, [:state, :treated_time]])
    policy_names, policy_code = policy_coding(results, agg, policy.labels)
    L_observed, P_observed = assignment_row(engine.states, Dict(zip(original_treated.state, original_treated.treated_time)),
                                            policy.labels, policy_code)
    rng = Xoshiro(seed)
    L, P, exact = ri_assignment_draws(L_observed, P_observed, nperm, rng)
    if size(L, 1) < nperm
        @warn "'nperm' was set to $nperm but only $(size(L, 1)) unique permutations exist. \n
                Setting 'nperm' to $(size(L, 1))."
        nperm = size(L, 1)
    end
    if nperm < 399
        @warn "'nperm' is less than 399!"
    end

    # PART TWO: COMPUTE THE PLACEBO SUB-AGGREGATE AND TOP-LEVEL ATTS FOR EVERY DRAW
    # att_ri holds the placebo aggregate ATT, or for the policy aggregations the pooled ATT;
    # att_ri_policy holds each policy's placebo ATT_p. The statistic is the raw estimate, and a
    # placebo sub-aggregate is matched to the observed one by its key (§7.4). "state" and "sgt"
    # keep the stand-in procedure of §7.5.
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
    # The kernel-smoothed p-value (§8.2) when at most ri_kernel_cutoff alternatives were used,
    # whether they were all enumerated or sampled (author's choice; the spec smooths only
    # enumerated ones), otherwise the ordinary p-value (§8.1). With a single alternative there
    # is nothing to smooth.
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
    # nperm is the number of alternatives used; ri_procedure records how the p-values were
    # computed: kernel-smoothed or ordinary, from all alternatives (exact) or a sample of them
    results.nperm[1] = nperm
    results.ri_procedure .= (use_kernel ? "kernel" : "ordinary") * (exact ? "_exact" : "_sampled")
    return results
end

## The following functions compute the randomization inference p-values (update spec §8)

# Relative tolerance for ties between a placebo estimate and the observed estimate
const RI_TIE_RTOL = sqrt(eps(Float64))

function ri_pvalue(att_ri::AbstractVector, observed::Real)
    # RI p-value p = (1 + #{|T*| >= |T|}) / (S + 1) over the S finite placebo estimates T*,
    # the more conservative of the two versions in MacKinnon & Webb (2020) (update spec D8, §8.1).
    # Ties count as hits. A placebo assignment that leaves a group's treated states and
    # its controls unchanged (common with use_pre_controls = true) reproduces the observed
    # estimate exactly in exact arithmetic, but the two are computed along different paths
    # and can differ by rounding (~1e-16), which also varies with compiler flags such as
    # --check-bounds. Treating |T*| within a relative RI_TIE_RTOL of |T| as a tie keeps
    # the count deterministic. Returns missing if there is no finite placebo estimate.
    finite_att_ri = filter(isfinite, att_ri)
    n_placebo = length(finite_att_ri)
    n_placebo == 0 && return missing
    threshold = abs(observed) * (1 - RI_TIE_RTOL)
    n_hits = count(x -> abs(x) >= threshold, finite_att_ri)
    return (1 + n_hits) / (n_placebo + 1)
end

function ri_pvalue_kernel(att_ri::AbstractVector, observed::Real)
    # Kernel-smoothed RI p-value (Racine & MacKinnon 2007, "Inference via kernel smoothing of
    # bootstrap P values"; update spec D10, §8.2):
    #   p = 1 - (1/S) sum_j Phi((|T| - |T*_j|) / h),   h = 1.575 sd(T*) S^(-4/9),
    # over the S finite placebo estimates T*, with Phi the standard normal CDF and sd the sample
    # standard deviation of the signed T*. With few alternatives the ordinary p-value can only
    # take a few values; smoothing the indicator |T*| >= |T| lets it move continuously with T.
    # Falls back to the ordinary p-value if S < 2 or sd is zero or not finite.
    finite_att_ri = filter(isfinite, att_ri)
    n_placebo = length(finite_att_ri)
    n_placebo == 0 && return missing
    sd = n_placebo >= 2 ? std(finite_att_ri) : NaN
    if !(isfinite(sd) && sd > 0)
        return ri_pvalue(att_ri, observed)
    end
    h = 1.575 * sd * n_placebo^(-4 / 9)
    return 1 - mean(cdf.(Normal(), (abs(observed) .- abs.(finite_att_ri)) ./ h))
end

## The following functions are for the jackknife procedure
function jackknife_procedure(results::DataFrame, att_st::DataFrame, diff_df::DataFrame, agg::AbstractString,
                             weighting::AbstractString, treated_weight::AbstractString, use_pre_controls::Bool,
                             truejack::Bool, data_copy, ccc, covariates_to_include, treated_states, treatment_times,
                             time_to_index, match_to_these_dates;
                             integration::Union{CovariateFit, Nothing} = nothing,
                             policy = default_policy_options(treated_states))

    # Delete-one-state jackknife (CRV3, centered at the full-sample estimate; update spec §9.1).
    # Each state is deleted as a treated state and as a control, and the same weighted means are
    # recomputed, at every level of the output (sub-aggregates, the aggregate ATT or each
    # policy's ATT, and the pooled ATT). Deleting a state leaves the other states' lambdas and
    # slopes unchanged under ccc = "int" and "state", or is taken to under truejack = false, so
    # the batched engine evaluates all deletions in one pass (§9.2). With truejack = true (ccc
    # not "int" or "state") the first stage is re-estimated without the state, since its lambdas
    # then depend on the other states. Under "int", the only ccc that integrates, each cell's
    # slopes and covariate means are unaffected by deletions.

    # A jackknife needs at least 2 treated and 2 control states; otherwise every jackknife
    # SE and p-value stays missing
    n_treated = length(unique(diff_df[diff_df.treat .== 1, :state]))
    n_control = length(unique(diff_df[diff_df.treat .== 0, :state]))
    if n_control < 2 || n_treated < 2
        @warn "Need at least 2 treated states and 2 control states to do jackknife procedure."
        return results
    end

    # Delete-one estimates, one row per deleted state. top_jacks holds the aggregate ATT, or for
    # the policy aggregations the pooled ATT; policy_jacks holds each policy's ATT_p.
    cells = att_st_cells(att_st)
    observed_labels = Dict(zip(treated_states, treatment_times))
    all_states = unique(diff_df.state)
    n_states = length(all_states)
    policy_names, policy_code = policy_coding(results, agg, policy.labels)
    if !truejack
        # Batched fast jackknife (update spec §9.2): every row of the assignment matrix is the
        # observed assignment, and row h switches state h off as a treated state and as a
        # control. The weights renormalize automatically.
        engine = build_engine(diff_df, cells, all_states, weighting, treated_weight; integration = integration)
        L_observed, P_observed = assignment_row(all_states, observed_labels, policy.labels, policy_code)
        est = engine_estimates(engine, results, agg, weighting, use_pre_controls,
                               repeat(L_observed', n_states), repeat(P_observed', n_states);
                               active = .!Matrix{Bool}(I, n_states, n_states),
                               policy_names = policy_names, policy_weighting = policy.weighting)
        sub, att_p, top = est.sub, est.att_p, est.top
    else
        # Refit path: re-estimate the first stage without state h, rebuild the long differences,
        # and evaluate the observed assignment of the remaining states
        sub = fill(NaN, n_states, nrow(results))
        att_p = fill(NaN, n_states, length(policy_names))
        top = fill(NaN, n_states)
        for (h, without) in enumerate(all_states)
            lambda_h, _, fit_h = iterative_demean(data_copy[data_copy.state_71X9yTx .!= without, :], ccc, covariates_to_include, "hc0", false)
            diff_h = construct_diff_df(lambda_h, treated_states, treatment_times, time_to_index, match_to_these_dates,
                                       use_pre_controls, data_copy; jackknife = true, without = without)
            remaining = filter(!=(without), all_states)
            engine_h = build_engine(diff_h, cells, remaining, weighting, treated_weight;
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

    # Treated states behind a set of ATT(s,t)'s, and control states in their cells
    control_rows = diff_df[diff_df.treat .== 0, :]
    function state_counts(in_group)
        group_cells = Set(zip(att_st.t[in_group], att_st.r1[in_group]))
        n_group_control = length(unique(control_rows.state[[cell in group_cells for cell in zip(control_rows.t, control_rows.r1)]]))
        return length(unique(att_st.state[in_group])), n_group_control
    end

    # Sub-aggregate jackknife SEs and p-values. A sub-aggregate gets none unless its ATT(s,t)'s
    # involve at least 2 treated states and their cells at least 2 control states, and every
    # deletion leaves it defined. The single-state sub-aggregates of "state" and "sgt" never do.
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

    # Top-level jackknife SEs and p-values, over the deletions that leave the estimate defined
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
    # the pool defined over the other policies; the SE is still reported (update spec §9.1).
    if policy.pooled
        if any(ismissing(policy_jacks[h, j]) && !ismissing(top_jacks[h]) for h in 1:length(all_states), j in 1:length(policy_names))
            @warn "Some jackknife deletions remove a whole policy. The jackknife standard error of the pooled ATT then compares pools over different sets of policies."
        end
        results.jknifese_att_pool[1], results.jknifepval_att_pool[1] = jackknife_se_pval(collect(skipmissing(top_jacks)), results.att_pool[1])
    end
    return results
end

function jackknife_se_pval(estimates::AbstractVector, full_sample::Real)

    # CRV3 jackknife SE sqrt((n-1)/n * sum_h (theta_h - theta)^2) and its two-sided p-value
    # from a t(n-1) distribution (update spec §9.1; MacKinnon, Nielsen & Webb 2023)
    n = length(estimates)
    se = n >= 2 ? sqrt(sum((estimates .- full_sample).^2) * ((n - 1) / n)) : missing
    pval = n - 1 > 0 && !ismissing(se) ? 2 * (1 - cdf(TDist(n - 1), abs(full_sample / se))) : missing
    return se, pval
end

function construct_diff_df(lambda_df, treated_states, treatment_times, time_to_index, match_to_these_dates, use_pre_controls,
                           data_copy; jackknife = false, without = nothing)

    if jackknife
        mask = treated_states .!= without
        treated_states = treated_states[mask]
        treatment_times = treatment_times[mask]
    end

    # Initialize diff_df
    unique_states = unique(lambda_df.state)
    diff_df = DataFrame(state = String[], treated_time = Date[],
                        t = Date[], r1 = Date[], diff = Float64[],
                        treat = Int[], n = Int[], n_t = Int[],
                        lambda_idx = Union{Tuple{Int,Int}, Missing}[])

    for i in eachindex(treated_states)
        idx = time_to_index[treatment_times[i]]
        one_period_prior_treatment = match_to_these_dates[idx - 1]
        temp = lambda_df[(lambda_df.state .== treated_states[i]) .& (lambda_df.time .>= one_period_prior_treatment), :]
        lambda_r1 = temp[temp.time .== one_period_prior_treatment, [:lambda, :lambda_index]]
        empty_check = isempty(lambda_r1.lambda)
        lambda_r1_var_idx =  empty_check ? missing : lambda_r1.lambda_index[1]
        lambda_r1 = empty_check ? missing : lambda_r1.lambda[1]
        temp = temp[temp.time .> one_period_prior_treatment, :]
        sort!(temp, :time)
        years = temp.time
        diffs = Vector{Union{Float64, Missing}}(undef, nrow(temp))
        lambda_indices = Vector{Union{Tuple{Int,Int}, Missing}}(undef, nrow(temp))
        n = Vector{Int}(undef, nrow(temp))
        n_t = Vector{Int}(undef, nrow(temp))
        for j in eachindex(diffs)
            lambda_j = temp[j, [:lambda, :lambda_index]]
            empty_check = isempty(lambda_j.lambda)
            lambda_j_var_idx = empty_check ? missing : lambda_j.lambda_index[1]
            lambda_j = empty_check ? missing : lambda_j.lambda[1]
            diffs[j] = lambda_j - lambda_r1
            # (r1, t) order, as for control rows: edge_case_se assumes both rows of a
            # two-row regression list their lambdas in the same order
            lambda_indices[j] = ismissing(lambda_j_var_idx) || ismissing(lambda_r1_var_idx) ? missing : (Int(lambda_r1_var_idx), Int(lambda_j_var_idx))
            post = temp[j, "time"]
            pre = one_period_prior_treatment
            n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([post, pre])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            n_t[j] = sum((data_copy.time_dmG5fpM .== post) .&& (data_copy.state_71X9yTx .== treated_states[i]))
        end 
        temp_df = DataFrame(state = treated_states[i], treated_time = treatment_times[i],
                            t = years, r1 = one_period_prior_treatment, diff = diffs, treat = 1,
                            n = n, n_t = n_t, lambda_idx = lambda_indices)
        diff_df = vcat(diff_df, temp_df) 
    end
    
    # Compute diffs for control states
    unique_diffs = unique(select(diff_df, :t, :r1, :treated_time))
    control_states = setdiff(unique_states, treated_states)     
    for i in eachindex(control_states)
        temp = lambda_df[(lambda_df.state .== control_states[i]), :]
        diffs = Vector{Union{Float64, Missing}}(undef, nrow(unique_diffs))
        lambda_indices = Vector{Union{Tuple{Int,Int}, Missing}}(undef, nrow(unique_diffs))
        n = Vector{Int}(undef, nrow(unique_diffs))
        n_t = Vector{Int}(undef, nrow(unique_diffs))
        for j in 1:nrow(unique_diffs)
            t = unique_diffs[j,"t"]
            r1 = unique_diffs[j,"r1"]

            lambda_t = temp[temp.time .== t, [:lambda, :lambda_index]]
            empty_check = isempty(lambda_t.lambda)
            lambda_t_idx = empty_check ? missing : lambda_t.lambda_index[1]
            lambda_t = empty_check ? missing : lambda_t.lambda[1]

            lambda_r1 = temp[temp.time .== r1, [:lambda, :lambda_index]]
            empty_check = isempty(lambda_r1.lambda)
            lambda_r1_idx = empty_check ? missing : lambda_r1.lambda_index[1]
            lambda_r1 = empty_check ? missing : lambda_r1.lambda[1]

            diffs[j] = lambda_t - lambda_r1
            lambda_indices[j] = ismissing(lambda_r1_idx) || ismissing(lambda_t_idx) ? missing : (Int(lambda_r1_idx), Int(lambda_t_idx))

            n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== control_states[i]))
            n_t[j] = sum((data_copy.time_dmG5fpM .== t) .&& (data_copy.state_71X9yTx .== control_states[i]))
        end 
        trtd_time = [match_to_these_dates[findfirst(==(t), match_to_these_dates) + 1] for t in unique_diffs.r1]
        temp_df = DataFrame(state = control_states[i], treated_time = trtd_time,
                            t = unique_diffs.t, r1 = unique_diffs.r1, diff = diffs, treat = 0,
                            n = n, n_t = n_t, lambda_idx = lambda_indices)
        diff_df = vcat(diff_df, temp_df)
    end

    # Add the periods since treatment column to diff_df
    ordered_t = match_to_these_dates   
    t_period = indexin(diff_df.t, ordered_t)             
    treated_time_period = indexin(diff_df.treated_time, ordered_t) 
    diff_df.time_since_treatment = t_period .- treated_time_period 

    # Compute diffs for randomization inference
    if !jackknife || use_pre_controls
        ri_diff_df = DataFrame(state = String[], treated_time = Date[],
                               t = Date[], r1 = Date[], diff = Float64[],
                               treat = Int[], n = Int[], n_t = Int[],
                               lambda_idx = Union{Tuple{Int,Int}, Missing}[])
        for i in eachindex(treated_states)
            ri_diffs = unique(select(diff_df[(diff_df.state .!= treated_states[i]) .& (diff_df.treated_time .!= treatment_times[i]), :], :t, :r1))
            temp = lambda_df[(lambda_df.state .== treated_states[i]), :]
            diffs = Vector{Union{Float64, Missing}}(undef, nrow(ri_diffs))
            lambda_indices = Vector{Union{Tuple{Int,Int}, Missing}}(undef, nrow(ri_diffs))
            n = Vector{Int}(undef, nrow(ri_diffs))
            n_t = Vector{Int}(undef, nrow(ri_diffs))
            for j in 1:nrow(ri_diffs)
                t = ri_diffs[j,"t"]
                r1 = ri_diffs[j,"r1"]

                lambda_t = temp[temp.time .== t, [:lambda, :lambda_index]]
                empty_check = isempty(lambda_t.lambda)
                lambda_t_idx = empty_check ? missing : lambda_t.lambda_index[1]
                lambda_t = empty_check ? missing : lambda_t.lambda[1]

                lambda_r1 = temp[temp.time .== r1, [:lambda, :lambda_index]]
                empty_check = isempty(lambda_r1.lambda)
                lambda_r1_idx = empty_check ? missing : lambda_r1.lambda_index[1]
                lambda_r1 = empty_check ? missing : lambda_r1.lambda[1]

                diffs[j] = lambda_t - lambda_r1
                lambda_indices[j] = ismissing(lambda_r1_idx) || ismissing(lambda_t_idx) ? missing : (Int(lambda_r1_idx), Int(lambda_t_idx))

                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
                n_t[j] = sum((data_copy.time_dmG5fpM .== t) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            end
            trtd_time = [match_to_these_dates[findfirst(==(t), match_to_these_dates) + 1] for t in ri_diffs.r1]
            temp_df = DataFrame(state = treated_states[i], treated_time = trtd_time,
                                t = ri_diffs.t, r1 = ri_diffs.r1, diff = diffs, treat = -1,
                                n = n, n_t = n_t, lambda_idx = lambda_indices)
            ri_diff_df = vcat(ri_diff_df, temp_df)
        end
        t_period = indexin(ri_diff_df.t, ordered_t)             
        treated_time_period = indexin(ri_diff_df.treated_time, ordered_t) 
        ri_diff_df.time_since_treatment = t_period .- treated_time_period 
    end


    # If use_pre_controls is true re-structure the ri_diff_df and diff_df
    # to add the pre-treatment periods from treated silos as control contrasts
    if use_pre_controls
        for (t,s) in zip(treatment_times, treated_states)
            ri_diff_df[ri_diff_df.state .== s .&& ri_diff_df.t .< t, "treat"] .= 0
        end
        diff_df = vcat(diff_df, ri_diff_df[ri_diff_df.treat .== 0,:])
        ri_diff_df = ri_diff_df[ri_diff_df.treat .== -1,:]
        unique_states = unique(diff_df.state)
        control_states = setdiff(unique_states, treated_states)
        if isempty(control_states)
            # If there are no pure control states, do extra processing to remove cells
            # where there are no treat-control pairs
            diff_df = subset(
                        groupby(diff_df, [:t, :r1]),
                        :treat => x -> any(x .== 1) && any(x .== 0)
                      )
            ri_diff_df = semijoin(ri_diff_df, diff_df, on=[:t, :r1])
        end 
    end

    # In the case where the diff_df has missing values for diffs, drop those rows
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

function edge_case_se(lambda_indices, vcov_lambda)

    # The fully expanded term var(ATT) involves multiple cov terms, that under different ccc options
    # are zero. In order to not have a bunch of if else statements and clutter things up, only one 
    # formula is provided here, but note that depending on the ccc selection some (or all for ccc = int)
    # of the cov terms will be zero.
    # Each element of lambda_indices is (index of lambda_r1, index of lambda_t) for one state.
    # The cross-state term -2 Cov(lambda_1[2] - lambda_1[1], lambda_2[2] - lambda_2[1]) only
    # equals -2 Cov(diff_1, diff_2) if both states use this same order.
    # Var(ATT) = c' V_lambda c (update spec Appendix B.2-B.3) assumes that observations in
    # different cells are independent, as in repeated cross-sections. With panel data the same
    # individuals appear in a state's r1 and t cells, their estimates are correlated, and that
    # covariance is omitted here (spec Appendix B.5).

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
function edge_case_se_integrated(fit::CovariateFit, treated, control, t, r1, hc::AbstractString)

    # Edge-case SE of an integrated ATT(s,t) with a single control under ccc = "int" with
    # covariates (update spec A6, Appendix B.4). An OLS fit with an intercept passes through
    # the means, so lambda_{s,t} + gamma_{s,t}' xbar_{s,t} = ybar_{s,t} and
    #   ATT(s,t) = (ybar_{s,t} - delta' xbar_{s,t}) - (lambda_{s,r1} + lambda_{c,t} - lambda_{c,r1}),
    #   delta = gamma_{s,r1} + gamma_{c,t} - gamma_{c,r1}.
    # delta' xbar_{s,t} multiplies slopes estimated from three cells by a covariate mean estimated
    # from the treated cell. xbar_{s,t} is treated as random (it estimates the mean of the
    # distribution eqs. 46-47 integrate over): this is the unconditional variance
    #   Var(zbar_{s,t}) + sum over the 3 other cells k of a' V_k a + tr(V_delta Sigma),
    # with z_i = y_i - delta' x_i over the observations of cell (s,t), a = (1, xbar_{s,t}), V_k the
    # HC covariance of cell k's (lambda, gamma), V_delta the sum of the slope blocks of the three
    # V_k, and Sigma the HC covariance of xbar_{s,t}. With the trace term it equals the variance
    # when the four cells are independent, with estimated inputs.
    # That independence holds for repeated cross-sections. With panel data the same individuals
    # appear in a state's r1 and t cells, so those cells' estimates, and xbar_{s,t} and the
    # state's base-period estimates, are correlated; those covariances are omitted (Appendix B.5).
    other_cells = [(treated, r1), (control, t), (control, r1)]
    if !haskey(fit.cell_obs, (treated, t)) || !all(haskey(fit.cell_vcov, cell) for cell in other_cells)
        return missing
    end

    delta = fit.slopes[(treated, r1)] .+ fit.slopes[(control, t)] .- fit.slopes[(control, r1)]
    y, X = fit.cell_obs[(treated, t)]
    z = y .- X * delta
    var_zbar = hc_mean_covariance(reshape(z, :, 1), hc)[1, 1]

    a = vcat(1.0, fit.xbar[(treated, t)])
    var_other_cells = sum(a' * fit.cell_vcov[cell] * a for cell in other_cells)

    V_delta = sum(fit.cell_vcov[cell][2:end, 2:end] for cell in other_cells)
    Sigma = hc_mean_covariance(X, hc)

    variance = var_zbar + var_other_cells + tr(V_delta * Sigma)
    return isfinite(variance) && variance > 0 ? sqrt(variance) : missing
end
