###########################################################################################
# didint_plot.jl: data for parallel trends and event study plots                          #
#                                                                                         #
#   Section 1: didint_plot(), the main function                                           #
#   Section 2: Parallel trends data                                                       #
#   Section 3: Event study data                                                           #
###########################################################################################


###########################################################################################
# Section 1: didint_plot(), the main function                                             #
###########################################################################################

"""
    didint_plot(
                outcome::Union{AbstractString, Symbol},
                state::Union{AbstractString, Symbol},
                time::Union{AbstractString, Symbol},
                data::DataFrame;
                gvar::Union{AbstractString, Symbol, Nothing} = nothing,
                treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
                treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
                date_format::Union{AbstractString, Nothing} = nothing,
                covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing,
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                ccc::Union{AbstractString, Vector{<:AbstractString}} = "all",
                event::Bool = false,
                weights::Bool = true,
                ci::Number = 0.95,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                hc::Union{AbstractString, Number} = "hc1"
               )

The `didint_plot()` function produces a dataset in a long format that can easily be used for plotting parallel
trends or event study plots.

# Details
The plotted quantities are the first-stage estimates `λ(s,t)` of Karim & Webb (2025, eq. 45): the mean of the outcome
for each state and period, adjusted for covariates as each `ccc` option specifies. For event study plots they are
averaged over the treated states by period relative to treatment. `didint_plot()` is a visual aid for judging parallel
trends under each `ccc` option; it does not reproduce the `didint()` estimates, whose long differences are also
evaluated at the treated states' covariate means under `ccc = "int"`.

The arguments `treated_states` and `treatment_times` should be entered such that the first element in
`treated_states` refers to the state treated at the date entered as the first element in `treatment_times`,
the second element in `treated_states` refers to the state treated at the date entered as the second element
in `treatment_times`, and so on.

# Parameters

## Required Parameters
- `outcome::Union{AbstractString, Symbol}`
    Input the name of the column which identifies the outcome of interest.
- `state::Union{AbstractString, Symbol}`
    Input the name of the column which identifies the state membership of the observation.
- `time::Union{AbstractString, Symbol}`
    Input the name of the column which identifies the date of the observation.
- `data::DataFrame`
    The DataFrame to be used for the analysis.

## Treatment Specification
- `gvar::Union{AbstractString, Symbol, Nothing} = nothing`
    Name of the column which indicates time of first treatment for each state.
- `treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing`
    A vector (or single entry) denoting the associated treatment times of the 'treated_states'.

## Model Specifications
- `ccc::Union{AbstractString, Vector{<:AbstractString}} = "all"`
    Specify which versions of DID-INT should be used.
    Options are either `"all"`, or any combination of: `"none"`, `"hom"`, `"time"`, `"state"`, "`add`", and `"int"`.
    `"none"` uses the unadjusted means of the outcome.
- `covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing`
    A vector of covariates entered as either strings or symbols (or a single covariate string or symbol),
    or, `nothing` (default).
- `ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing`
    A dictionary specifying which category in a categorical variable should be used
    as the reference (baseline) category.

## Event Study Plot
- `event::Bool = false`
   Specify if data should be prepared for an event study plot as opposed to a parallel trends plot.
- `weights::Bool = true`
    Whether to use weighted means when computing event study estimates. If `true`, estimates are computed
    as weighted averages of state-level means for each period relative to treatment, weighted by the states'
    numbers of observations; if `false`, uses simple unweighted averages.
- `treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing`
    A vector of strings or numbers (or a single one) noting the treated state(s), matching the type of the `state` column.
- `ci::Number = 0.95`
    Define the size of confidence bands for the event study plot.
- `hc::Union{AbstractString, Number} = "hc1"`
    Specify which heteroskedasticity-consistent covariance matrix estimator (HCCME) should be used.
    Options are `0`, `1`, `2`, `3`, and `4` (or `"hc0"`, `"hc1"`, `"hc2"`, `"hc3"`, `"hc4"`).

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

# Returns
With `event = false`, a DataFrame with one row per state, period and `ccc` option: the state's mean outcome in that period
(`lambda`), residualized by the specified covariates as the `ccc` option specifies, with the period number (`period`).
Extra rows with only `treat_period` filled give the periods in which treatments start.

With `event = true`, a DataFrame with one row per `ccc` option and period relative to treatment (`time_since_treatment`):
the mean of the treated states' `lambda` values (`y`), its standard error (`se`), the bounds of its confidence interval
(`ci_lower`, `ci_upper`) and the number of treated states behind it (`ngroup`).

# Citations
- Karim & Webb (2025). "Good Controls Gone Bad: Difference-in-Differences with Covariates". https://arxiv.org/abs/2412.14447

"""
function didint_plot(
         outcome::Union{AbstractString, Symbol},
         state::Union{AbstractString, Symbol},
         time::Union{AbstractString, Symbol},
         data;
         gvar::Union{AbstractString, Symbol, Nothing} = nothing,
         treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
         treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
         date_format::Union{AbstractString, Nothing} = nothing,
         covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing,
         ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
         ccc::Union{AbstractString, Vector{<:AbstractString}} = "all",
         event::Bool = false,
         weights::Bool = true,
         ci::Number = 0.95,
         freq::Union{AbstractString, Nothing} = nothing,
         freq_multiplier::Number = 1,
         start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
         end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
         hc::Union{AbstractString, Number} = "hc1",
         wrapper::Union{AbstractString, Nothing} = nothing)

    #--------------------------------------------------------------------------------------
    # Step 1: Check the arguments
    #--------------------------------------------------------------------------------------

    # The HC covariance estimator is only used for the event study
    if event
        hc = hc_checks(hc)
    end
    wrapper = init_wrapper_check(wrapper)
    data_type_check(data)
    date_format = check_date_format(date_format)
    freq_multiplier = check_positive_integer(freq_multiplier, "freq_multiplier")

    # The ccc options to plot
    if !(ccc isa AbstractVector)
        ccc = [ccc]
    end
    ccc = replace.(lowercase.(ccc), r"\s" => "")
    ccc_options = ["hom", "time", "state", "int", "add", "none"]
    if ccc == ["all"]
        ccc = ccc_options
    else
        for c in ccc
            if !(c in ccc_options)
                error("$c is not a valid 'ccc' option. Try any combination of: ", join(ccc_options, ", "))
            end
        end
    end

    # Column names as strings, whether they were given as strings or symbols
    outcome, state, time, gvar = symbol_to_string(outcome, state, time, gvar)
    covariates = isnothing(covariates) ? covariates : string.(covariates)

    # Treated states and their treatment times, from gvar or from treated_states and treatment_times
    treated_states, treatment_times = gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                                              event = event, check = "didint_plot")

    start_date = start_end_date_checks(start_date, "start_date", date_format)
    end_date = start_end_date_checks(end_date, "end_date", date_format)

    #--------------------------------------------------------------------------------------
    # Step 2: Clean the data and match its dates to periods
    #--------------------------------------------------------------------------------------

    # Work on a DataFrame that shares the columns of `data`; columns added here are not added to `data`
    data_copy = DataFrame(data; copycols = false)

    # Check that the specified columns exist and can be processed
    data_copy, treatment_times, covariates  = validate_data(data_copy, outcome, state, time, covariates, treatment_times, date_format;
                                                warn_missing_covariates = true)

    # Treated states (when given) and treatment times as vectors of the same length
    treated_states, treatment_times = check_treated_states(treated_states, treatment_times)

    # The treated states are only needed for the event study
    data_copy.state_71X9yTx = data_copy[!, state]
    if event == true
        if ci < 0 || ci >= 1
            error("'ci' must be between 0 and 1")
        end

        # The state column and treated_states must both be numbers or both be strings
        validate_state_types(treated_states, data_copy, state)
    end

    # Drop rows with missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # Check that the treated states exist
    if event == true
        treated_states, treatment_times = validate_treated_states(treated_states, treatment_times, data_copy)
    end

    # Convert dates to Date objects and keep the data between start_date and end_date
    data_copy, treatment_times, start_date, end_date, all_times, freq = validate_and_convert_dates(data_copy, time, treatment_times, date_format,
                                                                                                   freq, start_date, end_date)

    # Build the period grid and match every date and treatment time to a period on it
    data_copy, treatment_times, treated_states, _, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier, start_date,
                                                                                                 end_date, treatment_times, treated_states)

    # The event study only uses the treated states
    if event == true
        data_copy = filter(row -> row.state_71X9yTx in treated_states, data_copy)
        if nrow(data_copy) == 0
            error("No valid 'treated_states' were found in the data.")
        end
    end

    # States as strings from here on (as are treated_states)
    data_copy, treated_states = validate_string_treated_states(data_copy; treated_states = treated_states)

    # Keep the period of each observation as a Date (time_dmG5fpM) and as a string (time_71X9yTx)
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)

    # Convert factor covariates into numeric dummy variable columns
    data_copy, covariates_to_include = process_covariates(covariates, data_copy, ref)

    # A Float64 outcome column makes the regressions faster than a Union{Number, Missing} one
    data_copy.outcome_71X9yTx = convert(Vector{Float64}, data_copy.outcome_71X9yTx)

    #--------------------------------------------------------------------------------------
    # Step 3: First stage for each ccc option: lambda_{s,t}, the covariate-adjusted mean
    # outcome of each state and period (Karim & Webb 2025, eq. 45)
    #--------------------------------------------------------------------------------------

    master_lambda = DataFrame()
    for c in ccc
        lambda_df, _ = first_stage(data_copy, c, covariates_to_include, hc, false)
        master_lambda = [master_lambda;lambda_df]
    end

    # Number the periods from 0
    master_lambda.period = [time_to_index[date] for date in master_lambda.time] .- 1
    master_lambda = select(master_lambda, Not(:lambda_index))

    #--------------------------------------------------------------------------------------
    # Step 4: Arrange the lambdas for an event study plot or a parallel trends plot
    #--------------------------------------------------------------------------------------

    if event == true
        plot_data = event_study_data(master_lambda, data_copy, treated_states, treatment_times, time_to_index, weights, ci, hc)
    else
        if isnothing(date_format) || date_format == "yyyy"
            date_format = assume_date_format(period)
        end
        plot_data = parallel_trends_data(master_lambda, treatment_times, time_to_index, start_date, date_format)
    end
    plot_data.period_length .= string(period)
    return wrapper_check(plot_data, wrapper)
end


###########################################################################################
# Section 2: Parallel trends data                                                         #
###########################################################################################

function parallel_trends_data(master_lambda, treatment_times, time_to_index, start_date, date_format)

    # The lambdas of every state, period and ccc option, with dates as strings, followed by one
    # row per treatment period (only `treat_period` filled) so that plots can mark when
    # treatments start

    # Switch time column to string labels
    master_lambda.time = parse_date_to_string_didint.(master_lambda.time, date_format)
    master_lambda.start_date .= parse_date_to_string_didint(start_date, date_format)

    # The treatment periods, numbered like the periods. Which state is treated when does not
    # matter here.
    treatment_times = unique(treatment_times)
    treatment_periods = [time_to_index[date] for date in treatment_times] .-1

    # Rows for the treatment periods, with every other column missing
    treat_df = DataFrame(treat_period = treatment_periods)
    for col in names(master_lambda)
        if col != "treat_period"
            treat_df[!, col] .= missing
        end
    end
    max_period = maximum(master_lambda.period)
    treat_df = filter(row -> row.treat_period <= max_period, treat_df)

    master_lambda.treat_period .= missing
    return vcat(master_lambda, treat_df)
end


###########################################################################################
# Section 3: Event study data                                                             #
###########################################################################################

function event_study_data(master_lambda, data_copy, treated_states, treatment_times, time_to_index, weights, ci, hc)

    # For each ccc option and each period relative to treatment: the (weighted) mean of the
    # treated states' lambdas, with its HC standard error and confidence interval from the
    # regression of those lambdas on a constant

    # Period relative to treatment of every lambda
    state_to_treatment_time = Dict(zip(treated_states, treatment_times))
    master_lambda.treat_period = [haskey(state_to_treatment_time, row.state) ?
                                   time_to_index[state_to_treatment_time[row.state]] :
                                   missing
                                   for row in eachrow(master_lambda)] .- 1
    master_lambda.time_since_treatment = master_lambda.period - master_lambda.treat_period

    if weights == true
        # Weight of a state among the states observed at the same period relative to
        # treatment: its share of their observations
        data_copy.period = [time_to_index[date] for date in data_copy.time_dmG5fpM] .- 1
        data_copy.treat_period = [haskey(state_to_treatment_time, state) ?
                                   time_to_index[state_to_treatment_time[state]] :
                                   missing
                                   for state in data_copy.state_71X9yTx] .- 1
        data_copy.time_since_treatment = data_copy.period - data_copy.treat_period

        counts_df = combine(groupby(data_copy, [:state_71X9yTx, :time_since_treatment]), nrow => :count)
        time_totals = combine(groupby(counts_df, :time_since_treatment), :count => sum => :total_count)
        counts_df = leftjoin(counts_df, time_totals, on = :time_since_treatment)
        counts_df.weight = counts_df.count ./ counts_df.total_count
        master_lambda = leftjoin(master_lambda,
                                 select(counts_df, :state_71X9yTx => :state, :time_since_treatment, :weight),
                                 on = [:state, :time_since_treatment])
    elseif weights == false
        # Equal weights
        count_df = combine(groupby(master_lambda, [:time_since_treatment, :ccc]), nrow => :count_obs)
        master_lambda = leftjoin(master_lambda, count_df, on = [:time_since_treatment, :ccc])
        master_lambda.weight = 1 ./ master_lambda.count_obs
    end

    # Weighted least squares on a constant, as least squares on rows scaled by the square root
    # of their weights
    master_lambda.weighted_y = master_lambda.weight .* master_lambda.lambda
    master_lambda.intercept = sqrt.(master_lambda.weight)
    master_lambda.y = master_lambda.lambda .* master_lambda.intercept
    master_lambda.se = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
    master_lambda.ci_upper = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
    master_lambda.ci_lower = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
    master_lambda.fitted_y = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
    master_lambda.ngroup = Vector{Union{Int, Missing}}(missing, nrow(master_lambda))

    # One regression per ccc option and period relative to treatment; a standard error needs
    # at least two states
    unique_groups = unique(master_lambda[!, [:ccc, :time_since_treatment]])
    for row in eachrow(unique_groups)
        mask = (master_lambda.ccc .== row.ccc) .& (master_lambda.time_since_treatment .== row.time_since_treatment)
        group_data = master_lambda[mask, :]
        ngroup = nrow(group_data)
        master_lambda[mask, :ngroup] .= ngroup

        if ngroup >= 2

            X = reshape(group_data.intercept, :, 1)
            Y = convert(Vector{Float64}, group_data.y)
            beta_hat = safe_solve(X, Y)

            if !ismissing(beta_hat)
                resid = Y - X * beta_hat
                beta_hat_cov = compute_hc_covariance(X, resid, hc)
                beta_hat_var = diag(beta_hat_cov)
                beta_hat_se = sqrt(beta_hat_var[1])
                dof = length(Y) - 1
                t_crit = quantile(TDist(dof), 1 - ((1 - ci) / 2))
                beta_hat = first(beta_hat)
                master_lambda[mask, :fitted_y] .= beta_hat
                master_lambda[mask, :se] .= beta_hat_se
                master_lambda[mask, :ci_lower] .= beta_hat - t_crit * beta_hat_se
                master_lambda[mask, :ci_upper] .= beta_hat + t_crit * beta_hat_se
            else
                master_lambda[mask, :fitted_y] .= missing
                master_lambda[mask, :se] .= missing
                master_lambda[mask, :ci_lower] .= missing
                master_lambda[mask, :ci_upper] .= missing
            end

        end
    end

    # One row per ccc option and period relative to treatment: the weighted mean, with the
    # standard error and confidence interval of its group, sorted by ccc option and period
    event_data = combine(groupby(master_lambda, [:ccc, :time_since_treatment]),
                         :weighted_y => sum => :y,
                         :se => first => :se,
                         :ci_lower => first => :ci_lower,
                         :ci_upper => first => :ci_upper,
                         :ngroup => first => :ngroup)
    return sort!(event_data, [:ccc, :time_since_treatment])
end
