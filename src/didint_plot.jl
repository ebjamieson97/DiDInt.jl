###########################################################################################
# didint_plot.jl: data for parallel trends and event study plots (drawn as src/plots.jl   #
# describes)                                                                              #
#                                                                                         #
#   Section 1: didint_plot(), the main function                                           #
#   Section 2: Parallel trends data (conditional and integrated)                          #
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
                notyet::Bool = false,
                weights::Bool = true,
                ci::Number = 0.95,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                hc::Union{AbstractString, Number} = "hc1"
               )

The `didint_plot()` function prepares the data for parallel trends plots and event study plots. It returns them as a
[`DiDIntPlot`](@ref): the data are in a long format in its `data` field, and `plot` draws them once a Makie backend
such as CairoMakie is loaded.

# Details
`didint_plot()` prepares the data for three kinds of plots:

- Conditional parallel trends (`event = false`, every `ccc` option except `"int"`): the first-stage estimates `λ(s,t)`
  of Karim & Webb (2025, eq. 45), the mean of the outcome for each state and period, adjusted for covariates as the
  `ccc` option specifies.
- Integrated parallel trends (`event = false`, `ccc = "int"`), the diagnostic of Karim & Webb (2025, eqs. 36-38). The
  outcome is regressed on the covariates separately in every state-period cell. For each treated state, the change
  from one period to the next in the predicted outcome is computed for that state and for its comparison states, all
  evaluated at the treated state's covariate means in the later period.
- Event study (`event = true`): the `λ(s,t)` of the treated states, averaged by period relative to treatment, for
  each `ccc` option.

`didint_plot()` is a visual aid for judging parallel trends; it does not reproduce the `didint()` estimates.

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
- `treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing`
    A vector of strings or numbers (or a single one) noting the treated state(s), matching the type of the `state` column.
    Needed for the event study and for the integrated parallel trends data. For the conditional parallel trends data
    it is optional; with it, the plots tell the treated states from the control states.

## Model Specifications
- `ccc::Union{AbstractString, Vector{<:AbstractString}} = "all"`
    Specify which versions of DID-INT should be used.
    Options are either `"all"`, or any combination of: `"none"`, `"hom"`, `"time"`, `"state"`, "`add`", and `"int"`.
    `"none"` uses the unadjusted means of the outcome. With `event = false`, `"int"` gives the integrated parallel
    trends data, which need the treated states (`treated_states` or `gvar`), and every other option the
    conditional parallel trends data.
- `covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing`
    A vector of covariates entered as either strings or symbols (or a single covariate string or symbol),
    or, `nothing` (default).
- `ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing`
    A dictionary specifying which category in a categorical variable should be used
    as the reference (baseline) category.

## Control States
- `notyet::Bool = false`
    The comparison states of a treated state are the never-treated states and, with `notyet = true`, also the
    states not yet treated in the period in question. This applies to the integrated parallel trends data
    (`ccc = "int"`) and to the plots that compare a treated state with its control states.

## Event Study Plot
- `event::Bool = false`
   Specify if data should be prepared for an event study plot as opposed to a parallel trends plot.
- `weights::Bool = true`
    Whether to use weighted means when computing event study estimates. If `true`, estimates are computed
    as weighted averages of state-level means for each period relative to treatment, weighted by the states'
    numbers of observations; if `false`, uses simple unweighted averages.
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
A [`DiDIntPlot`](@ref), which is drawn with `plot` and holds the plot data as a DataFrame in its `data` field.

With `event = false`, the `ccc` column of the DataFrame tells the two kinds of rows apart:

- Conditional parallel trends (every `ccc` option except `"int"`): one row per state, period and `ccc` option, with
  the state's mean outcome in that period, residualized by the covariates as the `ccc` option specifies (`lambda`),
  and the number of observations of the state in that period (`n`).
- Integrated parallel trends (`ccc = "int"`): one row per treated state (`treated_state_covars`), state (`state`: that
  treated state itself or one of its comparison states) and period from the second one on, with the change in the
  state's predicted outcome since the previous period, evaluated at the covariate means of `treated_state_covars`
  (`delta`), the number of observations of the state in the two periods (`n`), and the period in which
  `treated_state_covars` is first treated (`treat_period`).

Both have the period number (`period`). `lambda` is missing in the integrated rows, and `delta` and
`treated_state_covars` in the conditional rows. When there are conditional rows, extra rows with only `treat_period`
filled give the periods in which treatments start.

With `event = true`, the DataFrame has one row per `ccc` option and period relative to treatment (`time_since_treatment`):
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
         notyet::Bool = false,
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
    all_ccc = ccc == ["all"]
    if all_ccc
        ccc = ccc_options
    else
        for c in ccc
            if !(c in ccc_options)
                error("$c is not a valid 'ccc' option. Try any combination of: ", join(ccc_options, ", "))
            end
        end
    end

    # For a parallel trends plot, ccc = "int" gives the integrated parallel trends data (Karim &
    # Webb 2025, eqs. 36-38) instead of the lambdas, and these need to know which states are
    # treated. Without the treated states, "int" is an error if it was asked for by name, and is
    # left out (with a warning) if it is only there through ccc = "all".
    integrated = !event && "int" in ccc
    if integrated && isnothing(gvar) && !isnothing(treatment_times) && isnothing(treated_states)
        if all_ccc
            @warn "The integrated parallel trends data of ccc = \"int\" need the treated states: specify 'treated_states' " *
                  "(or 'gvar') to include them. \"int\" is left out."
            ccc = filter(!=("int"), ccc)
            integrated = false
        else
            error("The integrated parallel trends data of ccc = \"int\" need the treated states: specify 'treated_states' " *
                  "(or 'gvar'), or leave \"int\" out of 'ccc'.")
        end
    end

    # Column names as strings, whether they were given as strings or symbols
    outcome, state, time, gvar = symbol_to_string(outcome, state, time, gvar)
    covariates = isnothing(covariates) ? covariates : string.(covariates)

    # Treated states and their treatment times, from gvar or from treated_states and treatment_times.
    # The event study and the integrated parallel trends data need to know which states are treated;
    # the conditional parallel trends data do not, but their plots mark the treated states when
    # they are given.
    needs_treated_states = event || integrated
    treated_states, treatment_times = gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                                              event = needs_treated_states, check = "didint_plot")
    has_treated_states = !isnothing(treated_states)

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

    data_copy.state_71X9yTx = data_copy[!, state]
    if event == true && (ci < 0 || ci >= 1)
        error("'ci' must be between 0 and 1")
    end
    if has_treated_states
        # The state column and treated_states must both be numbers or both be strings
        validate_state_types(treated_states, data_copy, state)
    end

    # Drop rows with missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # Check that the treated states exist
    if has_treated_states
        treated_states, treatment_times = validate_treated_states(treated_states, treatment_times, data_copy)
    end

    # Convert dates to Date objects and keep the data between start_date and end_date
    data_copy, treatment_times, start_date, end_date, all_times, freq = validate_and_convert_dates(data_copy, time, treatment_times, date_format,
                                                                                                   freq, start_date, end_date)

    # A state treated at or before the first date is not counted as a treated state below. It is
    # not a comparison state either (in the integrated parallel trends data and in the plots), as
    # it is treated throughout.
    always_treated = has_treated_states ? string.(treated_states[treatment_times .<= start_date]) : String[]

    # Build the period grid and match every date and treatment time to a period on it
    data_copy, treatment_times, treated_states, _, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier, start_date,
                                                                                                 end_date, treatment_times, treated_states)
    if integrated && isempty(treated_states)
        error("None of the treated states has a treatment time within the data used (after its first date and no later than its last date).")
    end

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

    # The lambdas of every ccc option, one below the other. For the integrated parallel trends
    # data, the cell-by-cell regressions of "int" are kept apart, with their covariate slopes and
    # covariate means.
    master_lambda = DataFrame()
    int_lambda, int_fit = nothing, nothing
    for c in ccc
        lambda_df, _, covariate_fit = first_stage(data_copy, c, covariates_to_include, hc, false)
        # Under "add" the plotted lambdas can depend on the least-squares solution used
        c == "add" && warn_additive_identification(covariate_fit)
        if integrated && c == "int"
            int_lambda, int_fit = lambda_df, covariate_fit
        else
            master_lambda = [master_lambda;lambda_df]
        end
    end

    #--------------------------------------------------------------------------------------
    # Step 4: Arrange the first-stage estimates for an event study plot or for parallel trends
    # plots (conditional and integrated)
    #--------------------------------------------------------------------------------------

    if isnothing(date_format) || date_format == "yyyy"
        date_format = assume_date_format(period)
    end
    # Number of observations of every state in every period
    cell_counts = combine(groupby(data_copy, [:state_71X9yTx, :time_dmG5fpM]), nrow => :n)
    n_obs = Dict(zip(zip(String.(cell_counts.state_71X9yTx), cell_counts.time_dmG5fpM), cell_counts.n))

    if event == true
        event_lambdas = event_study_lambdas(master_lambda, n_obs, treated_states, treatment_times, time_to_index)
        plot_data = event_study_data(event_lambdas, weights, ci, hc)
    else
        event_lambdas = DataFrame()
        conditional = conditional_trends_data(master_lambda, n_obs, time_to_index)
        integrated_rows = integrated ? integrated_trends_data(int_lambda, int_fit, n_obs, treated_states, treatment_times,
                                                              always_treated, time_to_index, notyet) : nothing
        plot_data = parallel_trends_data(conditional, integrated_rows, treatment_times, time_to_index, start_date, date_format)
    end
    plot_data.period_length .= string(period)

    # The Stata and R wrappers receive the table itself
    if wrapper in ("stata", "r")
        return wrapper_check(plot_data, wrapper)
    end

    #--------------------------------------------------------------------------------------
    # Step 5: Return the data as a DiDIntPlot, with what its plots need to know
    #--------------------------------------------------------------------------------------

    # Treatment periods, numbered from 0 like the `period` column: those of the treated states
    # when the states are known, otherwise the distinct ones
    treatment_periods = [time_to_index[date] for date in treatment_times] .- 1
    if !has_treated_states
        keep = unique(i -> treatment_periods[i], sortperm(treatment_periods))
        treatment_times, treatment_periods = treatment_times[keep], treatment_periods[keep]
    end
    specs = (outcome = outcome, event = event, ccc = ccc, covariates = isnothing(covariates) ? String[] : String.(covariates),
             treated_states = has_treated_states ? String.(treated_states) : String[],
             treatment_times = String[parse_date_to_string_didint(date, date_format) for date in treatment_times],
             treatment_periods = treatment_periods, always_treated = String.(always_treated),
             notyet = notyet, weights = weights, ci = Float64(ci), hc = event ? hc : nothing,
             period = string(period), start_date = parse_date_to_string_didint(start_date, date_format),
             end_date = parse_date_to_string_didint(end_date, date_format))
    return DiDIntPlot(plot_data, event_lambdas, specs)
end


###########################################################################################
# Section 2: Parallel trends data (conditional and integrated)                            #
###########################################################################################

# Columns of the parallel trends data, in their order (period_length is added by didint_plot)
const PARALLEL_TRENDS_COLUMNS = [:state, :time, :lambda, :delta, :n, :treated_state_covars, :ccc, :period,
                                 :start_date, :treat_period]

function conditional_trends_data(master_lambda, n_obs, time_to_index)

    # The conditional parallel trends data: the lambda of every state and period under every
    # ccc option other than "int" (Karim & Webb 2025, eq. 45), with the number of observations
    # of the state in that period and the period number (from 0). Empty if no such option was
    # asked for.
    nrow(master_lambda) == 0 && return DataFrame(state = String[], time = Date[], lambda = Float64[], n = Int[],
                                                 ccc = String[], period = Int[])
    return DataFrame(state = String.(master_lambda.state), time = master_lambda.time, lambda = master_lambda.lambda,
                     n = [n_obs[(String(s), t)] for (s, t) in zip(master_lambda.state, master_lambda.time)],
                     ccc = String.(master_lambda.ccc), period = [time_to_index[t] for t in master_lambda.time] .- 1)
end

function integrated_trends_data(lambda_df, fit::CovariateFit, n_obs, treated_states, treatment_times, always_treated,
                                time_to_index, notyet)

    # The data of the integrated parallel trends plot of Karim & Webb (2025, eqs. 36-38).
    # `lambda_df` and `fit` hold the intercept lambda_{j,t}, the covariate slopes gamma_{j,t} and
    # the covariate means xbar_{j,t} from the regression of the outcome on the covariates in
    # every state-period cell (eq. 36; the first stage under ccc = "int").
    # For every treated state s and every period t that has a preceding period t-1, the change
    # in the predicted outcome of state j, evaluated at the covariate means of s in period t, is
    #   delta_{j -> s,t} = (lambda_{j,t} + gamma_{j,t}' xbar_{s,t}) - (lambda_{j,t-1} + gamma_{j,t-1}' xbar_{s,t})
    # for j = s (eq. 37) and for each comparison state j of s (eq. 38): the never-treated states
    # and, with `notyet`, the states not yet treated in period t. A state that is treated
    # throughout the data (`always_treated`) is not a comparison state. These are the integrated
    # long differences of didint() under ccc = "int" (eqs. 46-47), taken from one period to the
    # next instead of from the period before treatment.
    # Returns one row per (treated state, state, period) for which both cells of the state and
    # the cell (s,t) of the treated state are observed, with n = n_{j,t} + n_{j,t-1}, the
    # observations behind the change (the weight of a control's long difference in didint()).
    lambda = Dict(zip(zip(String.(lambda_df.state), lambda_df.time), lambda_df.lambda))
    predicted(state, date, x) = lambda[(state, date)] + dot(fit.slopes[(state, date)], x)
    date_of_period = Dict(index => date for (date, index) in time_to_index)
    treated_states = String.(treated_states)
    treatment_time = Dict(zip(treated_states, treatment_times))
    never_treated = sort(String.(setdiff(unique(lambda_df.state), treated_states, always_treated)), by = custom_sort_order)
    dates = sort(unique(lambda_df.time))

    rows = NamedTuple{(:state, :time, :delta, :n, :treated_state_covars, :ccc, :period, :treat_period),
                      Tuple{String, Date, Float64, Int, String, String, Int, Int}}[]
    for s in sort(treated_states, by = custom_sort_order), t in dates
        previous = get(date_of_period, time_to_index[t] - 1, nothing)
        (isnothing(previous) || !haskey(fit.xbar, (s, t))) && continue
        x = fit.xbar[(s, t)]
        comparisons = never_treated
        if notyet
            not_yet = [j for j in treated_states if j != s && treatment_time[j] > t]
            comparisons = vcat(never_treated, sort(not_yet, by = custom_sort_order))
        end
        for j in vcat(s, comparisons)
            (haskey(lambda, (j, t)) && haskey(lambda, (j, previous))) || continue
            push!(rows, (state = j, time = t, delta = predicted(j, t, x) - predicted(j, previous, x),
                         n = n_obs[(j, t)] + n_obs[(j, previous)], treated_state_covars = s, ccc = "int",
                         period = time_to_index[t] - 1, treat_period = time_to_index[treatment_time[s]] - 1))
        end
    end
    return DataFrame(rows)
end

function parallel_trends_data(conditional, integrated, treatment_times, time_to_index, start_date, date_format)

    # The parallel trends data in one table: the conditional rows (`lambda`), then the integrated
    # rows of ccc = "int" (`delta`, `treated_state_covars` and the treatment period of that
    # treated state in `treat_period`), with dates as strings. The columns that a kind of row
    # does not have are missing. When there are conditional rows, they are followed by one row
    # per treatment period (only `treat_period` filled) so that plots can mark when treatments
    # start; which state is treated when does not matter there.
    pieces = DataFrame[conditional]
    isnothing(integrated) || push!(pieces, integrated)
    plot_data = vcat(pieces...; cols = :union)
    for col in PARALLEL_TRENDS_COLUMNS
        col in (:start_date,) && continue
        string(col) in names(plot_data) || (plot_data[!, col] = Vector{Missing}(missing, nrow(plot_data)))
    end
    plot_data.time = parse_date_to_string_didint.(plot_data.time, date_format)
    plot_data.start_date .= parse_date_to_string_didint(start_date, date_format)
    plot_data = select(plot_data, PARALLEL_TRENDS_COLUMNS)

    if nrow(conditional) > 0
        # Rows for the treatment periods, with every other column missing
        treatment_periods = [time_to_index[date] for date in unique(treatment_times)] .- 1
        treatment_periods = filter(<=(maximum(conditional.period)), treatment_periods)
        treat_df = DataFrame(treat_period = treatment_periods)
        for col in PARALLEL_TRENDS_COLUMNS
            col == :treat_period || (treat_df[!, col] .= missing)
        end
        plot_data = vcat(plot_data, select(treat_df, PARALLEL_TRENDS_COLUMNS))
    end
    return plot_data
end


###########################################################################################
# Section 3: Event study data                                                             #
###########################################################################################

function event_study_lambdas(master_lambda, n_obs, treated_states, treatment_times, time_to_index)

    # The lambda of every treated state and period under every ccc option, with the period
    # relative to the state's treatment (`time_since_treatment`) and the number of observations
    # of the state in that period (`n`). The event study averages these.
    treatment_period = Dict(zip(String.(treated_states), [time_to_index[date] for date in treatment_times]))
    states = String.(master_lambda.state)
    return DataFrame(ccc = String.(master_lambda.ccc), state = states,
                     time_since_treatment = [time_to_index[t] - treatment_period[s] for (s, t) in zip(states, master_lambda.time)],
                     lambda = master_lambda.lambda,
                     n = [n_obs[(s, t)] for (s, t) in zip(states, master_lambda.time)])
end

function event_study_data(event_lambdas, weights, ci, hc)

    # For each ccc option and each period relative to treatment: the mean of the treated
    # states' lambdas (`y`), with its HC standard error and confidence interval from the
    # regression of those lambdas on a constant, and the number of states behind it (`ngroup`).
    # With `weights`, a state is weighted by its share of the observations of the states
    # observed at the same period relative to treatment; otherwise the states count equally.
    # Rows are sorted by ccc option and period.
    rows = NamedTuple{(:ccc, :time_since_treatment, :y, :se, :ci_lower, :ci_upper, :ngroup),
                      Tuple{String, Int, Float64, Union{Float64, Missing}, Union{Float64, Missing}, Union{Float64, Missing}, Int}}[]
    for group in groupby(event_lambdas, [:ccc, :time_since_treatment])
        ngroup = nrow(group)
        weight = weights ? group.n ./ sum(group.n) : fill(1 / ngroup, ngroup)
        y = sum(weight .* group.lambda)
        se, ci_lower, ci_upper = missing, missing, missing

        # Weighted least squares on a constant, as least squares on rows scaled by the square
        # root of their weights; a standard error needs at least two states
        if ngroup >= 2
            X = reshape(sqrt.(weight), :, 1)
            Y = convert(Vector{Float64}, group.lambda .* sqrt.(weight))
            beta_hat = safe_solve(X, Y)
            if !ismissing(beta_hat)
                resid = Y - X * beta_hat
                se = sqrt(diag(compute_hc_covariance(X, resid, hc))[1])
                t_crit = quantile(TDist(ngroup - 1), 1 - ((1 - ci) / 2))
                ci_lower = first(beta_hat) - t_crit * se
                ci_upper = first(beta_hat) + t_crit * se
            end
        end
        push!(rows, (ccc = first(group.ccc), time_since_treatment = first(group.time_since_treatment), y = y, se = se,
                     ci_lower = ci_lower, ci_upper = ci_upper, ngroup = ngroup))
    end
    event_data = sort!(DataFrame(rows), [:ccc, :time_since_treatment])
    for column in (:se, :ci_lower, :ci_upper)
        any(ismissing, event_data[!, column]) || disallowmissing!(event_data, column)
    end
    return event_data
end
