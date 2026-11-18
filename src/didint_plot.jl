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
                plot_type::Union{AbstractString, Vector{<:AbstractString}} = "all",
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
[`DiDIntPlot`](@ref), which holds one table per kind of plot and is drawn with `plot` once a Makie backend such as
CairoMakie is loaded.

# Details
`didint_plot()` prepares the data for five kinds of plots (the `plot_type` argument), all from the first-stage
estimates `λ(s,t)` of Karim & Webb (2025, eq. 45), the mean of the outcome for each state and period, adjusted for
covariates as each `ccc` option specifies and estimated on all states as in `didint()`:

- `"conditional"`: conditional parallel trends, the `λ(s,t)` of every state and period under every `ccc` option
  except `"int"`.
- `"conditional2"`: for each treated state, the change in `λ` of that state and of its control states since the
  period before the treated state's treatment (the pre-treatment period). These are the differences from which
  `didint()` computes its estimates.
- `"integrated"`: integrated parallel trends (`ccc = "int"`), the diagnostic of Karim & Webb (2025, eqs. 36-38). The
  outcome is regressed on the covariates separately in every state-period cell. For each treated state, the change
  from one period to the next in the predicted outcome is computed for that state and for its control states, all
  evaluated at the treated state's covariate means in the later period.
- `"integrated2"`: the same changes measured from the pre-treatment period of the treated state, which are the
  integrated differences from which `didint()` computes its estimates under `ccc = "int"`.
- `"event"`: event study, the `λ(s,t)` of the treated states averaged by period relative to treatment, for each
  `ccc` option.

In `"conditional2"` and `"integrated2"`, the change of a treated state minus the average change of its control
states is its ATT in a treated period, and a placebo estimate in a period before treatment.

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
    Every kind of plot data except `"conditional"` needs the treated states. For `"conditional"` they are optional;
    with them, the plot tells the treated states from the control states.

## Plot Data
- `plot_type::Union{AbstractString, Vector{<:AbstractString}} = "all"`
    The kinds of plot data to prepare: `"all"`, or any combination of `"event"`, `"conditional"`, `"conditional2"`,
    `"integrated"` and `"integrated2"`. With `"all"`, the kinds that the other arguments do not allow are left out:
    without the treated states only `"conditional"` is prepared (with a warning).
- `ccc::Union{AbstractString, Vector{<:AbstractString}} = "all"`
    Specify which versions of DID-INT should be used.
    Options are either `"all"`, or any combination of: `"none"`, `"hom"`, `"time"`, `"state"`, "`add`", and `"int"`.
    `"none"` uses the unadjusted means of the outcome. `"int"` is the option of the integrated parallel trends data
    (`"integrated"` and `"integrated2"`), and every other option gives conditional parallel trends data
    (`"conditional"` and `"conditional2"`). The event study uses every option given.
- `covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol, Nothing} = nothing`
    A vector of covariates entered as either strings or symbols (or a single covariate string or symbol),
    or, `nothing` (default).
- `ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing`
    A dictionary specifying which category in a categorical variable should be used
    as the reference (baseline) category.

## Control States
- `notyet::Bool = false`
    The control states of a treated state are the never-treated states and, with `notyet = true`, also the states
    of later cohorts while they are not yet treated, as in `didint()`. This applies to `"conditional2"`,
    `"integrated"` and `"integrated2"`, and to the plots that compare a treated state with its control states.

## Event Study
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
A [`DiDIntPlot`](@ref), which is drawn with `plot`. It has one table (a DataFrame) for each kind of plot data that
was prepared, in the field of the same name; the others are empty. In every table `time` is the date of a period
and `period` its number, counted from 0.

- `conditional`: one row per state, period and `ccc` option, with the state's mean outcome in that period,
  residualized by the covariates as the `ccc` option specifies (`lambda`), and the number of observations of the
  state in that period (`n`).
- `conditional2`: one row per `ccc` option, treated state (`treated_state`), state (`state`: that treated state
  itself or one of its control states) and period other than the pre-treatment period of `treated_state`
  (`base_period`), with the change in the state's `lambda` since that period (`delta`), the number of observations
  of the state in the two periods (`n`), and the period in which `treated_state` is first treated (`treat_period`).
- `integrated`: one row per treated state (`treated_state_covars`), state (`state`: that treated state itself or
  one of its control states) and period from the second one on, with the change in the state's predicted outcome
  since the previous period, evaluated at the covariate means of `treated_state_covars` (`delta`), the number of
  observations of the state in the two periods (`n`), and `treat_period`.
- `integrated2`: as `integrated`, with the change measured from the pre-treatment period of `treated_state_covars`
  (`base_period`), for every other period.
- `event`: one row per `ccc` option and period relative to treatment (`time_since_treatment`), with the mean of the
  treated states' `lambda` values (`y`), its standard error (`se`), the bounds of its confidence interval
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
         plot_type::Union{AbstractString, Vector{<:AbstractString}} = "all",
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

    # The Stata and R wrappers receive the DiDIntPlot like every other caller; `wrapper` is only checked
    wrapper = init_wrapper_check(wrapper)
    data_type_check(data)
    date_format = check_date_format(date_format)
    freq_multiplier = check_positive_integer(freq_multiplier, "freq_multiplier")

    # The ccc options
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

    # The kinds of plot data to prepare. With "all", the kinds that the other arguments do not
    # allow are left out; a kind that is asked for by name and cannot be prepared is an error.
    plot_types, all_types = requested_plot_types(plot_type)
    conditional_ccc = filter(!=("int"), ccc)
    for type in copy(plot_types)
        if type in ("conditional", "conditional2") && isempty(conditional_ccc)
            all_types || error("plot_type = \"$type\" needs a 'ccc' option other than \"int\".")
            filter!(!=(type), plot_types)
        elseif type in ("integrated", "integrated2") && !("int" in ccc)
            all_types || error("plot_type = \"$type\" needs \"int\" among the 'ccc' options.")
            filter!(!=(type), plot_types)
        end
    end

    # Every kind except "conditional" needs to know which states are treated
    if isnothing(gvar) && !isnothing(treatment_times) && isnothing(treated_states)
        if all_types && "conditional" in plot_types
            @warn "Only the conditional parallel trends data are prepared: the other plot data need the treated states. " *
                  "Specify 'treated_states' (or 'gvar') to include them."
            plot_types = ["conditional"]
        elseif plot_types != ["conditional"]
            error("The plot data other than those of plot_type = \"conditional\" need the treated states: specify 'treated_states' " *
                  "(or 'gvar').")
        end
    end

    wants(type) = type in plot_types

    # The HC covariance estimator and the confidence level are those of the event study
    hc = hc_checks(hc)
    if wants("event") && (ci < 0 || ci >= 1)
        error("'ci' must be between 0 and 1")
    end

    # Column names as strings, whether they were given as strings or symbols
    outcome, state, time, gvar = symbol_to_string(outcome, state, time, gvar)
    covariates = isnothing(covariates) ? covariates : string.(covariates)

    # Treated states and their treatment times, from gvar or from treated_states and treatment_times
    needs_treated_states = plot_types != ["conditional"]
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
    # not a control state either, as it is treated throughout.
    always_treated = has_treated_states ? string.(treated_states[treatment_times .<= start_date]) : String[]

    # Build the period grid and match every date and treatment time to a period on it
    data_copy, treatment_times, treated_states, period_dates, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier, start_date,
                                                                                                             end_date, treatment_times, treated_states)
    if needs_treated_states && isempty(treated_states)
        error("None of the treated states has a treatment time within the data used (after its first date and no later than its last date).")
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
    # outcome of each state and period (Karim & Webb 2025, eq. 45), estimated on all states
    # as in didint()
    #--------------------------------------------------------------------------------------

    # The ccc options that the requested plot data use: the event study uses all of them, the
    # conditional data every option but "int", the integrated data "int"
    needs_conditional = wants("conditional") || wants("conditional2")
    needs_integrated = wants("integrated") || wants("integrated2")
    needed_ccc = [c for c in ccc if wants("event") || (c == "int" ? needs_integrated : needs_conditional)]

    # The lambdas of every ccc option, one below the other. The cell-by-cell regressions of
    # "int" are also kept apart, with their covariate slopes and covariate means, for the
    # integrated data.
    master_lambda = DataFrame()
    int_lambda, int_fit = nothing, nothing
    for c in needed_ccc
        lambda_df, _, covariate_fit = first_stage(data_copy, c, covariates_to_include, hc, false)
        # Under "add" the plotted lambdas can depend on the least-squares solution used
        c == "add" && warn_additive_identification(covariate_fit)
        if c == "int"
            int_lambda, int_fit = lambda_df, covariate_fit
        end
        master_lambda = [master_lambda;lambda_df]
    end

    #--------------------------------------------------------------------------------------
    # Step 4: Arrange the first-stage estimates into the table of each kind of plot data
    #--------------------------------------------------------------------------------------

    if isnothing(date_format) || date_format == "yyyy"
        date_format = assume_date_format(period)
    end
    date_label(date) = parse_date_to_string_didint(date, date_format)

    # Number of observations of every state in every period
    cell_counts = combine(groupby(data_copy, [:state_71X9yTx, :time_dmG5fpM]), nrow => :n)
    n_obs = Dict(zip(zip(String.(cell_counts.state_71X9yTx), cell_counts.time_dmG5fpM), cell_counts.n))

    # The control states of a treated state in a period
    controls_of = has_treated_states ? control_states_function(unique(String.(data_copy.state_71X9yTx)), treated_states, treatment_times,
                                                               always_treated, notyet) : nothing

    conditional_lambda = master_lambda[master_lambda.ccc .!= "int", :]
    conditional = wants("conditional") ? conditional_trends_data(conditional_lambda, n_obs, time_to_index, date_label) : DataFrame()
    conditional2 = wants("conditional2") ? conditional_long_differences(conditional_lambda, n_obs, treated_states, treatment_times,
                                                                        controls_of, time_to_index, date_label) : DataFrame()
    integrated = wants("integrated") ? integrated_trends_data(int_lambda, int_fit, n_obs, treated_states, treatment_times,
                                                              controls_of, time_to_index, date_label) : DataFrame()
    integrated2 = wants("integrated2") ? integrated_long_differences(int_lambda, int_fit, n_obs, treated_states, treatment_times,
                                                                     controls_of, time_to_index, date_label) : DataFrame()
    event_lambdas, event = DataFrame(), DataFrame()
    if wants("event")
        event_lambdas = event_study_lambdas(master_lambda, n_obs, treated_states, treatment_times, time_to_index)
        event = event_study_data(event_lambdas, weights, ci, hc)
    end

    #--------------------------------------------------------------------------------------
    # Step 5: Return the tables as a DiDIntPlot, with what its plots need to know
    #--------------------------------------------------------------------------------------

    # Treatment periods, numbered from 0 like the `period` columns: those of the treated states
    # when the states are known, otherwise the distinct ones
    treatment_periods = [time_to_index[date] for date in treatment_times] .- 1
    if !has_treated_states
        keep = unique(i -> treatment_periods[i], sortperm(treatment_periods))
        treatment_times, treatment_periods = treatment_times[keep], treatment_periods[keep]
    end
    specs = (outcome = outcome, plot_types = [type for type in PLOT_TYPES if wants(type)], ccc = needed_ccc,
             covariates = isnothing(covariates) ? String[] : String.(covariates),
             treated_states = has_treated_states ? String.(treated_states) : String[],
             treatment_times = String[date_label(date) for date in treatment_times],
             treatment_periods = treatment_periods, always_treated = String.(always_treated),
             period_labels = String[date_label(date) for date in period_dates],
             notyet = notyet, weights = weights, ci = Float64(ci), hc = hc,
             period = string(period), start_date = date_label(start_date), end_date = date_label(end_date))
    return DiDIntPlot(conditional, conditional2, integrated, integrated2, event, event_lambdas, specs)
end

function requested_plot_types(plot_type)

    # The kinds of plot data asked for with the `plot_type` argument, and whether that was "all"
    types = plot_type isa AbstractVector ? plot_type : [plot_type]
    types = unique(replace.(lowercase.(types), r"\s" => ""))
    types == ["all"] && return copy(PLOT_TYPES), true
    for type in types
        if !(type in PLOT_TYPES)
            error("\"$type\" is not a valid 'plot_type'. Try \"all\" or any combination of: ", join(PLOT_TYPES, ", "))
        end
    end
    return [type for type in PLOT_TYPES if type in types], false
end


###########################################################################################
# Section 2: Parallel trends data (conditional and integrated)                            #
###########################################################################################

function control_states_function(states, treated_states, treatment_times, always_treated, notyet)

    # The control states of a treated state s in a period, as a function of s and the date t of
    # the period: the never-treated states and, with `notyet`, the states of later cohorts that
    # are not yet treated in t. From the treatment of s on these are the controls of didint();
    # before it, they are the states that are untreated both in t and in the periods around the
    # treatment of s. A state treated throughout the data (`always_treated`) is never a control.
    treated_states = String.(treated_states)
    treatment_time = Dict(zip(treated_states, treatment_times))
    never_treated = sort(String.(setdiff(states, treated_states, always_treated)), by = custom_sort_order)
    later = sort(treated_states, by = custom_sort_order)
    return function (s, t)
        notyet || return never_treated
        g = treatment_time[s]
        return vcat(never_treated, [j for j in later if treatment_time[j] > max(t, g)])
    end
end

function conditional_trends_data(lambda_df, n_obs, time_to_index, date_label)

    # The conditional parallel trends data: the lambda of every state and period under every
    # ccc option other than "int" (Karim & Webb 2025, eq. 45), with the number of observations
    # of the state in that period and the period number (from 0)
    states = String.(lambda_df.state)
    return DataFrame(state = states, time = String[date_label(t) for t in lambda_df.time],
                     period = [time_to_index[t] for t in lambda_df.time] .- 1, ccc = String.(lambda_df.ccc),
                     lambda = Vector{Float64}(lambda_df.lambda), n = [n_obs[(s, t)] for (s, t) in zip(states, lambda_df.time)])
end

function conditional_long_differences(lambda_df, n_obs, treated_states, treatment_times, controls_of, time_to_index, date_label)

    # The long differences of didint() under the ccc options other than "int" (the difference of
    # a state's lambda between a period and the period before treatment), for every period
    # instead of the treated periods only. For every ccc option and treated state s, with r the
    # period before its treatment, and every other period t:
    #   delta_{j,t} = lambda_{j,t} - lambda_{j,r}
    # for j = s and for each control state j of s in t. For t at or after the treatment of s
    # these are the differences behind ATT(s,t); for earlier t they are the same comparison
    # where no effect is expected.
    # Returns one row per (ccc option, treated state, state, period) for which the state and
    # the treated state are observed in both periods, with n = n_{j,t} + n_{j,r} (the weight of
    # a control's long difference in didint()).
    date_of_period = Dict(index => date for (date, index) in time_to_index)
    treated_states = String.(treated_states)
    treatment_time = Dict(zip(treated_states, treatment_times))

    rows = NamedTuple{(:state, :time, :period, :ccc, :delta, :n, :treated_state, :base_period, :treat_period),
                      Tuple{String, String, Int, String, Float64, Int, String, Int, Int}}[]
    for option_rows in groupby(lambda_df, :ccc; sort = false)
        option = String(first(option_rows.ccc))
        lambda = Dict(zip(zip(String.(option_rows.state), option_rows.time), option_rows.lambda))
        dates = sort(unique(option_rows.time))
        for s in sort(treated_states, by = custom_sort_order)
            g = treatment_time[s]
            base = get(date_of_period, time_to_index[g] - 1, nothing)
            (isnothing(base) || !haskey(lambda, (s, base))) && continue
            for t in dates
                (t == base || !haskey(lambda, (s, t))) && continue
                for j in vcat(s, controls_of(s, t))
                    (haskey(lambda, (j, t)) && haskey(lambda, (j, base))) || continue
                    push!(rows, (state = j, time = date_label(t), period = time_to_index[t] - 1, ccc = option,
                                 delta = lambda[(j, t)] - lambda[(j, base)], n = n_obs[(j, t)] + n_obs[(j, base)],
                                 treated_state = s, base_period = time_to_index[base] - 1, treat_period = time_to_index[g] - 1))
                end
            end
        end
    end
    return DataFrame(rows)
end

function integrated_trends_data(lambda_df, fit::CovariateFit, n_obs, treated_states, treatment_times, controls_of, time_to_index,
                                date_label)

    # The data of the integrated parallel trends plot of Karim & Webb (2025, eqs. 36-38).
    # `lambda_df` and `fit` hold the intercept lambda_{j,t}, the covariate slopes gamma_{j,t} and
    # the covariate means xbar_{j,t} from the regression of the outcome on the covariates in
    # every state-period cell (eq. 36; the first stage under ccc = "int").
    # For every treated state s and every period t that has a preceding period t-1, the change
    # in the predicted outcome of state j, evaluated at the covariate means of s in period t, is
    #   delta_{j -> s,t} = (lambda_{j,t} + gamma_{j,t}' xbar_{s,t}) - (lambda_{j,t-1} + gamma_{j,t-1}' xbar_{s,t})
    # for j = s (eq. 37) and for each control state j of s in t (eq. 38). These are the
    # integrated long differences of didint() under ccc = "int" (eqs. 46-47), taken from one
    # period to the next instead of from the period before treatment.
    # Returns one row per (treated state, state, period) for which both cells of the state and
    # the cell (s,t) of the treated state are observed, with n = n_{j,t} + n_{j,t-1}, the
    # observations behind the change.
    lambda = Dict(zip(zip(String.(lambda_df.state), lambda_df.time), lambda_df.lambda))
    predicted(state, date, x) = lambda[(state, date)] + dot(fit.slopes[(state, date)], x)
    date_of_period = Dict(index => date for (date, index) in time_to_index)
    treated_states = String.(treated_states)
    treatment_time = Dict(zip(treated_states, treatment_times))
    dates = sort(unique(lambda_df.time))

    rows = NamedTuple{(:state, :time, :period, :delta, :n, :treated_state_covars, :treat_period),
                      Tuple{String, String, Int, Float64, Int, String, Int}}[]
    for s in sort(treated_states, by = custom_sort_order), t in dates
        previous = get(date_of_period, time_to_index[t] - 1, nothing)
        (isnothing(previous) || !haskey(fit.xbar, (s, t))) && continue
        x = fit.xbar[(s, t)]
        for j in vcat(s, controls_of(s, t))
            (haskey(lambda, (j, t)) && haskey(lambda, (j, previous))) || continue
            push!(rows, (state = j, time = date_label(t), period = time_to_index[t] - 1,
                         delta = predicted(j, t, x) - predicted(j, previous, x), n = n_obs[(j, t)] + n_obs[(j, previous)],
                         treated_state_covars = s, treat_period = time_to_index[treatment_time[s]] - 1))
        end
    end
    return DataFrame(rows)
end

function integrated_long_differences(lambda_df, fit::CovariateFit, n_obs, treated_states, treatment_times, controls_of,
                                     time_to_index, date_label)

    # The integrated long differences of didint() under ccc = "int" (Karim & Webb 2025, eqs.
    # 46-47), for every period instead of the treated periods only. For every treated state s,
    # with r the period before its treatment, and every other period t:
    #   delta_{j -> s,t} = (lambda_{j,t} + gamma_{j,t}' xbar_{s,t}) - (lambda_{j,r} + gamma_{j,r}' xbar_{s,t})
    # for j = s and for each control state j of s in t. For t at or after the treatment of s
    # these are the long differences behind ATT(s,t); for earlier t they are the same
    # comparison where no effect is expected.
    # Returns one row per (treated state, state, period) for which the state is observed in
    # both periods and the cell (s,t) of the treated state is observed, with
    # n = n_{j,t} + n_{j,r} (the weight of a control's long difference in didint()).
    lambda = Dict(zip(zip(String.(lambda_df.state), lambda_df.time), lambda_df.lambda))
    predicted(state, date, x) = lambda[(state, date)] + dot(fit.slopes[(state, date)], x)
    date_of_period = Dict(index => date for (date, index) in time_to_index)
    treated_states = String.(treated_states)
    treatment_time = Dict(zip(treated_states, treatment_times))
    dates = sort(unique(lambda_df.time))

    rows = NamedTuple{(:state, :time, :period, :delta, :n, :treated_state_covars, :base_period, :treat_period),
                      Tuple{String, String, Int, Float64, Int, String, Int, Int}}[]
    for s in sort(treated_states, by = custom_sort_order)
        g = treatment_time[s]
        base = get(date_of_period, time_to_index[g] - 1, nothing)
        (isnothing(base) || !haskey(lambda, (s, base))) && continue
        for t in dates
            (t == base || !haskey(fit.xbar, (s, t))) && continue
            x = fit.xbar[(s, t)]
            for j in vcat(s, controls_of(s, t))
                (haskey(lambda, (j, t)) && haskey(lambda, (j, base))) || continue
                push!(rows, (state = j, time = date_label(t), period = time_to_index[t] - 1,
                             delta = predicted(j, t, x) - predicted(j, base, x), n = n_obs[(j, t)] + n_obs[(j, base)],
                             treated_state_covars = s, base_period = time_to_index[base] - 1, treat_period = time_to_index[g] - 1))
            end
        end
    end
    return DataFrame(rows)
end


###########################################################################################
# Section 3: Event study data                                                             #
###########################################################################################

function event_study_lambdas(master_lambda, n_obs, treated_states, treatment_times, time_to_index)

    # The lambda of every treated state and period under every ccc option, with the period
    # relative to the state's treatment (`time_since_treatment`) and the number of observations
    # of the state in that period (`n`). The event study averages these.
    treatment_period = Dict(zip(String.(treated_states), [time_to_index[date] for date in treatment_times]))
    treated_rows = master_lambda[in.(String.(master_lambda.state), Ref(keys(treatment_period))), :]
    states = String.(treated_rows.state)
    return DataFrame(ccc = String.(treated_rows.ccc), state = states,
                     time_since_treatment = [time_to_index[t] - treatment_period[s] for (s, t) in zip(states, treated_rows.time)],
                     lambda = treated_rows.lambda,
                     n = [n_obs[(s, t)] for (s, t) in zip(states, treated_rows.time)])
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
