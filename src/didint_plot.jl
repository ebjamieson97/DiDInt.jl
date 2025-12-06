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
                covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                ccc::Union{AbstractString, Vector{<:AbstractString}} = "all",
                event::Bool = false,
                weights::Bool = true,
                ci::Number = 0.95,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                hc::Union{AbstractString, Number} = "hc3"
               )

The `didint_plot()` function produces a dataset in a long format that can easily be used for plotting parallel
trends or event study plots.

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
- `covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol} = nothing` 
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
    as weighted averages of state-level means for each period relative to treatment; if `false`, uses simple unweighted averages.
- `treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing`
    A vector of strings (or a single string) noting the treated state(s).
- `ci::Number = 0.95`
    Define the size of confidence bands for the event study plot.
- `hc::Union{AbstractString, Number} = "hc3"`
    Specify which heteroskedasticity-consistent covariance matrix estimator (HCCME) should be used.
    Options are `0`, `1`, `2`, `3`, and `4` (or `"hc0"`, `"hc1"`, `"hc2"`, `"hc3"`, `"hc4"`).

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

# Returns
A DataFrame of means and means residualized by the specified covariates for each of the specified common causal covariates (CCC) violations by period for each state,
or, a DataFrame of means of the treated states by periods before/after treatment (again, residualized by the specified covariates and for each of the specified
CCC violations).

# Citations
- Karim & Webb (2025). "Good Controls Gone Bad: Difference-in-Differences with Covariates". https://arxiv.org/abs/2412.14447
- MacKinnon & Webb (2020). "Randomization inference for difference-in-differences with few treated clusters". https://doi.org/10.1016/j.jeconom.2020.04.024 

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
         covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
         ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
         ccc::Union{AbstractString, Vector{<:AbstractString}} = "all",
         event::Bool = false,
         weights::Bool = true,
         ci::Number = 0.95,
         freq::Union{AbstractString, Nothing} = nothing,
         freq_multiplier::Number = 1,
         start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
         end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
         hc::Union{AbstractString, Number} = "hc3",
         wrapper::Union{AbstractString, Nothing} = nothing)

    # Check hc args
    if event
        hc = hc_checks(hc)
    end

    # Do initial wrapper check
    wrapper = init_wrapper_check(wrapper)

    # Check typeof(data)
    data_type_check(data)

   # Do check for ccc options
   if !(ccc isa AbstractVector)
       ccc = [ccc]
   end
   ccc = replace.(lowercase.(ccc), r"\s" => "")
   if ccc == ["all"]
       ccc = ["hom", "time", "state", "int", "add", "none"]
   else 
      ccc_options = ["hom", "time", "state", "int", "add", "none"]
      for c in ccc
          if !(c in ccc_options)
              error("$c is not a valid 'ccc' option. Try any combination of", join(ccc_options, ", "))
          end 
      end 
   end 

    # Turn outcome, state, time, and gvar to strings if they were inputted as symbols
    outcome, state, time, gvar = symbol_to_string(outcome, state, time, gvar)
    
    # Determine if the gvar method or the treatment_times & treated_states method is being used
    treated_states, treatment_times = gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                                              event = event, check = "didint_plot")

    # Do checks for start_date and end_date
    start_date = start_end_date_checks(start_date, "start_date", date_format)
    end_date = start_end_date_checks(end_date, "end_date", date_format)             

    # Create data copy
    data_copy = DataFrame(data; copycols = false)
    
    # Do validation checks to see if specified columns exist, and that they are able to be processed correctly
    data_copy, treatment_times, covariates  = validate_data(data_copy, outcome, state, time, covariates, treatment_times, date_format;
                                                warn_missing_covariates = true) 

    # Do checks related to treated_states
    data_copy.state_71X9yTx = data_copy[!, state]
    if event == true

        if ci < 0 || ci >= 1
            error("'ci' must be between 0 and 1")
        end 

        # Make sure treated_states is a vector
        treated_states = !(treated_states isa AbstractVector) ? [treated_states] : treated_states

        # Ensure the state column is a string or number and that the nonmissingtype(treated_states) == nonmissingtype(state column)
        validate_state_types(treated_states, data_copy, state)
    end

    # Check missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # Check for all the specified treated_states if doing event plot
    if event == true
        treated_states, treatment_times = validate_treated_states(treated_states, treatment_times, data_copy)
    end

    # Validate and convert dates to Date objects, filter data by date range
    data_copy, treatment_times, start_date, end_date, all_times, freq = validate_and_convert_dates(data_copy, time, treatment_times, date_format,
                                                                                                   freq, start_date, end_date) 

    # Do date matching procedure
    data_copy, treatment_times, treated_states, _, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier, start_date,
                                                                                                 end_date, treatment_times, treated_states)
    if event == true
        # Make sure the remaining treated_states are the only ones left in the df
        data_copy = filter(row -> row.state_71X9yTx in treated_states, data_copy)
        if nrow(data_copy) == 0 
            error("No valid 'treated_states' were found in the data.")
        end
    end

    # Ensure state column is a string
    data_copy, treated_states = validate_string_treated_states(data_copy; treated_states = treated_states, event = event)

    # Once any date matching procedures are done, convert time_71X9yTx back to a string for processing in `categorical()`
    # Also keep a column vector copy as a date
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)

    # Create dummies for each time and state interaction 
    data_copy.state_time = categorical(data_copy.state_71X9yTx .* "0IQR7q6Wei7Ejp4e" .* data_copy.time_71X9yTx)

    # Convert factor covariates into multiple numeric dummy variable columns
    data_copy, covariates_to_include = process_covariates(covariates, data_copy, ref)

    # Force outcome to float64 to speed up regression (runs faster if <:Number rather than <:Union{Number, Missing})
    data_copy.outcome_71X9yTx = convert(Vector{Float64}, data_copy.outcome_71X9yTx)

    # And now on to the actual computations
    master_lambda = DataFrame()
    for c in ccc

        # Construct formula, depending on DID-INT variation
        formula = construct_formula(c, covariates_to_include; forplot = true)

        # Run the fixed effects model and get back the dataframe of means (or means residualized by covariates) for each period at each state
        lambda_df, _, _ = run_fixed_effects_model(data_copy, formula, c, covariates)

        # Append master data
        master_lambda = [master_lambda;lambda_df]
    end

    # Add periods
    
    master_lambda.period = [time_to_index[date] for date in master_lambda.time] .- 1
    
    if event == true
        # If we are doing an event plot then we need to properly match the treatment times to states
        # Create a dictionary mapping each treated state to its treatment time
        state_to_treatment_time = Dict(zip(treated_states, treatment_times))
        
        # For each row, get the treatment period for that state
        master_lambda.treat_period = [haskey(state_to_treatment_time, row.state) ? 
                                       time_to_index[state_to_treatment_time[row.state]] : 
                                       missing 
                                       for row in eachrow(master_lambda)] .- 1
        master_lambda.time_since_treatment = master_lambda.period - master_lambda.treat_period

        if weights == true
            # Add weights
            # Add time_since_treatment to data_copy
            data_copy.period = [time_to_index[date] for date in data_copy.time_dmG5fpM] .- 1
            data_copy.treat_period = [haskey(state_to_treatment_time, state) ? 
                                       time_to_index[state_to_treatment_time[state]] : 
                                       missing 
                                       for state in data_copy.state_71X9yTx] .- 1
            data_copy.time_since_treatment = data_copy.period - data_copy.treat_period
        
        
            # Calculate counts from data_copy for each state and time_since_treatment combination
            counts_df = combine(groupby(data_copy, [:state_71X9yTx, :time_since_treatment]), nrow => :count)
            # Calculate total counts for each time_since_treatment
            time_totals = combine(groupby(counts_df, :time_since_treatment), :count => sum => :total_count)
            # Join to get totals, then calculate weights
            counts_df = leftjoin(counts_df, time_totals, on = :time_since_treatment)
            counts_df.weight = counts_df.count ./ counts_df.total_count
            # Now join the weights to master_lambda
            master_lambda = leftjoin(master_lambda, 
                                     select(counts_df, :state_71X9yTx => :state, :time_since_treatment, :weight),
                                     on = [:state, :time_since_treatment])
        elseif weights == false
            # When weights == false, just count observations per time_since_treatment
            count_df = combine(groupby(master_lambda, [:time_since_treatment, :ccc]), nrow => :count_obs)
            master_lambda = leftjoin(master_lambda, count_df, on = [:time_since_treatment, :ccc])
            master_lambda.weight = 1 ./ master_lambda.count_obs
        end

        # Prep data for regressions 
        master_lambda.weighted_y = master_lambda.weight .* master_lambda.lambda
        master_lambda.intercept = sqrt.(master_lambda.weight)
        master_lambda.y = master_lambda.lambda .* master_lambda.intercept
        master_lambda.se = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
        master_lambda.ci_upper = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
        master_lambda.ci_lower = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
        master_lambda.fitted_y = Vector{Union{Float64, Missing}}(missing, nrow(master_lambda))
        master_lambda.ngroup = Vector{Union{Int, Missing}}(missing, nrow(master_lambda))

        # Get unique combinations of ccc and time_since_treatment
        unique_groups = unique(master_lambda[!, [:ccc, :time_since_treatment]])
            
        # Loop through each ccc group and time_since_treatment
        for row in eachrow(unique_groups)
            ccc_val = row.ccc
            et = row.time_since_treatment
            mask = (master_lambda.ccc .== ccc_val) .& (master_lambda.time_since_treatment .== et)
            
            # Filter to this specific group
            group_data = filter(r -> r.ccc == ccc_val && r.time_since_treatment == et, master_lambda)
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
    
        # Collapse: sum weighted_y, take first se/ci for each ccc and time_since_treatment
        event_plot_data = combine(groupby(master_lambda, [:ccc, :time_since_treatment]),
                                  :weighted_y => sum => :y,
                                  :se => first => :se,
                                  :ci_lower => first => :ci_lower,
                                  :ci_upper => first => :ci_upper,
                                  :ngroup => first => :ngroup)

        event_plot_data.period_length .= string(period)

        event_plot_data = wrapper_check(event_plot_data, wrapper)
        return event_plot_data 

    elseif event == false

        # Switch time column to string labels
        if isnothing(date_format) || date_format == "yyyy"
            date_format = assume_date_format(period)
        end
        master_lambda.time = parse_date_to_string_didint.(master_lambda.time, date_format)
        master_lambda.start_date .= parse_date_to_string_didint(start_date, date_format)

        # If we're not doing an event plot, then it doesn't actually matter that we
        # assign the treat_periods to the right states, so long as we do add them
        treatment_times = unique(treatment_times)
        treatment_periods = [time_to_index[date] for date in treatment_times] .-1

        # Safely add the treatment_times - not dependent on the current size of master_lambda
        treat_df = DataFrame(treat_period = treatment_periods)
        for col in names(master_lambda)
            if col != "treat_period"
                treat_df[!, col] .= missing
            end
        end
        max_period = maximum(master_lambda.period)
        treat_df = filter(row -> row.treat_period <= max_period, treat_df)
        
        master_lambda.treat_period .= missing
        master_lambda = vcat(master_lambda, treat_df)
        master_lambda.period_length .= string(period)
        
        master_lambda = wrapper_check(master_lambda, wrapper)
        return master_lambda
    end

end