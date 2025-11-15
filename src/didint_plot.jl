"""
    didint_plot(outcome::AbstractString,
                state::AbstractString,
                time::AbstractString,
                data::DataFrame,
                treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number},
                treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date};
                date_format::Union{AbstractString, Nothing} = nothing,
                covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
                ccc::AbstractString = "int", agg::AbstractString = "state",
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                freq::Union{AbstractString, Nothing} = nothing, freq_multiplier::Number = 1,
                autoadjust::Bool = false, nperm::Number = 1000, verbose::Bool = true,
                seed::Number = rand(1:1000000))

The `didint_plot()` function produces a collection of parallel trends plots, each one correcting
for a different violation of the common causal covariates assumption.

# Details
The arguments `treated_states` and `treated_times` should be entered such that the first element in
`treated_states` refers to the state treated at the date entered as the first element in `treated_times`,
the second element in `treated_states` refers to the state treated at the date entered as the second element
in `treated_times`, and so on. 

# Parameters
- `outcome::AbstractString` 
    Input the name of the column which identifies the outcome of interest.
- `state::AbstractString` 
    Input the name of the column which identifies the state membership of the observation.
- `time::AbstractString`
    Input the name of the column which identifies the date of the observation.
- `data::DataFrame`
    The DataFrame to be used for the analysis.
- `treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number}`
    A vector of strings (or a single string) noting the treated state(s).
- `treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date}`
    A vector (or single entry) denoting the associated treatment times of the 'treated_states'.
- `covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing}` 
    A vector of covariates entered as strings (or a single covariate string),
    or, `nothing` (default).
- `ccc::AbstractString = "int"`
    Specify which version of DID-INT should be used.
    Options are: `"hom"`, `"time"`, `"state"`, "`add`", and `"int"` (default). 
- `agg::AbstractString = "cohort"` 
    Enter the weighting method as a string.
    Options are: `"cohort"` (default), `"simple"`, `"state"`, `"sgt"`, `"none"`.
- `ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing`
    A dictionary specifying which category in a categorical variable should be used
    as the reference (baseline) category.
- `freq::Union{AbstractString, Nothing} = nothing`
    A string indicating the desired timeframe of a period for the analysis for staggered adoption scenarios.
    Options are: `"year"`, `"month"`, `"week"`, `"day"`.
- `freq_multiplier::Number = 1`
    An integer by which the 'freq' argument should be multiplied in a staggered adoption scenario, e.g. if a two-year
    period is desired, set `freq = "year"` and `freq_multiplier = 2`.
- `autoadjust::Bool = false`
    A boolean option for the length of a period to be considered in staggered adoption scenario to be automatically
    determined by DiDInt.
- `nperm::Number = 1000`
    The number of unique permutations to be considered when performing the randomization inference.
- `verbose::Bool = true`
    A boolean option for displaying progress of the randomization procedure. 
- `seed::Number = rand(1:1000000)`
    An integer to set the random seed for the randomization inference procedure.

# Returns
A DataFrame of results including the estimate of the ATT as well as standard errors and p-values.

"""
function didint_plot(
         outcome::Union{AbstractString, Symbol},
         state::Union{AbstractString, Symbol},
         time::Union{AbstractString, Symbol},
         data::DataFrame;
         gvar::Union{AbstractString, Symbol, Nothing} = nothing,
         treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
         treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
         date_format::Union{AbstractString, Nothing} = nothing,
         covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
         ccc::Union{AbstractString, Vector{<:AbstractString}} = "all",
         event::Bool = false,
         weights::Bool = true,
         ci::Number = 0.95,
         freq::Union{AbstractString, Nothing} = nothing,
         freq_multiplier::Number = 1,
         start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
         end_date::Union{AbstractString, Number, Date, Nothing} = nothing)

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
    data_copy, treatment_times  = validate_data(data_copy, outcome, state, time, covariates, treatment_times, date_format;
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
        missing_states = setdiff(treated_states, data_copy.state_71X9yTx)
        if !isempty(missing_states)
            error("The following 'treated_states' could not be found in the data $(missing_states).\nOnly found the following states $(unique(data_copy.state_71X9yTx))")
        end

    end

    # Check missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # Validate and convert dates to Date objects, filter data by date range
    data_copy, treatment_times, start_date, end_date, all_times, freq = validate_and_convert_dates(data_copy, time, treatment_times, date_format,
                                                                                                   freq, start_date, end_date) 

    # Do date matching procedure
    data_copy, treatment_times, treated_states, _, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier, start_date,
                                                                                                 end_date, treatment_times, treated_states)
    if event == true
        # Make sure the remaining treated_states are the only ones left in the df
        data_copy = filter(row -> row.state_71X9yTx in treated_states, data_copy)
    end

    # Ensure state column is a string
    data_copy.state_71X9yTx = string.(data_copy.state_71X9yTx)

    # Once any date matching procedures are done, convert time_71X9yTx back to a string for processing in `categorical()`
    # Also keep a column vector copy as a date
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)

    # Create dummies for each time and state interaction 
    data_copy.state_time = categorical(data_copy.state_71X9yTx .* "0IQR7q6Wei7Ejp4e" .* data_copy.time_71X9yTx)

    # Convert factor covariates into multiple numeric dummy variable columns
    covariates_to_include = String[]
    if !isnothing(covariates)
        for cov in covariates
            cov_type = Base.nonmissingtype(eltype(data_copy[!, cov]))
            if cov_type <: AbstractString || cov_type <: CategoricalValue
                unique_categories = unique(data_copy[!, cov])
                if length(unique_categories) >= 2 

                    if !isnothing(ref) && haskey(ref, cov)
                        refcat = ref[cov]
                        if !(refcat in unique_categories)
                            error("$refcat $cov Reference category '$refcat' not found in column '$cov'.")
                        end
                    else
                        refcat = first(unique_categories)
                    end

                    if length(unique_categories) > 2
                        for category in unique_categories
                            if category != refcat
                                newcol = Symbol("$(cov)_$(category)")
                                data_copy[!, newcol] = (data_copy[!, cov] .== category) .|> Int
                                push!(covariates_to_include, string(newcol))
                            end
                        end
                    elseif length(unique_categories) == 2
                        newcol = Symbol("$(cov)_not_$(refcat)")
                        data_copy[!, cov] = (data_copy[!, cov] .!= refcat) .|> Int
                        push!(covariates_to_include, string(newcol))
                    end
                else    
                    error("$cov Only detected one unique factor ($unique_categories) in factor variable $cov.")
                end 
            elseif cov_type <:Number
                data_copy[!, cov] = convert(Vector{Float64}, data_copy[!, cov])
                push!(covariates_to_include, cov)
            else
                error("$cov column was found to be ($cov_type) neither of type Number, AbstractString, nor CategoricalValue!")
            end
        end
    end

    # Force outcome to float64 to speed up regression (runs faster is <:Number rather than <:Union{Number, Missing})
    data_copy.outcome_71X9yTx = convert(Vector{Float64}, data_copy.outcome_71X9yTx)

    # And now on to the actual computations
    master_lambda = DataFrame()
    for c in ccc
        # Construct formula, depending on DID-INT variation
        formula_str = "@formula(outcome_71X9yTx ~ 0 + state_time"

        # Set up formula
        if c == "int"
            for c in covariates_to_include
                formula_str *= " + fe(state_71X9yTx)&fe(time_71X9yTx)&$c"
            end
        elseif c == "time"
            for c in covariates_to_include
                formula_str *= " + fe(time_71X9yTx)&$c"
            end
        elseif c == "state"
            for c in covariates_to_include
                formula_str *= " + fe(state_71X9yTx)&$c"
            end
        elseif c == "add"
            for c in covariates_to_include
                formula_str *= " + fe(time_71X9yTx)&$c + fe(state_71X9yTx)&$c"
            end
        elseif c == "hom"
            for c in covariates_to_include
                formula_str *= " + $c"
            end
        elseif c == "none"
            formula_str *= ""
        end 
        formula_str *= ")"
        formula_expr = Meta.parse(formula_str)
        formula = eval(formula_expr)

        # Call GC.gc() before running big FixedEffectsModels regression
        GC.gc()

        stage1 = reg(data_copy, formula, contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
                         save = false)

        # Recover lambdas
        state_time = split.(replace.(coefnames(stage1), "state_time: " => ""), "0IQR7q6Wei7Ejp4e")
        coefs = coef(stage1)
        if c == "hom" && !isnothing(covariates)
            keep_mask = map(x -> !(x[1] in covariates), state_time)
            state_time = state_time[keep_mask]
            coefs = coefs[keep_mask]
        end
        lambda_df = DataFrame(state = first.(state_time), time = last.(state_time))
        lambda_df.lambda = coefs

        # Handle the edge case
        expected_df = unique(data_copy[:, [:state_71X9yTx, :time_71X9yTx]])
        expected_df = rename(expected_df, :state_71X9yTx => :state, :time_71X9yTx => :time)
        lambda_df = leftjoin(expected_df, lambda_df, on = [:state, :time])
        lambda_df.lambda = coalesce.(lambda_df.lambda, 0.0)

        # Parse the treatment_times and the time column to dates for staggered adoption
        lambda_df.time = Date.(lambda_df.time)

        # Note the ccc for this loop
        lambda_df.ccc .= c

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

        # Get unique combinations of ccc and time_since_treatment
        unique_groups = unique(master_lambda[!, [:ccc, :time_since_treatment]])
            
        # Loop through each ccc group and time_since_treatment
        for row in eachrow(unique_groups)
            ccc_val = row.ccc
            et = row.time_since_treatment
            
            # Filter to this specific group
            group_data = filter(r -> r.ccc == ccc_val && r.time_since_treatment == et, master_lambda)
        
            if nrow(group_data) > 1
                result = lm(@formula(y ~ 0 + intercept), group_data)
                
                # Get coefficient
                coef_val = coef(result)[1]

                # Manual HC1 robust standard error calculation
                X = hcat(group_data.intercept) 
                residuals = group_data.y .- X * coef_val 
                n = length(residuals)
                k = 1 

                # HC1 adjustment factor (Stata default)
                adj = n / (n - k)

                # Compute meat of sandwich estimator
                bread = inv(X' * X)
                meat = sum((residuals[i]^2) * (X[i, :]' * X[i, :]) for i in 1:n)

                # HC1 robust variance-covariance matrix
                vcov_robust = adj * bread * meat * bread

                # Standard error is the square root of the diagonal
                se_val = sqrt(vcov_robust[1, 1])
                
                # Calculate critical t-value for CI
                df = dof_residual(result)
                t_crit = quantile(TDist(df), 1 - (1 - ci) / 2)
                
                # Update master_lambda for this group
                mask = (master_lambda.ccc .== ccc_val) .& (master_lambda.time_since_treatment .== et)
                master_lambda[mask, :fitted_y] .= coef_val  
                master_lambda[mask, :se] .= se_val
                master_lambda[mask, :ci_lower] .= coef_val - t_crit * se_val
                master_lambda[mask, :ci_upper] .= coef_val + t_crit * se_val
            end
        end
    
    # Collapse: sum weighted_y, take first se/ci for each ccc and time_since_treatment
    event_plot_data = combine(groupby(master_lambda, [:ccc, :time_since_treatment]),
                              :weighted_y => sum => :y,
                              :se => first => :se,
                              :ci_lower => first => :ci_lower,
                              :ci_upper => first => :ci_upper)

    event_plot_data.period_length .= string(period)

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
        
        return master_lambda
    end

end