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
       ccc = ["hom", "time", "state", "int", "add"]
   else 
      ccc_options = ["hom", "time", "state", "int", "add"]
      for c in ccc
          if !(c in ccc_options)
              error("$c is not a valid 'ccc' option. Try any combination of", join(ccc_options, ", "))
          end 
      end 
   end 

    # Turn outcome, state, time, and gvar to strings if they were inputted as symbols
    if typeof(outcome) <: Symbol
        outcome = string(outcome)
    end 
    if typeof(state) <: Symbol
        state = string(state)
    end
    if typeof(time) <: Symbol
        time = string(time)
    end
    if typeof(gvar) <: Symbol
        gvar = string(gvar)
    end
    
    # Determine if the gvar method or the treatment_times & treated_states method is being used
    if isnothing(gvar)
        if (event == true  && (isnothing(treatment_times) || isnothing(treated_states))) ||
            (event == false && isnothing(treatment_times))
                error("If 'gvar' is not specified, then 'treatment_times'' must be specified!\nIf 'event = true' then 'treated_states' must also be specified!")
        end
    else
        if !isnothing(treatment_times) || !isnothing(treated_states)
            error("If 'gvar' is specified, do not specify the 'treatment_times' and 'treated_states'!")
        end
        mask = .!ismissing.(data[!, gvar]) .&& data[!, gvar] .!= 0 .&& data[!, gvar] .!= "0" .&& data[!, gvar] .!= "0.0"
        filtered = unique(data[mask, [state, gvar]])
        treated_states = filtered[!, state]
        treatment_times = filtered[!, gvar]
    end

    # Do checks for start_date and end_date
    if !isnothing(start_date)
        if eltype(start_date) <: Number
            if length(string(start_date)) != 4
                error("If 'start_date' is a number it must be 4 digits long.")
            else 
                start_date = Date(start_date)
            end 
        elseif start_date isa String 
            start_date = parse_string_to_date_didint(start_date, date_format)
        end 
    end
    if !isnothing(end_date)
        if eltype(end_date) <: Number
            if length(string(end_date)) != 4
                error("If 'end_date' is a number it must be 4 digits long.")
            else 
                end_date = Date(end_date)
            end 
        elseif end_date isa String
            end_date = parse_string_to_date_didint(end_date, date_format)
        end 
    end            

    # Create data copy
    data_copy = DataFrame(data; copycols = false)
    
    # Ensure the specified outcome, state, and time columns exist in the data
    missing_cols = [col for col in [outcome, state, time] if !(col in names(data_copy))]
    if !isempty(missing_cols)
        error("The following columns could not be found in the data ", join(missing_cols, ", "))
    end

    # Ensure the outcome variable is a numeric variable and create outcome column
    outcome_nonmissingtype = Base.nonmissingtype(eltype(data_copy[!, outcome]))
    if !(outcome_nonmissingtype <: Number)
        error("Column '$outcome' must be numeric, but found $(eltype(data_copy[!, outcome]))")
    end
    data_copy.outcome_71X9yTx = data_copy[!, outcome]

    # Check that the specified covariates exist
    if !isnothing(covariates)
        if isa(covariates, AbstractString)
            covariates = [covariates]
        end
        missing_cov = [col for col in covariates if !(col in names(data_copy))]
        if !isempty(missing_cov)
            error("$(join(missing_cov, ", ")) The preceding covariates could not be found in the data.")
        end
    end

    # Ensure that treatment_times, if a number, are all 4 digit entries
    if Base.nonmissingtype(eltype(treatment_times)) <: Number
        treatment_times = round.(Int, treatment_times)
        treatment_times_numeric = true
        if isnothing(date_format) 
            date_format = "yyyy"
        end 
        if lowercase(date_format) != "yyyy"
            error("If 'treatment_times' are entered as a number, the 'date_format' must be \"yyyy\".")
        end 
        if treatment_times isa AbstractVector
            unique_lengths = unique(length.(string.(treatment_times)))
        else 
            unique_lengths = [length(string(treatment_times))]
            treatment_times = [treatment_times]
        end
        if length(unique_lengths) == 1
            if unique_lengths[1] != 4
                error("If 'treatment_times' are entered as numbers, they must all be 4 digits long in 'yyyy' date_format.'treatment_times' found: $(join(treatment_times, ", "))")
            end 
        else
            error("Detected multiple unique date_formats in 'treatment_times'.")
        end
        treatment_times = Date.(treatment_times)
    else 
        treatment_times_numeric = false
    end
    nonmissing_time_type = Base.nonmissingtype(eltype(data_copy[!, time]))
    if nonmissing_time_type <: Number
        data_copy[!, time] = round.(Int, data_copy[!, time])
        time_column_numeric = true
        unique_time_lengths = unique(length.(string.(data_copy[!, time])))
        if length(unique_time_lengths) > 1 || unique_time_lengths[1] != 4
            error("The 'time' column was found to be numeric but consisting of values of ambiguous date formatting (i.e. not consistent 4 digit entries.)")
        end 
    else
        time_column_numeric = false
    end 
    if xor(time_column_numeric, treatment_times_numeric)
        error("If 'time' column is numeric or 'treatment_times' is numeric, then both must be numeric.")
    end 

    # Detect if staggered adoption or common treatment time, also check if treatment_times is a vector
    if !(treatment_times isa AbstractVector)
        treatment_times = [treatment_times]
    end

    # Do checks related to treated_states
    data_copy.state_71X9yTx = data_copy[!, state]
    if event == true

        if ci < 0 || ci >= 1
            error("'ci' must be between 0 and 1")
        end 

        # Make sure treated_states and control_states are vectors, make sure control_states exist
        if !(treated_states isa AbstractVector)
            treated_states = [treated_states]
        end 
        unique_states = unique(data_copy[!, state])
        control_states = setdiff(unique_states, treated_states)
        if !(control_states isa AbstractVector)
            control_states = [control_states]
        end 

        # Ensure the state column is a string or number and that the nonmissingtype(treated_states) == nonmissingtype(state column)
        treated_states_type = Base.nonmissingtype(eltype(treated_states))
        state_column_type = Base.nonmissingtype(eltype(data_copy[!, state]))
        if !((treated_states_type <: Number && state_column_type <: Number) || 
            (treated_states_type <: AbstractString && state_column_type <: AbstractString))
            error("'treated_states' and the 'state' column ($state) must both be numerical or both be strings.\nInstead, found 'treated_states' $treated_states_type and '$state' $state_column_type.")
        end
        missing_states = setdiff(treated_states, data_copy.state_71X9yTx)
        if !isempty(missing_states)
            error("The following 'treated_states' could not be found in the data $(missing_states).\nOnly found the following states $(unique(data_copy.state_71X9yTx))")
        end

        # We can also filter the df down to just the treated_states
        filter!(row -> row.state_71X9yTx in treated_states, data_copy)
    end

    # Check for missing/nothing/NaN values and drop those rows
    if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy.outcome_71X9yTx)
        @warn "Found missing values in the 'outcome' column. Dropping those rows."
        data_copy = filter(row -> !(row.outcome_71X9yTx === missing || row.outcome_71X9yTx === nothing || 
                                    (row.outcome_71X9yTx isa AbstractFloat && isnan(row.outcome_71X9yTx))), data_copy)
    end 

    if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy[!, state])
        @warn "Found missing values in the 'state' column. Dropping those rows."
        data_copy = filter(row -> !(row[state] === missing || row[state] === nothing || 
                                    (row[state] isa AbstractFloat && isnan(row[state]))), data_copy)
    end 

    if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy[!, time])
        @warn "Found missing values in the 'time' column. Dropping those rows."
        data_copy = filter(row -> !(row[time] === missing || row[time] === nothing || 
                                    (row[time] isa AbstractFloat && isnan(row[time]))), data_copy)
    end
    if !(isnothing(covariates))
        for cov in covariates
            if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy[!, cov])
                @warn "Found missing values in the '$cov' column. Dropping those rows."
                data_copy = filter(row -> !(row[cov] === missing || row[cov] === nothing ||
                                            (row[cov] isa AbstractFloat && isnan(row[cov]))), data_copy)
            end 
        end
    end

    # Force outcome to float64 to speed up regression (runs faster is <:Number rather than <:Union{Number, Missing})
    data_copy.outcome_71X9yTx = convert(Vector{Float64}, data_copy.outcome_71X9yTx)

    # Check that time column entries are all in the same date format
    nonmissing_time_type = Base.nonmissingtype(eltype(data_copy[!, time]))
    if nonmissing_time_type <: AbstractString || nonmissing_time_type <:Number
        dates_str = string.(data_copy[!, time])
        if length(unique(length.(dates_str))) != 1
            error("Dates in the 'time' column are not all the same length!")
        end
        if unique(length.(dates_str))[1] != 4 && nonmissing_time_type <:Number
            error("If 'time' is a numeric column, dates must be 4 digits long.")
        end
    end
    
    if nonmissing_time_type <: AbstractString
        ref_positions, ref_sep_types = get_sep_info(dates_str[1])
        if length(ref_sep_types) > 1
            error("First date in 'time' column uses mixed separators: - and /.")
        end
        i = 0
        for date in dates_str
            i += 1
            sep_positions, sep_types = get_sep_info(date)
            if sep_positions != ref_positions
                error("$i Separator positions differs from the first date in 'time' column in date entry $i $date")
            end
            if length(sep_types) > 1
                error("$i Date found in 'time' column (entry $i) which uses multiple separator types $date")
            end
            if sep_types != ref_sep_types
                error("$i Date found in 'time' column (entry $i $date) which uses different separator types from first date.")
            end
        end
    end

    # Check that treatment_times, if strings, are entered in the same date format
    if eltype(treatment_times) <: AbstractString
        ref_positions, ref_sep_types = get_sep_info(treatment_times[1])
        if length(ref_sep_types) > 1
            error("First date in 'treatment_times' uses mixed separators - and /.")
        end
        i = 0
        for date in treatment_times
            i += 1
            sep_positions, sep_types = get_sep_info(date)
            if sep_positions != ref_positions
                error("$i Separator positions differs from the first date in 'treatment_times' in the $i'th entry $date")
            end
            if length(sep_types) > 1
                error("$i Date found in 'treatment_times' column (entry $i) which uses multiple separator types $date")
            end
            if sep_types != ref_sep_types
                error("$i Date found in 'treatment_times' (entry $i $date) which uses different separator types from first date.")
            end
        end
    end

    # Check that time column and treatment_times have the same date format if they are both strings
    if nonmissing_time_type <:AbstractString && eltype(treatment_times) <:AbstractString
        treatment_times_length = length(treatment_times[1])
        date_column_entry_length = length(data_copy[!, time][1])
        if treatment_times_length != date_column_entry_length
            error("'treatment_times' and 'time' column were found to have date strings with different lengths.")
        end
        sep_positions_time, sep_types_time = get_sep_info(data_copy[!, time][1])
        sep_positions_treatment, sep_types_treatment = get_sep_info(treatment_times[1])
        if sep_positions_time != sep_positions_treatment
            error("'treatment_times' and 'time' column have different separator positions.")
        end
        if sep_types_time != sep_types_treatment
            error("'treatment_times' and 'time' column have different separator types.")
        end
    end

    # Make sure the time column is a Date object, especially relevant for staggered adoption
    if nonmissing_time_type <: Number
        data_copy.time_71X9yTx = Date.(data_copy[!, time])
    elseif nonmissing_time_type <: AbstractString
        if isnothing(date_format)
            error("If 'time' column is a String column, must specify the 'date_format' argument.")
        end 
        data_copy.time_71X9yTx = parse_string_to_date_didint.(data_copy[!, time], date_format)
    elseif nonmissing_time_type <: Date
        data_copy.time_71X9yTx = data_copy[!, time]
    else
        error("'time' column must be a String, Date, or Number column.")
    end

    # Convert treatment_times to Date objects
    if !(typeof(treatment_times) <: Vector{Date})
        treatment_times = parse_string_to_date_didint.(treatment_times, date_format)
    end 

    # Check freq args
    if !isnothing(freq)
        freq = lowercase(freq)
        freq_options = ["week", "weeks", "weekly", "day", "days", "daily", "month", "months", "monthly", "year", "years", "yearly"]
        if !(freq in freq_options)
            error("'freq' was not set to a valid option. Try one of $freq_options")
        end 
    end 

    # Do date matching procedure
    if !isnothing(freq)
        if isnothing(start_date) 
            start_date = minimum(data_copy.time_71X9yTx)
        end
        if isnothing(end_date)
            end_date   = maximum(data_copy.time_71X9yTx)
        end
        period = parse_freq(string(freq_multiplier)*" "*freq)
        match_to_these_dates = collect(start_date:period:end_date)
        matched = [match_date(t, match_to_these_dates, treatment_times) for t in data_copy.time_71X9yTx]
        data_copy.time_71X9yTx = matched
        matched_treatment = [match_treatment_time(t, match_to_these_dates) for t in treatment_times]
        treatment_times = matched_treatment
    else
        should_you_specify_freq_and_start_end_dates = check_dates_by_state(data_copy)
    end 

    # Check that treatment_times actually exist in all_times from the data
    all_times = sort(unique(data_copy.time_71X9yTx))
    missing_dates = setdiff(treatment_times, all_times)
    if !isempty(missing_dates)
        error("The following 'treatment_times' are not found in the data $(missing_dates).\nTry defining an argument for 'freq' (and 'start_date' and 'end_date') in order to activate the date matching procedure.")
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
        if c == "hom"
            state_time = state_time[1:end - length(covariates_to_include)]
            coefs = coefs[1:end - length(covariates_to_include)]
        end
        lambda_df = DataFrame(state = first.(state_time), time = last.(state_time))
        lambda_df.lambda = coefs


        # Parse the treatment_times and the time column to dates for staggered adoption
        lambda_df.time = Date.(lambda_df.time)

        # Note the ccc for this loop
        lambda_df.ccc .= c
        master_lambda = [master_lambda;lambda_df]
    end

    time_to_index = Dict(time => idx for (idx, time) in enumerate(all_times))
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

    return event_plot_data 

    elseif event == false
        # If we're not doing an event plot, then it doesn't actually
        # matter where we assign the treat_periods, so long as we do
        treatment_times = unique(treatment_times)
        master_lambda.treat_period = [date in treatment_times ? time_to_index[date] : missing for date in master_lambda.time] .- 1
        master_lambda.time = parse_date_to_string_didint.(master_lambda.time, date_format)
        return master_lambda
    end

end