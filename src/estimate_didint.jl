"""
    didint(outcome::AbstractString, state::AbstractString,
           time::AbstractString, data::DataFrame,
           treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number},
           treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date};
           date_format::Union{AbstractString, Nothing} = nothing,
           covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
           ccc::AbstractString = "int", agg::AbstractString = "state",
           ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
           freq::Union{AbstractString, Nothing} = nothing, freq_multiplier::Number = 1,
           autoadjust::Bool = false, nperm::Number = 1000, verbose::Bool = true,
           seed::Number = rand(1:1000000))

The `didint()` function estimates the average effect of treatment on the treated (ATT)
while accounting for covariates that may vary by state, time, or by both state and time simultaneously.

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
function didint(outcome::Union{AbstractString, Symbol},
                state::Union{AbstractString, Symbol},
                time::Union{AbstractString, Symbol},
                data::DataFrame;
                gvar::Union{AbstractString, Symbol, Nothing} = nothing,
                treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
                treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
                date_format::Union{AbstractString, Nothing} = nothing,
                covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
                ccc::AbstractString = "int",
                agg::AbstractString = "cohort",
                weighting::AbstractString = "both",
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                nperm::Number = 1001,
                verbose::Bool = true,
                seed::Number = rand(1:1000000),
                use_pre_controls::Bool = false)

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
        if isnothing(treatment_times) || isnothing(treated_states)
            error("If 'gvar' is not specified, then 'treatment_times' and 'treated_states' must be specified!")
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

    # Check that seed is set correctly
    seed = Int(round(seed))
    if seed <= 0 
        error("'seed' must be a positive integer.")
    end 

    # Check that agg args are passed correctly
    agg = lowercase(agg)
    agg_options = ["cohort", "state", "simple", "none", "sgt", "time"]
    if !(agg in agg_options)
        error("'agg' must be one of $(agg_options)")
    end 

    # Check that weighting arg is passed correctly
    weighting = lowercase(weighting)
    weighting_options = ["both", "none", "att", "diff"]
    if !(weighting in weighting_options) 
        error("'weighting' must be one of $(weighting_options)")
    end

    # Round nperm
    nperm = Int(round(nperm))
    if nperm < 1
        error("'nperm' must be a positive integer > 0.")
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
            error("The following covariates could not be found in the data $(join(missing_cov, ", ")) ")
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
    if length(unique(treatment_times)) == 1
        common_adoption = true
        staggered_adoption = false
    elseif length(unique(treatment_times)) > 1 
        common_adoption = false
        staggered_adoption = true
    else
        error("'treatment_times' must have at least one entry.")
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
    missing_control_states = isempty(control_states)
    if missing_control_states && !use_pre_controls
        error("No control states were found.")
    elseif missing_control_states && use_pre_controls
        # If there are at least two unique treatment times, then should still be able to compute
        # some ATTs if using the notyet treated cells as controls 
        if !staggered_adoption
            error("No valid control states were found.")
        end 
    end

    # Check missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # Ensure the state column is a string or number and that the nonmissingtype(treated_states) == nonmissingtype(state column)
    treated_states_type = Base.nonmissingtype(eltype(treated_states))
    state_column_type = Base.nonmissingtype(eltype(data_copy[!, state]))
    if !((treated_states_type <: Number && state_column_type <: Number) || 
        (treated_states_type <: AbstractString && state_column_type <: AbstractString))
        error("'treated_states' and the 'state' column ($state) must both be numerical or both be strings. \n 
Instead, found 'treated_states' $treated_states_type and '$state' $state_column_type.")
    end
    data_copy.state_71X9yTx = data_copy[!, state]

    # Do some checks after (potentially) dropping rows to make sure that the analysis
    # is still valid
    if (length(unique(data_copy.state_71X9yTx)) == 1)
        error("Not enough non-missing observations.\nFound", length(unique(data_copy.state_71X9yTx)), "states after dropping missing observations.")
    end 

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

    # Grab all existing times in the data
    all_times = sort(unique(data_copy.time_71X9yTx))
    if isnothing(start_date)
        start_date = minimum(all_times)
    end
    if isnothing(end_date)
        end_date = maximum(all_times)
    end

    data_copy = filter(row -> row.time_71X9yTx >= start_date && row.time_71X9yTx <= end_date, data_copy)

    if nrow(data_copy) == 0
        error("After filtering by 'start_date' and 'end_date' found only 0 rows of data!")
    end 

    # In the case of staggered adoption, check if date matching procedure should be done
    if staggered_adoption
        if !isnothing(freq)
            period = parse_freq(string(freq_multiplier)*" "*freq)
        else
            period = get_freq_approx(all_times)
        end
        match_to_these_dates = collect(start_date:period:end_date)
        one_past = end_date + period
        matched = [match_date(t, match_to_these_dates, treatment_times) for t in data_copy.time_71X9yTx]
        data_copy.time_71X9yTx = matched
        treated_states = treated_states[(treatment_times .< one_past) .&& (treatment_times .> start_date)]
        treatment_times = treatment_times[(treatment_times .< one_past) .&& (treatment_times .> start_date)]
        treatment_times = [match_treatment_time(t, match_to_these_dates) for t in treatment_times]
    end

    # Check that treatment_times actually exist in all_times from the data
    if staggered_adoption
        all_times = sort(unique(data_copy.time_71X9yTx))
        missing_dates = setdiff(treatment_times, all_times)
        if !isempty(missing_dates)
            error("The following 'treatment_times' are not found in the data $(missing_dates). \n 
Try defining an argument for 'freq' (and 'start_date' and 'end_date') in order to activate the date matching procedure.")
        end
    end 

    # Check that treated states actually exist
    missing_states = setdiff(treated_states, data_copy.state_71X9yTx)
    if !isempty(missing_states)
        error("The following 'treated_states' could not be found in the data $(missing_states). \n 
Only found the following states $(unique(data_copy.state_71X9yTx))")
    end

    # Do some checks for treatment_times vector length
    if common_adoption && length(treatment_times) != length(treated_states)
        if length(treatment_times) != 1
            error("'treatment_times' should either be the same length as the 'treated_states' vector or of length 1.")
        else
            treatment_times = fill(treatment_times[1], length(treated_states))
        end 
    elseif staggered_adoption && length(treatment_times) != length(treated_states)
        error("'treatment_times' should be the same length as the 'treated_states'.")
    end 

    # Also need to make sure that start_times < treatment_times < end_times is true for each state
    for i in eachindex(treated_states)
        s = treated_states[i]
        treat_time = treatment_times[i]
        state_dates = data_copy[data_copy.state_71X9yTx .== s, :time_71X9yTx]
        earliest = minimum(state_dates)
        latest = maximum(state_dates)
        if !(earliest < treat_time)
            error("For state $s, the earliest date ($earliest) is not strictly less than the treatment time ($treat_time).")
        end
        if !(treat_time <= latest)
            error("For state $s, the treatment time ($treat_time) is greater than the last date ($latest).")
        end
    end

    # Ensure state column is a string (as are treated_states)
    data_copy.state_71X9yTx = string.(data_copy.state_71X9yTx)
    treated_states = string.(treated_states)
    check_states = unique(data_copy.state_71X9yTx)
    missing_states = setdiff(treated_states, check_states)
    if !isempty(missing_states)
        error("The states $missing_states could not be found among the states in the data $check_states")
    end

    # Once any date matching procedures are done, convert time_71X9yTx back to a string for processing in `categorical()`
    # Also keep a column vector copy as a date
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    if common_adoption
        data_copy.time_71X9yTx = ifelse.(data_copy.time_71X9yTx .>= treatment_times[1], "post", "pre")
        data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)
        if length(unique(data_copy.time_71X9yTx)) == 1
            error("Only", unique(data_copy.time_71X9yTx), "treatment periods were found in the data!")
        end
    elseif staggered_adoption
        data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)
    else
        error("A non-common or non-staggered adoption scenario was discovered!?")
    end 

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

    # Construct formula, depending on DID-INT variation
    formula_str = "@formula(outcome_71X9yTx ~ 0 + state_time"

    # Set up formula
    ccc = lowercase(ccc)
    ccc = replace(ccc, r"\s" => "")
    if ccc == "int"
        for c in covariates_to_include
            formula_str *= " + fe(state_71X9yTx)&fe(time_71X9yTx)&$c"
        end
    elseif ccc == "time"
        for c in covariates_to_include
            formula_str *= " + fe(time_71X9yTx)&$c"
        end
    elseif ccc == "state"
        for c in covariates_to_include
            formula_str *= " + fe(state_71X9yTx)&$c"
        end
    elseif ccc == "add"
        for c in covariates_to_include
            formula_str *= " + fe(time_71X9yTx)&$c + fe(state_71X9yTx)&$c"
        end
    elseif ccc == "hom"
        for c in covariates_to_include
            formula_str *= " + $c"
        end
    else 
        error("'ccc' must be set to one of \"int\", \"time\", \"state\", \"add\", or \"hom\".")
    end 
    formula_str *= ")"
    formula_expr = Meta.parse(formula_str)
    formula = eval(formula_expr)

    # Call GC.gc() before running big FixedEffectsModels regression
    GC.gc()

    # Determine if corner case 
    cornercase = false
    if common_adoption && length(unique(data_copy.state_71X9yTx)) == 2
        cornercase = true
    end 

    # Run FixedEffectsModels regression (stage 1)
    if cornercase
        stage1 = reg(data_copy, formula, Vcov.robust(), contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
                     save = false)
    else
        stage1 = reg(data_copy, formula, contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
                     save = false)
    end

    # Recover lambdas
    state_time = split.(replace.(coefnames(stage1), "state_time: " => ""), "0IQR7q6Wei7Ejp4e")
    coefs = coef(stage1)
    if ccc == "hom"
        state_time = state_time[1:end - length(covariates_to_include)]
        coefs = coefs[1:end - length(covariates_to_include)]
    end 
    lambda_df = DataFrame(state = first.(state_time), time = last.(state_time))
    lambda_df.lambda = coefs
    if cornercase
        vcov_matrix = vcov(stage1)[1:4,1:4]

        treat_post_var = vcov_matrix[1,1]
        treat_pre_var = vcov_matrix[2,2]
        treat_cov = vcov_matrix[2,1]
        control_post_var = vcov_matrix[3,3]
        control_pre_var = vcov_matrix[4,4]
        control_cov = vcov_matrix[4,3]

        cornercase_se = sqrt(treat_post_var + treat_pre_var + control_post_var + control_pre_var - 2*treat_cov - 2*control_cov)

    end 

    # Compute diff for each treated state
    unique_states = unique(lambda_df.state)
    diff_df = DataFrame(state = String[], treated_time = Date[],
                        t = Date[], r1 = Date[], diff = Float64[],
                        treat = Int[], n = Int[], n_t = Int[])

    # Parse the treatment_times and the time column to dates for staggered adoption
    if staggered_adoption
        lambda_df.time = Date.(lambda_df.time)
        time_to_index = Dict(time => idx for (idx, time) in enumerate(all_times))
    end    

    # Define a function that initializes the results dataframe columns
    init_column() = Vector{Union{Missing, Float64}}(missing, nrows)

    # Compute ATT
    if staggered_adoption
        # Compute diffs for each treated_state
        for i in eachindex(treated_states)
            idx = time_to_index[treatment_times[i]]
            one_period_prior_treatment = all_times[idx - 1]
            temp = lambda_df[(lambda_df.state .== treated_states[i]) .& (lambda_df.time .>= one_period_prior_treatment), :]
            lambda_r1 = temp[temp.time .== one_period_prior_treatment, "lambda"][1]
            temp = temp[temp.time .> one_period_prior_treatment, :]
            sort!(temp, :time)
            years = temp.time
            diffs = Vector{Float64}(undef, nrow(temp))
            n = Vector{Int}(undef, nrow(temp))
            n_t = Vector{Int}(undef, nrow(temp))
            for j in eachindex(diffs)
                diffs[j] = temp[j, "lambda"][1] - lambda_r1
                post = temp[j, "time"]
                pre = one_period_prior_treatment
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([post, pre])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
                n_t[j] = sum((data_copy.time_dmG5fpM .== post) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            end 
            temp_df = DataFrame(state = treated_states[i], treated_time = treatment_times[i],
                                t = years, r1 = one_period_prior_treatment, diff = diffs, treat = 1,
                                n = n, n_t = n_t)
            diff_df = vcat(diff_df, temp_df) 
        end
        
        # Compute diffs for control states
        unique_diffs = unique(select(diff_df, :t, :r1, :treated_time))
        control_states = setdiff(unique_states, treated_states)     
        for i in eachindex(control_states)
            temp = lambda_df[(lambda_df.state .== control_states[i]), :]
            diffs = Vector{Float64}(undef, nrow(unique_diffs))
            n = Vector{Int}(undef, nrow(unique_diffs))
            n_t = Vector{Int}(undef, nrow(unique_diffs))
            for j in 1:nrow(unique_diffs)
                t = unique_diffs[j,"t"]
                r1 = unique_diffs[j,"r1"]
                diffs[j] = temp[temp.time .== t, "lambda"][1] - temp[temp.time .== r1, "lambda"][1]
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== control_states[i]))
                n_t[j] = sum((data_copy.time_dmG5fpM .== t) .&& (data_copy.state_71X9yTx .== control_states[i]))
            end 
            trtd_time = [all_times[findfirst(==(t), all_times) + 1] for t in unique_diffs.r1]
            temp_df = DataFrame(state = control_states[i], treated_time = trtd_time,
                                t = unique_diffs.t, r1 = unique_diffs.r1, diff = diffs, treat = 0,
                                n = n, n_t = n_t)
            diff_df = vcat(diff_df, temp_df)
        end

        # Compute diffs for randomization inference
        ri_diff_df = DataFrame(state = String[], treated_time = Date[],
                               t = Date[], r1 = Date[], diff = Float64[],
                               treat = Int[], n = Int[], n_t = Int[])
        for i in eachindex(treated_states)
            ri_diffs = unique(select(diff_df[(diff_df.state .!= treated_states[i]) .& (diff_df.treated_time .!= treatment_times[i]), :], :t, :r1))
            temp = lambda_df[(lambda_df.state .== treated_states[i]), :]
            diffs = Vector{Float64}(undef, nrow(ri_diffs))
            n = Vector{Int}(undef, nrow(ri_diffs))
            n_t = Vector{Int}(undef, nrow(ri_diffs))
            for j in 1:nrow(ri_diffs)
                t = unique_diffs[j,"t"]
                r1 = unique_diffs[j,"r1"]
                diffs[j] = temp[temp.time .== t, "lambda"][1] - temp[temp.time .== r1, "lambda"][1]
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
                n_t[j] = sum((data_copy.time_dmG5fpM .== t) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            end
            trtd_time = [all_times[findfirst(==(t), all_times) + 1] for t in ri_diffs.r1]
            temp_df = DataFrame(state = treated_states[i], treated_time = trtd_time,
                                t = ri_diffs.t, r1 = ri_diffs.r1, diff = diffs, treat = -1,
                                n = n, n_t = n_t)
            ri_diff_df = vcat(ri_diff_df, temp_df)
        end

        # Add the periods since treatment column to ri_diff_df and diff_df
        ordered_t = sort(unique(diff_df.t))        
        t_period = indexin(diff_df.t, ordered_t)             
        treated_time_period = indexin(diff_df.treated_time, ordered_t) 
        diff_df.time_since_treatment = t_period .- treated_time_period 
        t_period = indexin(ri_diff_df.t, ordered_t)             
        treated_time_period = indexin(ri_diff_df.treated_time, ordered_t) 
        ri_diff_df.time_since_treatment = t_period .- treated_time_period 

        # If use_pre_controls is true (default) re-structure the ri_diff_df and diff_df
        # to add the pre-treatment periods from treated silos as control contrasts
        if use_pre_controls
            for (t,s) in zip(treatment_times, treated_states)
                ri_diff_df[ri_diff_df.state .== s .&& ri_diff_df.t .< t, "treat"] .= 0
            end
            diff_df = vcat(diff_df, ri_diff_df[ri_diff_df.treat .== 0,:])
            ri_diff_df = ri_diff_df[ri_diff_df.treat .== -1,:]
            if missing_control_states
                # If there are no pure control states, do extra processing to remove cells
                # where there are no treat-control pairs
                diff_df = subset(
                            groupby(diff_df, [:t, :r1]),
                            :treat => x -> any(x .== 1) && any(x .== 0)
                          )
                ri_diff_df = semijoin(ri_diff_df, diff_df, on=[:t, :r1])
            end 
        end

        # Run final regression to compute ATT based on weighting/aggregation method
        if agg == "cohort"

            # Define nrows and vector to iterate through for sub aggregate ATTs
            unique_treatment_times = sort!(unique(diff_df.treated_time))
            nrows = length(unique_treatment_times)

            # Create results df
            results = DataFrame(treatment_time = Vector{Date}(undef, nrows),
                                att_cohort = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_cohort = init_column(),
                                pval_att_cohort = init_column(),
                                jknifese_att_cohort = init_column(),
                                jknifepval_att_cohort = init_column(),
                                ri_pval_att_cohort = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Compute the sub aggregate ATTs
            for i in eachindex(unique_treatment_times)
                trt = unique_treatment_times[i]
                temp = diff_df[diff_df.treated_time .== trt, :]
                X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end            
                results.treatment_time[i] = trt
                result_dict = final_regression_results(X, Y, W = W)
                results.att_cohort[i] = result_dict["beta_hat"]
                results.se_att_cohort[i] = result_dict["beta_hat_se"]
                results.pval_att_cohort[i] = result_dict["pval_att"] 
                results.jknifese_att_cohort[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_cohort[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Compute aggregate ATT and return results
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_cohort)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "cohort",
                                                      verbose, seed, data_copy, weighting, use_pre_controls)
            return results

        elseif agg == "simple"

            # Define nrows and object to iterate through for sub aggregate ATTs (unique_diffs)
            nrows = nrow(unique_diffs)

            # Create results df
            results = DataFrame(r1 = Vector{Date}(undef, nrows),
                                time = Vector{Date}(undef, nrows),
                                gvar = Vector{Date}(undef, nrows),
                                att_gt = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_gt = init_column(),
                                pval_att_gt = init_column(),
                                jknifese_att_gt = init_column(),
                                jknifepval_att_gt = init_column(),
                                ri_pval_att_gt = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))
            
            # Compute the sub aggregate ATTs
            for i in 1:nrow(unique_diffs)
                t = unique_diffs[i,"t"]
                r1 = unique_diffs[i,"r1"]
                gvar = unique_diffs[i, "treated_time"]
                temp = diff_df[(diff_df.t .== t) .& (diff_df.r1 .== r1), :]
                X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.time[i] = t
                results.r1[i] = r1
                results.gvar[i] = gvar
                result_dict = final_regression_results(X, Y, W = W)
                results.att_gt[i] = result_dict["beta_hat"]
                results.se_att_gt[i] = result_dict["beta_hat_se"]
                results.pval_att_gt[i] = result_dict["pval_att"] 
                results.jknifese_att_gt[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_gt[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Sort results dataframe and compute aggregate ATT, return results
            sort!(results, [:time])
            sort!(results, [:r1])
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_gt)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "simple",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)
            return results

        elseif agg == "state"

            # Define nrows and identify object to iterate thru for the sub aggregate ATTs
            nrows = length(treated_states)

            # Create results df
            results = DataFrame(state = Vector{String}(undef, nrows),
                                att_s = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_s = init_column(),
                                pval_att_s = init_column(),
                                jknifese_att_s = init_column(),
                                jknifepval_att_s = init_column(),
                                ri_pval_att_s = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Compute sub aggregate ATTs
            for i in eachindex(treated_states)
                state = treated_states[i]
                trt = treatment_times[i]
                temp_treated = diff_df[diff_df.state .== state, :]
                temp_control = diff_df[(diff_df.treat .== 0) .&& (diff_df.treated_time .== trt), :]
                temp = vcat(temp_control, temp_treated)
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.state[i] = state
                result_dict = final_regression_results(X, Y, W = W)
                results.att_s[i] = result_dict["beta_hat"]
                results.se_att_s[i] = result_dict["beta_hat_se"]
                results.pval_att_s[i] = result_dict["pval_att"] 
                results.jknifese_att_s[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_s[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Sort results dataframe and compute aggregate ATT, return results
            results.tuple_state = custom_sort_order.(results.state)
            sort!(results,[order(:tuple_state)])
            select!(results, Not([:tuple_state]))
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_s)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "state",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)
            return results

        elseif agg == "none"  
            
            # Create results df
            results = DataFrame(agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                nperm = Vector{Union{Missing, Float64}}(missing, 1))
            
            # Compute aggregate ATT, return results
            X = convert(Matrix{Float64}, hcat(ones(nrow(diff_df)), diff_df.treat))
            Y = convert(Vector{Float64}, diff_df.diff)
            if in(weighting, ["diff", "both"])
                W = convert(Vector{Float64}, diff_df.n)
                W ./= sum(W)
            else
                W = fill(nothing, nrow(diff_df))
            end 
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "none",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)
            return results

        elseif agg == "sgt"
            
            # Compute nrows and define object to iterate thru for the sub aggregate ATTs
            unique_sgt = unique(select(diff_df[diff_df.treat .== 1,:], :state, :t, :treated_time))
            nrows = nrow(unique_sgt)

            # Create results df
            results = DataFrame(state = Vector{String}(undef, nrows),
                                gvar = Vector{Date}(undef, nrows),
                                t = Vector{Date}(undef, nrows),
                                att_sgt = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_sgt = init_column(),
                                pval_att_sgt = init_column(),
                                jknifese_att_sgt = init_column(),
                                jknifepval_att_sgt = init_column(),
                                ri_pval_att_sgt = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))
            
            # Compute the sub aggregate ATTs 
            for i in 1:nrow(unique_sgt)
                state = unique_sgt[i, "state"]
                gvar = unique_sgt[i, "treated_time"]
                t = unique_sgt[i, "t"]
                temp_treated = diff_df[(diff_df.state .== state) .&& (diff_df.treated_time .== gvar) .&& (diff_df.t .== t),:]
                temp_control = diff_df[(diff_df.treat .== 0) .&& (diff_df.treated_time .== gvar) .&& (diff_df.t .== t),:]
                temp = vcat(temp_treated, temp_control)
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.state[i] = state
                results.gvar[i] = gvar
                results.t[i] = t
                result_dict = final_regression_results(X, Y, W = W)
                results.att_sgt[i] = result_dict["beta_hat"]
                results.se_att_sgt[i] = result_dict["beta_hat_se"]
                results.pval_att_sgt[i] = result_dict["pval_att"] 
                results.jknifese_att_sgt[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_sgt[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end 

            # Compute aggregate ATT, return results
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_sgt)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "sgt",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)
            return results

        elseif agg == "time"

            # Compute nrows, define times to iterate thru to get sub aggregate ATTs
            nrows = length(unique(diff_df.time_since_treatment))
            times = sort(unique(diff_df.time_since_treatment))

            # Create results df
            results = DataFrame(periods_post_treat = init_column(),
                                att_t = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_t = init_column(),
                                pval_att_t = init_column(),
                                jknifese_att_t = init_column(),
                                jknifepval_att_t = init_column(),
                                ri_pval_att_t = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Special procedure for adding time or gvar dummies (yields equivalent results)
            diff = vcat(diff_df, ri_diff_df)
            treated_category = categorical(string.(diff.treated_time); ordered = true)
            levels_tt   = levels(treated_category)         
            dummy_cols  = Symbol[]                         
            for i in eachindex(levels_tt) 
                tt = levels_tt[i]                    
                colname = Symbol("tt_", i)
                diff[!, colname] = Int.(treated_category .== tt)
                push!(dummy_cols, colname)
            end

            # Compute sub aggregate ATTs
            for i in eachindex(times)
                t = times[i]
                temp = diff[(diff.time_since_treatment .== t) .&& (diff.treat .!= -1), :]
                X = design_matrix_time_agg(temp, dummy_cols, Symbol("treat"))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.periods_post_treat[i] = t
                result_dict = final_regression_results(X, Y, W = W)
                results.att_t[i] = result_dict["beta_hat"]
                results.se_att_t[i] = result_dict["beta_hat_se"]
                results.pval_att_t[i] = result_dict["pval_att"] 
                results.jknifese_att_t[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_t[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Compute aggregate ATT, return results
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_t)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(diff, nperm, results, "time",
                                                 verbose, seed, data_copy, weighting, use_pre_controls,
                                                 dummy_cols = dummy_cols)
            return results
        end

    elseif common_adoption

        # Create diff matrix
        states = unique(lambda_df.state)
        diff_df = DataFrame(state = states,
                            diff = Vector{Float64}(undef, length(states)),
                            treat = Vector{Float64}(undef, length(states)),
                            treated_time = fill(unique(treatment_times)[1], length(states)),
                            n = Vector{Float64}(undef, length(states)),
                            n_t = Vector{Float64}(undef, length(states)))
        for i in eachindex(states)
            state = states[i]
            if state in treated_states
                trt = 1
            else 
                trt = 0
            end 
            diff = lambda_df[(lambda_df.state .== state) .&& (lambda_df.time .== "post"), "lambda"][1] - lambda_df[(lambda_df.state .== state) .&& (lambda_df.time .== "pre"), "lambda"][1]
            diff_df.diff[i] = diff
            diff_df.treat[i] = trt
            diff_df.n[i] = sum(data_copy.state_71X9yTx .== state)
            diff_df.n_t[i] = sum((data_copy.state_71X9yTx .== state) .&& (data_copy.time_71X9yTx .== "post"))
        end 

        # Compute results
        if in(agg, ["cohort", "gt", "none"]) || length(treated_states) == 1
            results = DataFrame(agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                nperm = Vector{Union{Missing, Float64}}(missing, 1))
            X = convert(Matrix{Float64}, hcat(ones(nrow(diff_df)), diff_df.treat))
            Y = convert(Vector{Float64}, diff_df.diff)
            if in(weighting, ["diff", "both"])
                W = convert(Vector{Float64}, diff_df.n)
                W ./= sum(W)
            else
                W = fill(nothing, nrow(diff_df))
            end 
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            if cornercase
                results.se_agg_att[1] = cornercase_se
            else
                results.se_agg_att[1] = result_dict["beta_hat_se"]
            end 
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(diff_df, nperm, results, agg,
                                                 verbose, seed, data_copy, weighting,
                                                 use_pre_controls)

            return results
        elseif in(agg, ["state", "sgt"])

            # Define nrows and identify object to iterate thru for the sub aggregate ATTs
            nrows = length(treated_states)

            # Create results df
            results = DataFrame(state = Vector{String}(undef, nrows),
                                att_s = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_s = init_column(),
                                pval_att_s = init_column(),
                                jknifese_att_s = init_column(),
                                jknifepval_att_s = init_column(),
                                ri_pval_att_s = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Compute sub aggregate ATTs
            for i in eachindex(treated_states)
                state = treated_states[i]
                trt = treatment_times[i]
                temp_treated = diff_df[diff_df.state .== state, :]
                temp_control = diff_df[(diff_df.treat .== 0) .&& (diff_df.treated_time .== trt), :]
                temp = vcat(temp_control, temp_treated)
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.state[i] = state
                result_dict = final_regression_results(X, Y, W = W)
                results.att_s[i] = result_dict["beta_hat"]
                results.se_att_s[i] = result_dict["beta_hat_se"]
                results.pval_att_s[i] = result_dict["pval_att"] 
                results.jknifese_att_s[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_s[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Sort results dataframe and compute aggregate ATT, return results
            results.tuple_state = custom_sort_order.(results.state)
            sort!(results,[order(:tuple_state)])
            select!(results, Not([:tuple_state]))
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_s)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(diff_df, nperm, results, "state",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)
            return results

        end 
    end 



end 