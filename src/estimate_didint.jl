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
function didint(outcome::AbstractString,
                state::AbstractString,
                time::AbstractString,
                data::DataFrame,
                treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number},
                treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date};
                date_format::Union{AbstractString, Nothing} = nothing,
                covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
                ccc::AbstractString = "int", agg::AbstractString = "cohort",
                weighting::AbstractString = "att",
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                autoadjust::Bool = false,
                nperm::Number = 1001,
                verbose::Bool = true,
                seed::Number = rand(1:1000000))

    # Check that seed is set correctly
    seed = Int(round(seed))
    if seed <= 0 
        error("Er47: 'seed' must be a positive integer.")
    end 

    # Check that agg args are passed correctly
    agg = lowercase(agg)
    agg_options = ["cohort", "state", "simple", "none", "sgt"]
    if !(agg in agg_options)
        error("Er03: 'agg' must be one of: $(agg_options)")
    end 

    # Check that weighting arg is passed correctly
    weighting = lowercase(weighting)
    weighting_options = ["both", "none", "att", "diff"]
    if !(weighting in weighting_options) 
        error("Er99: 'weighting must be one of: $(weighting_options)'")
    end

    # Round nperm
    nperm = Int(round(nperm))
    if nperm < 1
        error("Er04: 'nperm' must be a positive integer > 0.")
    end 

    # Create data copy
    data_copy = DataFrame(data; copycols = false)    

    # Ensure the specified outcome, state, and time columns exist in the data
    missing_cols = [col for col in [outcome, state, time] if !(col in names(data_copy))]
    if !isempty(missing_cols)
        error("Er05: The following columns could not be found in the data: ", join(missing_cols, ", "))
    end

    # Ensure the outcome variable is a numeric variable and create outcome column
    outcome_nonmissingtype = Base.nonmissingtype(eltype(data_copy[!, outcome]))
    if !(outcome_nonmissingtype <: Number)
        error("Er06: Column '$outcome' must be numeric, but found $(eltype(data_copy[!, outcome]))")
    end
    data_copy.outcome_71X9yTx = data_copy[!, outcome]

    # Check that the specified covariates exist
    if !isnothing(covariates)
        if isa(covariates, AbstractString)
            covariates = [covariates]
        end
        missing_cov = [col for col in covariates if !(col in names(data_copy))]
        if !isempty(missing_cov)
            error("Er07$(join(missing_cov, ", ")): The preceding covariates could not be found in the data.")
        end
    end

    # Ensure that treatment_times, if a number, are all 4 digit entries
    if eltype(treatment_times) <: Number
        treatment_times_numeric = true
        if isnothing(date_format) 
            date_format = "yyyy"
        end 
        if lowercase(date_format) != "yyyy"
            error("Er08: If 'treatment_times' are entered as a number, the 'date_format' must be \"yyyy\".")
        end 
        if treatment_times isa AbstractVector
            unique_lengths = unique(length.(string.(treatment_times)))
        else 
            unique_lengths = [length(string(treatment_times))]
            treatment_times = [treatment_times]
        end
        if length(unique_lengths) == 1
            if unique_lengths[1] != 4
                error("Er09: If 'treatment_times' are entered as numbers, they must all be 4 digits long in 'yyyy' date_format.")
            end 
        else
            error("Er10: Detected multiple unique date_formats in 'treatment_times'.")
        end
        treatment_times = Date.(treatment_times)
    else 
        treatment_times_numeric = false
    end
    nonmissing_time_type = Base.nonmissingtype(eltype(data_copy[!, time]))
    if nonmissing_time_type <: Number
        time_column_numeric = true
        unique_time_lengths = unique(length.(string.(data_copy[!, time])))
        if length(unique_time_lengths) > 1 || unique_time_lengths[1] != 4
            error("Er11: The 'time' column was found to be numeric but consisting of values of ambiguous date formatting (i.e. not consistent 4 digit entries.)")
        end 
    else
        time_column_numeric = false
    end 
    if xor(time_column_numeric, treatment_times_numeric)
        error("Er12: If 'time' column is numeric or 'treatment_times' is numeric, then both must be numeric.")
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
        error("Er13: 'treatment_times' must have at least one entry.")
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
    if missing_control_states
        error("Er14: No control states were found.")
    end 

    # Ensure the state column is a string or number and that the nonmissingtype(treated_states) == nonmissingtype(state column)
    treated_states_type = Base.nonmissingtype(eltype(treated_states))
    state_column_type = Base.nonmissingtype(eltype(data_copy[!, state]))
    if !((treated_states_type <: Number && state_column_type <: Number) || 
        (treated_states_type <: AbstractString && state_column_type <: AbstractString))
        error("Er15: 'treated_states' and the 'state' column ($state) must both be numerical or both be strings. \n 
Instead, found: 'treated_states': $treated_states_type and '$state': $state_column_type.")
    end
    data_copy.state_71X9yTx = data_copy[!, state]
    missing_states = setdiff(treated_states, data_copy.state_71X9yTx)
    if !isempty(missing_states)
        error("Er16: The following 'treated_states' could not be found in the data: $(missing_states). \n 
Only found the following states: $(unique(data_copy.state_71X9yTx))")
    end

    # Check for missing/nothing/NaN values
    if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy.outcome_71X9yTx)
        error("Er17: Found missing values in the 'outcome' column.")
    end 
    if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy[!, state])
        error("Er18: Found missing values in the 'state' column.")
    end 
    if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy[!, time])
        error("Er19: Found missing values in the 'time' column.")
    end 
    if !(isnothing(covariates))
        for cov in covariates
            if any(x -> x === missing || x === nothing || (x isa AbstractFloat && isnan(x)), data_copy[!, cov])
                error("Er20: Found missing values in the '$cov' column.")
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
            error("Er21: Dates in the 'time' column are not all the same length!")
        end
        if unique(length.(dates_str))[1] != 4 && nonmissing_time_type <:Number
            error("Er22: If 'time' is a numeric column, dates must be 4 digits long.")
        end
    end
    
    if nonmissing_time_type <: AbstractString
        ref_positions, ref_sep_types = get_sep_info(dates_str[1])
        if length(ref_sep_types) > 1
            error("Er23: First date in 'time' column uses mixed separators: - and /.")
        end
        i = 0
        for date in dates_str
            i += 1
            sep_positions, sep_types = get_sep_info(date)
            if sep_positions != ref_positions
                error("Er24 $i: Separator positions differs from the first date in 'time' column in date entry $i: $date")
            end
            if length(sep_types) > 1
                error("Er25 $i: Date found in 'time' column (entry $i) which uses multiple separator types: $date")
            end
            if sep_types != ref_sep_types
                error("Er26 $i: Date found in 'time' column (entry $i: $date) which uses different separator types from first date.")
            end
        end
    end

    # Check that treatment_times, if strings, are entered in the same date format
    if eltype(treatment_times) <: AbstractString
        ref_positions, ref_sep_types = get_sep_info(treatment_times[1])
        if length(ref_sep_types) > 1
            error("Er27: First date in 'treatment_times' uses mixed separators: - and /.")
        end
        i = 0
        for date in treatment_times
            i += 1
            sep_positions, sep_types = get_sep_info(date)
            if sep_positions != ref_positions
                error("Er28 $i: Separator positions differs from the first date in 'treatment_times' in the $i'th entry: $date")
            end
            if length(sep_types) > 1
                error("Er29 $i: Date found in 'treatment_times' column (entry $i) which uses multiple separator types: $date")
            end
            if sep_types != ref_sep_types
                error("Er30 $i: Date found in 'treatment_times' (entry $i: $date) which uses different separator types from first date.")
            end
        end
    end

    # Check that time column and treatment_times have the same date format if they are both strings
    if nonmissing_time_type <:AbstractString && eltype(treatment_times) <:AbstractString
        treatment_times_length = length(treatment_times[1])
        date_column_entry_length = length(data_copy[!, time][1])
        if treatment_times_length != date_column_entry_length
            error("Er31: 'treatment_times' and 'time' column were found to have date strings with different lengths.")
        end
        sep_positions_time, sep_types_time = get_sep_info(data_copy[!, time][1])
        sep_positions_treatment, sep_types_treatment = get_sep_info(treatment_times[1])
        if sep_positions_time != sep_positions_treatment
            error("Er32: 'treatment_times' and 'time' column have different separator positions.")
        end
        if sep_types_time != sep_types_treatment
            error("Er33: 'treatment_times' and 'time' column have different separator types.")
        end
    end

    # Make sure the time column is a Date object, especially relevant for staggered adoption
    if nonmissing_time_type <: Number
        data_copy.time_71X9yTx = Date.(data_copy[!, time])
    elseif nonmissing_time_type <: AbstractString
        if isnothing(date_format)
            error("Er34: If 'time' column is a String column, must specify the 'date_format' argument.")
        end 
        data_copy.time_71X9yTx = parse_string_to_date_didint.(data_copy[!, time], date_format)
    elseif nonmissing_time_type <: Date
        data_copy.time_71X9yTx = data_copy[!, time]
    else
        error("Er35: 'time' column must be a String, Date, or Number column.")
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
            error("Er36: 'freq' was not set to a valid option. Try one of: $freq_options")
        end 
    end 

    # In the case of staggered adoption, check if date matching procedure should be done
    if staggered_adoption
        if !isnothing(freq) || autoadjust
            start_date = minimum(data_copy.time_71X9yTx)
            end_date   = maximum(data_copy.time_71X9yTx)
            if autoadjust
                intervals = diff([start_date; sort(treatment_times); end_date])
                days_int = map(Dates.value, intervals)
                med = median(days_int)
                mea = mean(days_int)
                date_interval = Day(round(mean([med, mea])))
                match_to_these_dates = collect(start_date:date_interval:end_date)
            elseif !isnothing(freq)
                period = parse_freq(string(freq_multiplier)*" "*freq)
                match_to_these_dates = collect(start_date:period:end_date)
            end
            matched = [match_date(t, match_to_these_dates, treatment_times) for t in data_copy.time_71X9yTx]
            data_copy.time_71X9yTx = matched
            matched_treatment = [match_treatment_time(t, match_to_these_dates) for t in treatment_times]
            treatment_times = matched_treatment
        end 
    end

    # Check that treatment_times actually exist in all_times from the data
    if staggered_adoption
        all_times = sort(unique(data_copy.time_71X9yTx))
        missing_dates = setdiff(treatment_times, all_times)
        if !isempty(missing_dates)
            error("Er37: The following 'treatment_times' are not found in the data: $(missing_dates). \n 
Try defining an argument for 'freq' or set 'autoadjust = true' in order to activate the date matching procedure.")
        end
    end 

    # Do some checks for treatment_times vector length
    if common_adoption && length(treatment_times) != length(treated_states)
        if length(treatment_times) != 1
            error("Er38: 'treatment_times' should either be the same length as the 'treated_states' vector or of length 1.")
        else
            treatment_times = fill(treatment_times[1], length(treated_states))
        end 
    elseif staggered_adoption && length(treatment_times) != length(treated_states)
        error("Er39: 'treatment_times' should be the same length as the 'treated_states'.")
    end 

    # Also need to make sure that start_times < treatment_times < end_times is true for each state
    for i in eachindex(treated_states)
        s = treated_states[i]
        treat_time = treatment_times[i]
        state_dates = data_copy[data_copy.state_71X9yTx .== s, :time_71X9yTx]
        earliest = minimum(state_dates)
        latest = maximum(state_dates)
        if !(earliest < treat_time)
            error("Er40 $s: For state $s, the earliest date ($earliest) is not strictly less than the treatment time ($treat_time).")
        end
        if !(treat_time <= latest)
            error("Er41 $s: For state $s, the treatment time ($treat_time) is greater than the last date ($latest).")
        end
    end

    # Ensure state column is a string (as are treated_states)
    data_copy.state_71X9yTx = string.(data_copy.state_71X9yTx)
    treated_states = string.(treated_states)
    check_states = unique(data_copy.state_71X9yTx)
    missing_states = setdiff(treated_states, check_states)
    if !isempty(missing_states)
        error("Er88: The states: $missing_states could not be found among the states in the data: $check_states")
    end

    # Once any date matching procedures are done, convert time_71X9yTx back to a string for processing in `categorical()`
    # Also keep a column vector copy as a date
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    if common_adoption
        data_copy.time_71X9yTx = ifelse.(data_copy.time_71X9yTx .>= treatment_times[1], "post", "pre")
        data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)
    elseif staggered_adoption
        data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)
    else
        error("Er42: A non-common or non-staggered adoption scenario was discovered!?")
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
                            error("Er43 $refcat $cov: Reference category '$refcat' not found in column '$cov'.")
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
                    error("Er44 $cov: Only detected one unique factor ($unique_categories) in factor variable $cov.")
                end 
            elseif cov_type <:Number
                data_copy[!, cov] = convert(Vector{Float64}, data_copy[!, cov])
                push!(covariates_to_include, cov)
            else
                error("Er45 $cov: column was found to be ($cov_type) neither of type Number, AbstractString, nor CategoricalValue!")
            end
        end
    end 

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
        error("Er46: 'ccc' must be set to one of: \"int\", \"time\", \"state\", \"add\", or \"hom\".")
    end 
    formula_str *= ")"
    formula_expr = Meta.parse(formula_str)
    formula = eval(formula_expr)

    # Call GC.gc() before running big FixedEffectsModels regression
    GC.gc()

    # Run FixedEffectsModels regression (stage 1)
    stage1 = reg(data_copy, formula; contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
                 save = false)

    # Recover lambdas
    state_time = split.(replace.(coefnames(stage1), "state_time: " => ""), "0IQR7q6Wei7Ejp4e")
    lambda_df = DataFrame(state = first.(state_time), time = last.(state_time))
    lambda_df.lambda = coef(stage1)
    cornercase = false
    if common_adoption && length(unique(lambda_df.state)) == 2
        lambda_df.stderr = stderror(stage1)

        cornercase = true
        treat_post_var = (lambda_df[lambda_df.state .== treated_states[1] .&& lambda_df.time .== "post", "stderr"][1])^2
        treat_pre_var = (lambda_df[lambda_df.state .== treated_states[1] .&& lambda_df.time .== "pre", "stderr"][1])^2
        control_post_var = (lambda_df[lambda_df.state .!= treated_states[1] .&& lambda_df.time .== "post", "stderr"][1])^2
        control_pre_var = (lambda_df[lambda_df.state .!= treated_states[1] .&& lambda_df.time .== "pre", "stderr"][1])^2

        # using fixedeffectsmodels means the off-diagonal entries of the covariance matrix will all be zero
        # so when covariates are included this estimate will not be accurate! fix later
        cornercase_se = sqrt(treat_post_var + treat_pre_var + control_post_var + control_pre_var)

        select!(lambda_df, Not(:stderr))
    end 

    # Compute diff for each treated state
    unique_states = unique(lambda_df.state)
    diff_df = DataFrame(state = String[], treated_time = Date[],
                        t = Date[], r1 = Date[], diff = Float64[],
                        treat = Int[], n = Int[])

    # Parse the treatment_times and the time column to dates for staggered adoption
    if staggered_adoption
        lambda_df.time = Date.(lambda_df.time)
        time_to_index = Dict(time => idx for (idx, time) in enumerate(all_times))
    end    

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
            for j in eachindex(diffs)
                diffs[j] = temp[j, "lambda"][1] - lambda_r1
                post = temp[j, "time"]
                pre = one_period_prior_treatment
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([post, pre])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            end 
            temp_df = DataFrame(state = treated_states[i], treated_time = treatment_times[i],
                                t = years, r1 = one_period_prior_treatment, diff = diffs, treat = 1,
                                n = n)
            diff_df = vcat(diff_df, temp_df) 
        end
        
        # Compute diffs for control states
        unique_diffs = unique(select(diff_df, :t, :r1, :treated_time))
        control_states = setdiff(unique_states, treated_states)     
        for i in eachindex(control_states)
            temp = lambda_df[(lambda_df.state .== control_states[i]), :]
            diffs = Vector{Float64}(undef, nrow(unique_diffs))
            n = Vector{Int}(undef, nrow(unique_diffs))
            for j in 1:nrow(unique_diffs)
                t = unique_diffs[j,"t"]
                r1 = unique_diffs[j,"r1"]
                diffs[j] = temp[temp.time .== t, "lambda"][1] - temp[temp.time .== r1, "lambda"][1]
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== control_states[i]))
            end 
            trtd_time = [all_times[findfirst(==(t), all_times) + 1] for t in unique_diffs.r1]
            temp_df = DataFrame(state = control_states[i], treated_time = trtd_time,
                                t = unique_diffs.t, r1 = unique_diffs.r1, diff = diffs, treat = 0,
                                n = n)
            diff_df = vcat(diff_df, temp_df)
        end

        # Compute diffs for randomization inference
        ri_diff_df = DataFrame(state = String[], treated_time = Date[],
                               t = Date[], r1 = Date[], diff = Float64[],
                               treat = Int[], n = Int[])
        for i in eachindex(treated_states)
            ri_diffs = unique(select(diff_df[(diff_df.state .!= treated_states[i]) .& (diff_df.treated_time .!= treatment_times[i]), :], :t, :r1))
            temp = lambda_df[(lambda_df.state .== treated_states[i]), :]
            diffs = Vector{Float64}(undef, nrow(ri_diffs))
            n = Vector{Int}(undef, nrow(ri_diffs))
            for j in 1:nrow(ri_diffs)
                t = unique_diffs[j,"t"]
                r1 = unique_diffs[j,"r1"]
                diffs[j] = temp[temp.time .== t, "lambda"][1] - temp[temp.time .== r1, "lambda"][1]
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            end
            trtd_time = [all_times[findfirst(==(t), all_times) + 1] for t in ri_diffs.r1]
            temp_df = DataFrame(state = treated_states[i], treated_time = trtd_time,
                                t = ri_diffs.t, r1 = ri_diffs.r1, diff = diffs, treat = -1,
                                n = n)
            ri_diff_df = vcat(ri_diff_df, temp_df)
        end
        
        # Run final regression to compute ATT based on weighting/aggregation method
        if agg == "cohort"
            # Run final regression -- weighted by treatment time
            unique_treatment_times = sort!(unique(diff_df.treated_time))
            results = DataFrame(treatment_time = Vector{Date}(undef, length(unique_treatment_times)),
                                att_cohort = Vector{Float64}(undef, length(unique_treatment_times)),
                                agg_att = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)),
                                nperm = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)),
                                se_att_cohort = Vector{Float64}(undef, length(unique_treatment_times)),
                                pval_att_cohort = Vector{Float64}(undef, length(unique_treatment_times)),
                                jknifese_att_cohort = Vector{Union{Missing, Float64}}(undef, length(unique_treatment_times)),
                                jknifepval_att_cohort = Vector{Union{Missing, Float64}}(undef, length(unique_treatment_times)),
                                ri_pval_att_cohort = Vector{Union{Missing, Float64}}(missing, length(unique_treatment_times)))
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
            end
            results = compute_weights(results, data_copy, "cohort", weighting, treated_states, treatment_times)
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
                                                      verbose, seed, data_copy, weighting)
            return results
        elseif agg == "simple"
            # Run final regression -- weighted by time and r - 1 groups
            results = DataFrame(r1 = Vector{Date}(undef, nrow(unique_diffs)),
                                time = Vector{Date}(undef, nrow(unique_diffs)),
                                gvar = Vector{Date}(undef, nrow(unique_diffs)),
                                att_gt = Vector{Float64}(undef, nrow(unique_diffs)),
                                agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)),
                                nperm = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)),
                                se_att_gt = Vector{Float64}(undef, nrow(unique_diffs)),
                                pval_att_gt = Vector{Float64}(undef, nrow(unique_diffs)),
                                jknifese_att_gt = Vector{Float64}(undef, nrow(unique_diffs)),
                                jknifepval_att_gt = Vector{Float64}(undef, nrow(unique_diffs)),
                                ri_pval_att_gt = Vector{Union{Missing, Float64}}(missing, nrow(unique_diffs)))
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
            end
            sort!(results, [:time])
            sort!(results, [:r1])
            results = compute_weights(results, data_copy, "simple", weighting, treated_states, treatment_times)
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
                                                 verbose, seed, data_copy, weighting)
            return results
        elseif agg == "state"
            # Run final regression -- weighted by state
            results = DataFrame(state = Vector{String}(undef, length(treated_states)),
                                att_s = Vector{Float64}(undef, length(treated_states)),
                                agg_att = Vector{Union{Missing, Float64}}(missing, length(treated_states)),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, length(treated_states)),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, length(treated_states)),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, length(treated_states)),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, length(treated_states)),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, length(treated_states)),
                                nperm = Vector{Union{Missing, Float64}}(missing, length(treated_states)),
                                se_att_s = Vector{Float64}(undef, length(treated_states)),
                                pval_att_s = Vector{Float64}(undef, length(treated_states)),
                                jknifese_att_s = Vector{Float64}(undef, length(treated_states)),
                                jknifepval_att_s = Vector{Float64}(undef, length(treated_states)),
                                ri_pval_att_s = Vector{Union{Missing, Float64}}(missing, length(treated_states)))
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
            end
            results.tuple_state = custom_sort_order.(results.state)
            sort!(results,[order(:tuple_state)])
            select!(results, Not([:tuple_state]))
            results = compute_weights(results, data_copy, "state", weighting, treated_states, treatment_times)
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
                                                 verbose, seed, data_copy, weighting)
            return results
        elseif agg == "none"  
            # Run final regrsesion -- each diff weighted equally
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
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "none",
                                                 verbose, seed, data_copy, weighting)
            return results
        elseif agg == "sgt"
            # Run final regrsesion -- each ATT_sgt weighted equally
            unique_sgt = unique(select(diff_df[diff_df.treat .== 1,:], :state, :t, :treated_time))
            #treated_states = unique(unique_sgt.state)
            results = DataFrame(state = Vector{String}(undef, nrow(unique_sgt)),
                                gvar = Vector{Date}(undef, nrow(unique_sgt)),
                                t = Vector{Date}(undef, nrow(unique_sgt)),
                                att_sgt = Vector{Float64}(undef, nrow(unique_sgt)),
                                agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)),
                                nperm = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)),
                                se_att_sgt = Vector{Float64}(undef, nrow(unique_sgt)),
                                pval_att_sgt = Vector{Float64}(undef, nrow(unique_sgt)),
                                jknifese_att_sgt = Vector{Float64}(undef, nrow(unique_sgt)),
                                jknifepval_att_sgt = Vector{Float64}(undef, nrow(unique_sgt)),
                                ri_pval_att_sgt = Vector{Union{Missing, Float64}}(missing, nrow(unique_sgt)))
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
            end 
            results = compute_weights(results, data_copy, "sgt", weighting, treated_states, treatment_times)
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
                                                 verbose, seed, data_copy, weighting)
        end
    elseif common_adoption
        states = unique(lambda_df.state)
        diff_df = DataFrame(state = states,
                            diff = Vector{Float64}(undef, length(states)),
                            treat = Vector{Float64}(undef, length(states)),
                            treated_time = fill(unique(treatment_times)[1], length(states)),
                            n = Vector{Float64}(undef, length(states)))
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
        end 
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
                                             verbose, seed, data_copy, weighting)

        return results
    end 



end 