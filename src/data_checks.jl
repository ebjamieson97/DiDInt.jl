
## The following functions have to do with data checks/validation/processing common to both didint_plot and didint_estimate
function hc_checks(hc)

    # Check hc args
    if hc isa Number 
        hc = round(hc)
        if !(hc in [0, 1, 2, 3, 4])
            error("'hc' must be one of 0, 1, 2, 3, or 4.")
        end
        hc = string("hc", hc)
    end
    hc = lowercase(replace(hc, r"\s" => ""))
    if !(hc in ["hc0", "hc1", "hc2", "hc3", "hc4"])
        error("'hc' must be one of $(join(["hc0", "hc1", "hc2", "hc3", "hc4"], ","))")
    end

    return hc
end

function get_sep_info(date::String)
    sep_positions = findall(c -> c in ['/', '-'], date)
    sep_types = unique(date[i] for i in sep_positions)
    return (sep_positions, sep_types)
end

function check_missing_vals(data_copy, state, time, covariates)

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

    return data_copy
end

function symbol_to_string(outcome, state, time, gvar)

    outcome = typeof(outcome) <: Symbol ? string(outcome) : outcome
    state = typeof(state) <: Symbol ? string(state) : state
    time = typeof(time) <: Symbol ? string(time) : time
    gvar = typeof(gvar) <: Symbol ? string(gvar) : gvar

    return outcome, state, time, gvar
end

function gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                  event = false, check = "estimate_didint")
    if isnothing(gvar)
        # treatment_times method being used
        if isnothing(treatment_times)
            error("If 'gvar' is not specified, then 'treatment_times' must be specified!")
        end
        
        # Check treated_states only when needed
        if isnothing(treated_states) && (event || check == "estimate_didint")
            error("If 'gvar' is not specified, then 'treated_states' must also be specified!")
        end
        
        return treated_states, treatment_times
    else
        # gvar method is being used
        if !isnothing(treatment_times) || !isnothing(treated_states)
            error("If 'gvar' is specified, do not specify 'treatment_times' and 'treated_states'!")
        end
        
        mask = .!ismissing.(data[!, gvar]) .&& data[!, gvar] .!= 0 .&& data[!, gvar] .!= "0" .&& data[!, gvar] .!= "0.0"
        filtered = unique(data[mask, [state, gvar]])
        treated_states = filtered[!, state]
        treatment_times = filtered[!, gvar]
        
        return treated_states, treatment_times
    end
end

function start_end_date_checks(date_value, date_name, date_format)

    if isnothing(date_value)
        return nothing
    end
    
    # If it's a number or string, convert to Date, otherwise the function inputs
    # already require to be a date, so just leave it as is
    if eltype(date_value) <: Number
        if length(string(date_value)) != 4
            error("If '$date_name' is a number it must be 4 digits long.")
        end
        return Date(date_value)
    elseif date_value isa String
        if isnothing(date_format)
            error("If '$date_name' is entered as a string, please specify 'date_format'.")
        end
        return parse_string_to_date_didint(date_value, date_format)
    else
        return date_value
    end
end

function validate_data(data_copy, outcome, state, time, covariates, treatment_times, date_format; 
                                   warn_missing_covariates = true)
    
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
    elseif warn_missing_covariates
        @warn "No covariates specified!"
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
                error("If 'treatment_times' are entered as numbers, they must all be 4 digits long in 'yyyy' date_format. 'treatment_times' found are $(join(treatment_times, ", "))")
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
        data_copy[!, time] = [ismissing(t) ? missing : round(Int, t) for t in data_copy[!, time]]
        time_column_numeric = true
        unique_time_lengths = unique(length.(string.(skipmissing(data_copy[!, time]))))
        if isempty(unique_time_lengths) || length(unique_time_lengths) > 1 || unique_time_lengths[1] != 4
            error("The 'time' column was found to be numeric but consisting of values of ambiguous date formatting (i.e. not consistent 4 digit entries.\nFound entries of the following lengths $(join(unique_time_lengths, ", "))")
        end 
    else
        time_column_numeric = false
    end 
    
    if xor(time_column_numeric, treatment_times_numeric)
        error("If 'time' column is numeric or 'treatment_times' is numeric, then both must be numeric.")
    end
    
    return data_copy, treatment_times, covariates
end

function validate_state_types(treated_states, data_copy, state)
    
    # Validate that treated_states and the state column have compatible types.
    # Both must be either numeric or string types.
    
    treated_states_type = Base.nonmissingtype(eltype(treated_states))
    state_column_type = Base.nonmissingtype(eltype(data_copy[!, state]))
    
    if !((treated_states_type <: Number && state_column_type <: Number) || 
         (treated_states_type <: AbstractString && state_column_type <: AbstractString))
        error("'treated_states' and the 'state' column ($state) must both be numerical or both be strings.\n" *
              "Instead, found 'treated_states' $treated_states_type and '$state' $state_column_type.")
    end
end

function validate_and_convert_dates(data_copy, time, treatment_times, date_format, freq, start_date, end_date)

    # Validate date formats in time column and treatment_times, convert to Date objects,
    # validate freq argument, and filter data by date range.
    
    # Check that time column entries are all in the same date format
    nonmissing_time_type = Base.nonmissingtype(eltype(data_copy[!, time]))
    if nonmissing_time_type <: AbstractString || nonmissing_time_type <: Number
        dates_str = string.(data_copy[!, time])
        if length(unique(length.(dates_str))) != 1
            error("Dates in the 'time' column are not all the same length!")
        end
        if unique(length.(dates_str))[1] != 4 && nonmissing_time_type <: Number
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
    if nonmissing_time_type <: AbstractString && eltype(treatment_times) <: AbstractString
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
    
    return data_copy, treatment_times, start_date, end_date, all_times, freq
end

function validate_treated_states(treated_states, treatment_times, data_copy)
    missing_states = setdiff(treated_states, data_copy.state_71X9yTx)
    if !isempty(missing_states)
        @warn "The following 'treated_states' could not be found in the data $(join(missing_states, ", ")).\nOnly found the following states $(join(unique(data_copy.state_71X9yTx), ", "))"
        keep_mask = treated_states .∉ Ref(missing_states)
        treated_states = treated_states[keep_mask]
        treatment_times = treatment_times[keep_mask]
        if length(treated_states) == 0 || length(treatment_times) == 0
            error("No valid 'treated_states' were found in the data.")
        end
    end
    return treated_states, treatment_times
end

function validate_string_treated_states(data_copy; treated_states = nothing, event = false)

    data_copy.state_71X9yTx = string.(data_copy.state_71X9yTx)
    if event == true
        treated_states = string.(treated_states)
        check_states = unique(data_copy.state_71X9yTx)
        missing_states = setdiff(treated_states, check_states)
        if !isempty(missing_states)
            error("The states $missing_states could not be found among the states in the data $check_states")
        end
    end
    return data_copy, treated_states

end

function process_covariates(covariates, data_copy, ref)
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
                            error("Reference category '$refcat' not found in column '$cov'.")
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

    return data_copy, covariates_to_include

end