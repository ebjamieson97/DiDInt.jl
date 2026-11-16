###########################################################################################
# data_checks.jl: argument checks and data cleaning shared by didint() and didint_plot()  #
#                                                                                         #
#   Section 1: Checks of single arguments                                                 #
#   Section 2: Treatment assignment (gvar, treated_states, treatment_times, policies)     #
#   Section 3: Columns, missing values and states                                         #
#   Section 4: Dates                                                                      #
#   Section 5: Covariates                                                                 #
#   Section 6: Output for the Stata and R wrappers                                        #
###########################################################################################


###########################################################################################
# Section 1: Checks of single arguments                                                   #
###########################################################################################

function data_type_check(data)

    if !(data isa DataFrame)
        if !((string(typeof(data)) == "Main.RConnector.RDataFrame") || (string(typeof(data)) == "RConnector.RDataFrame"))
	    error("'data' needs to be of type DataFrame (or RConnector.RDataFrame if running from R via JuliaConnectoR)")
        end
    end

end

function hc_checks(hc)

    # The HC covariance estimator as one of "hc0", ..., "hc4", from a number or a string
    if hc isa Number
        hc = round(Int, hc)
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

function check_option(value::AbstractString, name::AbstractString, options::Vector{String})

    # A string argument in lower case, after checking that it is one of the options
    value = lowercase(value)
    if !(value in options)
        error("'$name' must be one of $(options)")
    end
    return value
end

function check_positive_integer(value::Number, name::AbstractString)

    # A numeric argument rounded to an integer, which must be at least 1
    value = Int(round(value))
    if value < 1
        error("'$name' must be a positive integer.")
    end
    return value
end

function check_ri_kernel_cutoff(ri_kernel_cutoff)

    # The kernel smoothing cutoff as a non-negative integer; nothing or false turn kernel
    # smoothing off (the same as 0)
    if ri_kernel_cutoff === true
        error("'ri_kernel_cutoff' must be a non-negative integer, or nothing or false to turn kernel smoothing off.")
    end
    ri_kernel_cutoff = (isnothing(ri_kernel_cutoff) || ri_kernel_cutoff === false) ? 0 : Int(round(ri_kernel_cutoff))
    if ri_kernel_cutoff < 0
        error("'ri_kernel_cutoff' must be a non-negative integer, or nothing or false to turn kernel smoothing off.")
    end
    return ri_kernel_cutoff
end

function check_date_format(date_format)

    # The date format in lower case (as the date parsing and printing expect it; in upper case
    # "MM" would mean minutes), after checking that it is a supported format or at least only
    # made of years, months, days and separators. Nothing if none was given.
    isnothing(date_format) && return nothing
    date_format = lowercase(date_format)
    if !(date_format in DATE_FORMATS) && !occursin(r"^[ymd/\-. ]+$", date_format)
        error("'date_format' was not recognized. Try one of $DATE_FORMATS.")
    end
    return date_format
end

function check_treated_states(treated_states, treatment_times)

    # `treated_states` and `treatment_times` as vectors of the same length, with each state
    # listed once. A single treatment time applies to every treated state: common adoption is
    # the one-cohort case of the same estimator. `treated_states` may be nothing (didint_plot
    # does not always need it).
    treatment_times = !(treatment_times isa AbstractVector) ? [treatment_times] : treatment_times
    if isempty(treatment_times)
        error("'treatment_times' must have at least one entry.")
    end
    isnothing(treated_states) && return treated_states, treatment_times

    treated_states = !(treated_states isa AbstractVector) ? [treated_states] : treated_states
    if length(treatment_times) == 1 && length(treated_states) > 1
        treatment_times = fill(treatment_times[1], length(treated_states))
    end
    if length(treatment_times) != length(treated_states)
        error("'treatment_times' should either be the same length as the 'treated_states' vector or of length 1.")
    end
    if !allunique(treated_states)
        repeated = unique(s for s in treated_states if count(isequal(s), treated_states) > 1)
        error("Each treated state can only have one treatment time, but the following states are listed more than once in 'treated_states' (or have more than one value in the 'gvar' column): $(join(repeated, ", ")).")
    end
    return treated_states, treatment_times
end

function symbol_to_string(outcome, state, time, gvar)

    outcome = typeof(outcome) <: Symbol ? string(outcome) : outcome
    state = typeof(state) <: Symbol ? string(state) : state
    time = typeof(time) <: Symbol ? string(time) : time
    gvar = typeof(gvar) <: Symbol ? string(gvar) : gvar

    return outcome, state, time, gvar
end

function start_end_date_checks(date_value, date_name, date_format)

    if isnothing(date_value)
        return nothing
    end

    # If it's a number or string, convert to Date, otherwise the function inputs
    # already require to be a date, so just leave it as is
    if eltype(date_value) <: Number
        date_value = round(Int, date_value)
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


###########################################################################################
# Section 2: Treatment assignment (gvar, treated_states, treatment_times, policies)       #
###########################################################################################

function gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                  event = false, check = "estimate_didint")

    # The treated states and their treatment times, from the `gvar` column or from the
    # `treated_states` and `treatment_times` arguments (one or the other)
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
        # gvar method is being used: a state is treated if its gvar is not missing and not 0
        if !isnothing(treatment_times) || !isnothing(treated_states)
            error("If 'gvar' is specified, do not specify 'treatment_times' and 'treated_states'!")
        end

        mask = .!ismissing.(data[!, state]) .&& .!ismissing.(data[!, gvar]) .&& data[!, gvar] .!= 0 .&& data[!, gvar] .!= "0" .&& data[!, gvar] .!= "0.0"
        filtered = unique(data[mask, [state, gvar]])
        treated_states = filtered[!, state]
        treatment_times = filtered[!, gvar]

        return treated_states, treatment_times
    end
end

# Policy label of every treated state when neither `pvar` nor `treated_policies` is given
const DEFAULT_POLICY = "1"

function normalize_policy_label(value)

    # A policy label as a string, or nothing for "no policy": missing, 0, "0" or ""
    ismissing(value) && return nothing
    if value isa Real
        value == 0 && return nothing
        return isinteger(value) ? string(Int(value)) : string(value)
    end
    label = strip(string(value))
    return label in ("", "0", "0.0") ? nothing : String(label)
end

function policy_assignments(pvar, treated_policies, treated_states, gvar, data, state)

    # Policy of each treated state, keyed by the state as a string. Mirrors gvar /
    # treatment_times: `pvar` names a column holding each state's policy, `treated_policies`
    # gives one policy per element of `treated_states`; use one or the other. Without either,
    # every treated state has the same policy.
    states = treated_states isa AbstractVector ? treated_states : [treated_states]
    if !isnothing(pvar) && !isnothing(treated_policies)
        error("Specify either 'pvar' or 'treated_policies', not both!")
    end
    if isnothing(pvar) && isnothing(treated_policies)
        return Dict(string(s) => DEFAULT_POLICY for s in states)
    end

    if !isnothing(treated_policies)
        if !isnothing(gvar)
            error("'treated_policies' is matched to 'treated_states' by position, so it cannot be used with 'gvar'. Use 'pvar' instead.")
        end
        policies = treated_policies isa AbstractVector ? treated_policies : [treated_policies]
        if length(policies) == 1 && length(states) > 1
            policies = fill(policies[1], length(states))
        end
        if length(policies) != length(states)
            error("'treated_policies' should either be the same length as the 'treated_states' vector or of length 1.")
        end
        labels = normalize_policy_label.(policies)
        if any(isnothing, labels)
            error("Every treated state needs a policy: 'treated_policies' cannot contain missing, 0, \"0\" or \"\".")
        end
        return Dict(string(s) => label for (s, label) in zip(states, labels))
    end

    # pvar: the policy must be constant within a state and present for every treated state
    pvar = string(pvar)
    if !(pvar in names(data))
        error("The 'pvar' column $pvar could not be found in the data.")
    end
    assignments = Dict{String, String}()
    for s in states
        labels = unique(normalize_policy_label.(data[isequal.(data[!, state], s), pvar]))
        labels = filter(!isnothing, labels)
        if isempty(labels)
            error("No policy found in the 'pvar' column $pvar for treated state $s.")
        elseif length(labels) > 1
            error("The 'pvar' column $pvar is not constant within state $s: found $(join(labels, ", ")).")
        end
        assignments[string(s)] = labels[1]
    end
    return assignments
end

function validate_treated_states(treated_states, treatment_times, data_copy)

    # Drop, with a warning, the treated states that are not in the data
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

function check_treatment_windows(treated_states, treatment_times, data_copy)

    # Every treated state needs an observation before its treatment time (for the base period)
    # and one at or after it
    for (s, treat_time) in zip(treated_states, treatment_times)
        state_dates = data_copy[data_copy.state_71X9yTx .== s, :time_71X9yTx]
        earliest = minimum(state_dates)
        latest = maximum(state_dates)
        if !(earliest < treat_time)
            error("For state $s, the earliest date for which there is non-missing value ($earliest) is not strictly less than the treatment time ($treat_time).")
        end
        if !(treat_time <= latest)
            error("For state $s, the treatment time ($treat_time) is greater than the last date ($latest) for which there is a non-missing observation.")
        end
    end
end


###########################################################################################
# Section 3: Columns, missing values and states                                           #
###########################################################################################

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
            error("The following covariates could not be found in the data $(join(missing_cov, ", "))")
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

    # A numeric time column must hold 4 digit years
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

# A value that counts as missing: missing, nothing or NaN
is_missing_value(x) = x === missing || x === nothing || (x isa AbstractFloat && isnan(x))

function drop_missing_rows(data_copy, column, label)

    # Drop the rows where `column` is missing, with a warning that names the column as `label`
    is_missing = is_missing_value.(data_copy[!, column])
    if any(is_missing)
        @warn "Found missing values in the $label column. Dropping those rows."
        data_copy = data_copy[.!is_missing, :]
    end
    return data_copy
end

function check_missing_vals(data_copy, state, time, covariates)

    # Drop the rows with a missing outcome, state, time or covariate
    data_copy = drop_missing_rows(data_copy, "outcome_71X9yTx", "'outcome'")
    data_copy = drop_missing_rows(data_copy, state, "'state'")
    data_copy = drop_missing_rows(data_copy, time, "'time'")
    if !(isnothing(covariates))
        for cov in covariates
            data_copy = drop_missing_rows(data_copy, cov, "'$cov'")
        end
    end

    return data_copy
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

function validate_string_treated_states(data_copy; treated_states = nothing)

    # Convert the states (and the treated states) to strings
    data_copy.state_71X9yTx = string.(data_copy.state_71X9yTx)
    if isnothing(treated_states)
        return data_copy, treated_states
    end
    treated_states = string.(treated_states)
    check_states = unique(data_copy.state_71X9yTx)
    missing_states = setdiff(treated_states, check_states)
    if !isempty(missing_states)
        error("The states $missing_states could not be found among the states in the data $check_states")
    end
    return data_copy, treated_states

end


###########################################################################################
# Section 4: Dates                                                                        #
###########################################################################################

function get_sep_info(date::String)

    # Positions and types of the separators (/ or -) in a date string
    sep_positions = findall(c -> c in ['/', '-'], date)
    sep_types = unique(date[i] for i in sep_positions)
    return (sep_positions, sep_types)
end

function check_date_strings(dates, name)

    # Check that every date string has the separators of the first one, in the same positions
    ref_positions, ref_sep_types = get_sep_info(dates[1])
    if length(ref_sep_types) > 1
        error("First date in $name uses mixed separators: - and /.")
    end
    for (i, date) in enumerate(dates)
        sep_positions, sep_types = get_sep_info(date)
        if sep_positions != ref_positions
            error("$i Separator positions differs from the first date in $name in date entry $i $date")
        end
        if length(sep_types) > 1
            error("$i Date found in $name (entry $i) which uses multiple separator types $date")
        end
        if sep_types != ref_sep_types
            error("$i Date found in $name (entry $i $date) which uses different separator types from first date.")
        end
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
        check_date_strings(dates_str, "'time' column")
    end

    # Check that treatment_times, if strings, are entered in the same date format
    if eltype(treatment_times) <: AbstractString
        check_date_strings(treatment_times, "'treatment_times'")
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

    # Make sure the time column is a Date object
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
    if Base.nonmissingtype(eltype(treatment_times)) <: Date
        treatment_times = convert(Vector{Date}, treatment_times)
    else
        if isnothing(date_format)
            error("If 'treatment_times' are entered as strings, please specify 'date_format'.")
        end
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

    # Keep the data between start_date and end_date (by default, the first and last dates in the data)
    all_times = sort(unique(data_copy.time_71X9yTx))
    if isnothing(start_date)
        start_date = minimum(all_times)
    end
    if isnothing(end_date)
        end_date = maximum(all_times)
    end
    data_copy = data_copy[(data_copy.time_71X9yTx .>= start_date) .& (data_copy.time_71X9yTx .<= end_date), :]

    if nrow(data_copy) == 0
        error("After filtering by 'start_date' and 'end_date' found only 0 rows of data!")
    end

    return data_copy, treatment_times, start_date, end_date, all_times, freq
end


###########################################################################################
# Section 5: Covariates                                                                   #
###########################################################################################

function process_covariates(covariates, data_copy, ref)

    # Convert factor covariates into numeric dummy variable columns, one per category other
    # than the reference category (`ref`, or the first category seen), and numeric covariates
    # to Float64. Returns the data and the names of the covariate columns to use.
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
                        # One dummy for the non-reference category, in the column it is registered under
                        newcol = Symbol("$(cov)_not_$(refcat)")
                        data_copy[!, newcol] = (data_copy[!, cov] .!= refcat) .|> Int
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


###########################################################################################
# Section 6: Output for the Stata and R wrappers                                          #
###########################################################################################

function init_wrapper_check(wrapper)

    if !isnothing(wrapper)
        wrapper = lowercase(replace(wrapper, r"\s" => ""))
        wrapper_options = ["julia", "stata", "r"]
        if !(wrapper in wrapper_options)
            error("'wrapper' must be one of 'julia', 'stata', or 'r'.")
        end
    end
    return wrapper

end

function wrapper_check(results, wrapper)

    # Dates are returned to R as strings
    if isnothing(wrapper) || wrapper == "julia"
        return results
    elseif wrapper == "r"
	    for col in names(results)
            if eltype(results[!, col]) <: Date
           	    results[!, col] = string.(results[!, col])
       	    end
        end
    end
    return results

end
