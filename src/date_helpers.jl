###########################################################################################
# date_helpers.jl: parsing dates and matching them to periods                             #
#                                                                                         #
#   Section 1: Converting between strings and dates                                       #
#   Section 2: The period grid and the matching of dates to periods                       #
###########################################################################################


###########################################################################################
# Section 1: Converting between strings and dates                                         #
###########################################################################################

# Date formats that the 'date_format' argument accepts. The last two are how Stata prints dates
# as strings, for example 25dec1991 and 1991m12.
const DATE_FORMATS = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm",
                      "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy",
                      "mm-dd-yyyy", "mmddyyyy", "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy",
                      "ddmonyyyy", "yyyym00"]

function parse_string_to_date_didint(date::Union{Vector{<:AbstractString}, AbstractString},
                                     date_format::AbstractString)
    
    # Mapping from month abbreviations to numbers
    month_map = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04",
                     "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08",
                     "sep" => "09", "oct" => "10", "nov" => "11", "dec" => "12")

    # Convert strings to date objects depending on the specified date_format
    if date_format == "ddmonyyyy"
        # The first case is dealing with formats such as 25dec1991 
        # which is what dates sometimes look like when converting to strings in Stata
        date = lowercase(date)
        day = date[1:2]
        month_str = date[3:5]
        year = date[6:end]
        month = month_map[month_str]
        output = Date("$day/$month/$year", "dd/mm/yyyy")
    elseif date_format == "yyyym00"
        # This is for another common Stata formatting of dates when converted to strings
        date = lowercase(date)
        info = split(date, "m")        
        output = Date("$(info[2])/$(info[1])", "mm/yyyy")
    elseif date_format in DATE_FORMATS
        # Other formats are handled natively by Dates      
        output = Date(date, lowercase(date_format))
    else
        error("Please specify a date_format listed here $DATE_FORMATS.")
    end 

    return output 
end

function parse_date_to_string_didint(date, date_format::AbstractString)
    
    # This function is basically a wrapper Dates.format()
    # except this adds a bit more functionality so that it can 
    # take date strings in Stata formats (e.g. 25dec2020 or 2020m12) and return those as strings
    if date_format == "ddmonyyyy"
        month_dict = Dict("01" => "jan", "02" => "feb", "03" => "mar", "04" => "apr", "05" => "may", 
        "06" => "jun", "07" => "jul", "08" => "aug", "09" => "sep", "10" => "oct", 
        "11" => "nov", "12" => "dec")
        vectorized_date_object = split(string(date), "-")
        return "$(vectorized_date_object[3])$(month_dict[vectorized_date_object[2]])$(vectorized_date_object[1])"        
    elseif date_format == "yyyym00"
        month_dict = Dict("01" => "m1", "02" => "m2", "03" => "m3", "04" => "m4", "05" => "m5", 
        "06" => "m6", "07" => "m7", "08" => "m8", "09" => "m9", "10" => "m10", 
        "11" => "m11", "12" => "m12")
        vectorized_date_object = split(string(date), "-")
        return "$(vectorized_date_object[1])$(month_dict[vectorized_date_object[2]])"
    else
        return Dates.format(date, date_format)
    end 
end


###########################################################################################
# Section 2: The period grid and the matching of dates to periods                         #
#                                                                                         #
# The analysis runs on a grid of evenly spaced periods from start_date to end_date. Each  #
# observation is matched to the latest grid date at or before its own date, except that   #
# an observation at or after a treatment time that falls between two grid dates is        #
# matched to the later one. Treatment times are matched to the first grid date at or      #
# after them.                                                                             #
###########################################################################################

function perform_date_matching(data_copy, all_times, freq, freq_multiplier, start_date, end_date,
                               treatment_times, treated_states; drop_already_treated = false)

    # Build the period grid (of the given `freq`, or of the longest gap between dates in the data)
    # and match the dates of the observations and the treatment times to it.
    # A state whose treatment time is not within the sample is not a treated state of the
    # analysis and is removed from `treated_states`, with a warning:
    #   - treated at or before the first date: it is treated throughout, so it has no
    #     pre-treatment period and is not a valid control either. With `drop_already_treated`
    #     its observations are removed from the data.
    #   - treated after the last date: it is untreated throughout and serves as a control.

    # Do date matching procedure
    if !isnothing(freq)
        max_dist = Dates.value(maximum(diff(all_times)))
        period = parse_freq(string(freq_multiplier)*" "*freq)
        match_to_these_dates = collect(start_date:period:end_date)
        max_dist_grid = Dates.value(maximum(diff(match_to_these_dates)))
        if max_dist > max_dist_grid
            period_str = string(period)
            @warn "The specified period length $period_str is less than the maximum observed period length ($max_dist days)."
        end
    else
        period = get_max_period(all_times)
        match_to_these_dates = collect(start_date:period:end_date)
    end 
    
    # Treatment times outside the sample
    already_treated = treatment_times .<= start_date
    treated_later = treatment_times .> end_date
    if !isnothing(treated_states)
        if any(already_treated)
            states = join(treated_states[already_treated], ", ")
            if drop_already_treated
                @warn "The following states are treated at or before the first date of the data used ($start_date), so they have no pre-treatment period and cannot serve as controls. Their observations are dropped: $states"
                data_copy = data_copy[.!in.(data_copy.state_71X9yTx, Ref(treated_states[already_treated])), :]
            else
                @warn "The following states are treated at or before the first date of the data used ($start_date) and are not counted as treated states: $states"
            end
        end
        if any(treated_later)
            @warn "The following states are treated after the last date of the data used ($end_date), so they are untreated throughout and are used as controls: $(join(treated_states[treated_later], ", "))"
        end
    end
    in_sample = .!already_treated .& .!treated_later

    # Match each distinct date once, then every observation through the lookup
    matched_date = Dict(t => match_date(t, match_to_these_dates, treatment_times[in_sample]) for t in unique(data_copy.time_71X9yTx))
    data_copy.time_71X9yTx = [matched_date[t] for t in data_copy.time_71X9yTx]
    treated_states = isnothing(treated_states) ? nothing : treated_states[in_sample]
    treatment_times = [match_treatment_time(t, match_to_these_dates) for t in treatment_times[in_sample]]
    time_to_index = Dict(time => idx for (idx, time) in enumerate(match_to_these_dates))
    
    return data_copy, treatment_times, treated_states, match_to_these_dates, time_to_index, period
end

function parse_freq(period_str::AbstractString)
    parts = split(period_str)
    value = parse(Int, parts[1])
    period_type = parts[2]
    
    if period_type in ["week", "weeks", "weekly"]
        return Week(value)
    elseif period_type in ["day", "days", "daily"]
        return Day(value)
    elseif period_type in ["month", "months", "monthly"]
        return Month(value)
    elseif period_type in ["year", "years", "yearly"]
        return Year(value)
    else
        error("Unsupported period type $period_type, try day(s), week(s), month(s), or year(s).")
    end
end

function match_date(t::Date, periods::Vector{Date}, treatment_times::Vector{Date})
    # Find the index of the last period that is ≤ t.
    i = nothing
    try
        i = findlast(p -> p <= t, periods)
    catch e
        i = nothing
    end
    if i === nothing
        return periods[1]
    end
    if i == length(periods)
        return periods[end]
    end
    current_period = periods[i]
    next_period = periods[i+1]
    candidate_T = filter(x -> current_period < x < next_period, treatment_times)
    if isempty(candidate_T)
        return current_period
    else
        T = minimum(candidate_T)
        return (t < T) ? current_period : next_period
    end
end

function match_treatment_time(t::Date, periods::Vector{Date})
    i = searchsortedfirst(periods, t)
    if i > length(periods)
        return periods[end]
    else
        return periods[i]
    end
end

function get_max_period(times)

    # This function simply returns the maximum period length found in the data

    # find max distances between dates
    max_dist = maximum(diff(times))

    # Convert to days for comparison
    days = Dates.value(max_dist)
    
    # Return a period that is AT LEAST as large as the max gap
    # This ensures no "empty" periods in the grid 
    
    if days >= 365
        # A gap of k whole years is 365k days plus the leap days that fall in it, so it is
        # within a day of 365.2425k. Such a gap is k years; any other gap is rounded up.
        whole_years = round(Int, days / 365.2425)
        if abs(days - 365.2425 * whole_years) < 1
            return Year(whole_years)
        end
        return Year(ceil(Int, days / 365))
    elseif days >= 28
        months = ceil(Int, days / 31)
        return Month(months)
    elseif days >= 7
        weeks = ceil(Int, days / 7)
        return Week(weeks)
    else
        return Day(days)
    end
end

function assume_date_format(period)

    # Date format for the output when none was specified: years for yearly periods, full dates otherwise
    if typeof(period) <: Year
        return "yyyy"
    else
        return "yyyy/mm/dd"
    end
end