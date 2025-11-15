
### These two functions pass string and date objects back and forth
function parse_string_to_date_didint(date::Union{Vector{<:AbstractString}, AbstractString},
                                     date_format::AbstractString)
    
    # Define mapping from months to numbers and possible date formats
    month_map = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04",
                     "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08",
                     "sep" => "09", "oct" => "10", "nov" => "11", "dec" => "12")
    possible_formats = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm",
                        "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy",
                        "mm-dd-yyyy", "mmddyyyy", "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy",
                        "ddmonyyyy", "yyyym00"]

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
    elseif date_format in possible_formats  
        # Other formats are handled natively by Dates      
        output = Date(date, lowercase(date_format))
    else
        error("Please specify a date_format listed here $possible_formats.")
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

### These functions are all apart of the date matching / period sorting suite of functions:
function perform_date_matching(data_copy, all_times, freq, freq_multiplier, start_date, end_date, 
                               treatment_times, treated_states)

    # Meta function for the date matching procedure to align dates to the constructed grid
    
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
    
    one_past = end_date + period
    matched = [match_date(t, match_to_these_dates, treatment_times) for t in data_copy.time_71X9yTx]
    data_copy.time_71X9yTx = matched
    treated_states = treated_states[(treatment_times .< one_past) .&& (treatment_times .> start_date)]
    treatment_times = treatment_times[(treatment_times .< one_past) .&& (treatment_times .> start_date)]
    treatment_times = [match_treatment_time(t, match_to_these_dates) for t in treatment_times]
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
        # If it is 1 day off from being cleanly divided into years, assume leap year
        if days > 365
            if days % 365 == 1 
                years = floor(days / 365)
                return Year(years)
            end
            # Days is > 1461, it is bound to have one leap year
            if days > 1460
                leap_count = days ÷ (365*4)
                if days % 365 == leap_count
                    years = floor(days / 365)
                    return Year(years)
                end 
            end
        end
        years = ceil(Int, days / 365)
        return Year(years)
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

    # For edge cases in didint_plot where date_format wasn't specified, but period length has between
    if typeof(period) <: Year
        return "yyyy"
    else
        return "yyyy/mm/dd"
    end
end

# This one is not currently used
function get_freq_approx(times)

    # This function finds the maximum period_length such that
    # period_{t-1} + period_length <= period_{t+1} for all periods in the data

    # find max distances between dates
    max_dist = maximum(diff(times))

    # Test each period type to find the largest that fits
    # Start with largest periods and work down

    max_years = ceil(max_dist / Day(365))
    max_months = ceil(max_dist / Day(28))
    max_weeks = ceil(max_dist / Day(7))
    max_days = ceil(max_dist / Day(1))

    # Test years
    for n in max_years:-1:1
        period = Year(n)
        if all(times[i] + period <= times[i+1] for i in 1:length(times)-1)
            return period
        end
    end
    
    # Test months
    for n in max_months:-1:1
        period = Month(n)
        if all(times[i] + period <= times[i+1] for i in 1:length(times)-1)
            return period
        end
    end
    
    # Test weeks
    for n in max_weeks:-1:1
        period = Week(n)
        if all(times[i] + period <= times[i+1] for i in 1:length(times)-1)
            return period
        end
    end
    
    # Test days
    for n in max_days:-1:1
        period = Day(n)
        if all(times[i] + period <= times[i+1] for i in 1:length(times)-1)
            return period
        end
    end


end