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
        error("Please specify a date_format listed here: $possible_formats. Format 'ddmonyyyy' should look like '25dec2020' and format yyyym00 should look like '2020m12'.")
    end 

    return output 
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
        throw(ArgumentError("Er11: Unsupported period type: $period_type, try day(s), week(s), month(s), or year(s)."))
    end
end

function compute_jknife_se(X::Matrix{<:Number}, Y::Vector{<:Number}, original_att::Number)
    n = length(Y)
    jknife_beta = Vector{Float64}(undef, length(Y))
    ncolx = size(X,2)
    for i in eachindex(Y)
        idx = [1:i-1; i+1:size(X, 1)]
        X_sub = X[idx, :]
        Y_sub = Y[idx]
        jknife_beta[i] = (X_sub \ Y_sub)[ncolx]
    end 
    jknife_se = sqrt(sum((jknife_beta .- original_att).^2) * ((n - 1) / n))
    return jknife_se
end 

function final_regression_results(X::Matrix{<:Number}, Y::Vector{<:Number})
    beta_hat = nothing
    try
        beta_hat = (X \ Y) #
    catch e
        @warn "Direct solve failed, using pseudoinverse: $e"
        beta_hat = pinv(X' * X) * X' * Y
    end 
    ncolx = size(X, 2)
    resid = Y - X * beta_hat
    sigma_sq = resid .^ 2
    omega = Diagonal(sigma_sq)
    beta_hat_cov = nothing
    try
        beta_hat_cov = inv(X' * X) * (X' * omega * X) * inv(X' * X)
    catch e
        @warn "Direct solve failed, using pseudoinverse: $e"
        beta_hat_cov = pinv(X' * X) * (X' * omega * X) * pinv(X' * X)
    end 
    beta_hat_var = diag(beta_hat_cov)
    beta_hat_se = sqrt(beta_hat_var[ncolx]) #
    dof = length(Y) - ncolx
    pval_att = dof > 0 ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se))) : missing #
    beta_hat_se_jknife = compute_jknife_se(X, Y, beta_hat[ncolx]) #
    pval_att_jknife = dof > 0 ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se_jknife))) : missing
    result_dict = Dict("beta_hat" => beta_hat[ncolx], "beta_hat_se" => beta_hat_se, "pval_att" => pval_att,
                       "beta_hat_se_jknife" => beta_hat_se_jknife, "pval_att_jknife" => pval_att_jknife)
    return result_dict
end 

function randomization_inference_didint(diff_df::DataFrame, agg::AbstractString,
                                        original_att::Number, nperm::Integer, control_states::AbstractVector,
                                        treated_states::AbstractVector, verbose::Bool)

    # Check total number of possible permutations
    n = length(control_states) + length(treated_states)
    k = length(treated_states)
    n_unique_perms = binomial(n, k)
    if nperm > n_unique_perms
        @warn "'nperm' was set to $nperm but only $n_unique_perms unique permutations exist. \n 
Setting 'nperm' to $n_unique_perms."
        nperm = n_unique_perms
    end 
    if nperm < 500
        @warn "'nperm' is less than 500!"
    end 

    # Create df which keeps track of treatment_time randomizations and associated states
    init = copy(diff_df)
    init = unique(select(init[init.treat .!= -1,:], :state, :treated_time, :treat ))
    init.treated_time[init.treat .== 0] .= missing
    init = unique(init)
    treated_time = init.treated_time
    seen = Set{String}()
    key = join(treated_time, "")
    push!(seen, key)
    i = 1
    while i < nperm
        new_perm = shuffle(treated_time)
        key = join(new_perm, "")
        if !(key in seen)
            # Add a new column with the randomized permutation
            init[!, Symbol("treated_time_randomized_$i")] = new_perm
            push!(seen, key)
            i += 1
        end
    end
    
    # Create preallocation vector to store agg_ATT from each randomization
    ri_att = Vector{Float64}(undef, nperm - 1)

    # Make sure unique_diffs is defined, no need to repeatedly redfine during the loop
    if agg == "simple"
        unique_diffs = unique(select(diff_df, :t, :r1))
    end

    # Loop thru all the treated_time randomizations
    for j in 1:(nperm - 1)
        # Construct the ri_df for this iteration
        colname = Symbol("treated_time_randomized_" * string(j))
        mask_treated = ismissing.(init[!, colname])
        new_treated = init[.!mask_treated, colname]
        new_treated_states = init[.!mask_treated, :state]
        new_control_states = init[mask_treated, :state]
        ri_df = filter(row -> row.state in new_control_states, diff_df)
        ri_df.treat .= 0
        for i in eachindex(new_treated)
            state = new_treated_states[i]
            treatment_time = new_treated[i]
            rows_to_add = filter(row -> row.state == state && row.treated_time == treatment_time, diff_df)
            ri_df = vcat(ri_df, rows_to_add)
            ri_mask = (ri_df.state .== state) .& (ri_df.treated_time .== treatment_time)
            ri_df.treat[ri_mask] .= 1
        end

        # Compute agg_att for this ri_df conditional on aggregation method
        if agg == "cohort"
            unique_treatment_times = unique(ri_df.treated_time)
            att_cohort = Vector{Float64}(undef, length(unique_treatment_times))
            for i in eachindex(unique_treatment_times)
                trt = unique_treatment_times[i]
                temp = ri_df[ri_df.treated_time .== trt, :]
                X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
                Y = convert(Vector{Float64}, temp.diff)
                β = nothing
                try
                    β = (X \ Y) 
                catch e
                    @warn "Direct solve failed, using pseudoinverse: $e"
                    β = pinv(X' * X) * X' * Y
                end
                att_cohort[i] = β[2]
            end
            ri_att[j] = mean(att_cohort)       
        elseif agg == "simple"
            att_simple = Vector{Float64}(undef, nrow(unique_diffs))
            for i in 1:nrow(unique_diffs)
                t = unique_diffs[i,"t"]
                r1 = unique_diffs[i,"r1"]
                temp = ri_df[(ri_df.t .== t) .& (ri_df.r1 .== r1), :]
                X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.diff)))
                Y = convert(Vector{Float64}, temp.treat)
                β = nothing
                try
                    β = (X \ Y) 
                catch e
                    @warn "Direct solve failed, using pseudoinverse: $e"
                    β = pinv(X' * X) * X' * Y
                end
                att_simple[i] = β[2]
            end
            ri_att[j] = mean(att_simple)
        elseif agg == "state"
            att_state = Vector{Float64}(undef, length(new_treated_states))
            for i in eachindex(new_treated_states)
                state = new_treated_states[i]
                trt = new_treated[i]
                temp_treated = ri_df[ri_df.state .== state, :]
                temp_control = ri_df[(ri_df.treat .== 0) .& (ri_df.treated_time .== trt), :]
                temp = vcat(temp_control, temp_treated)
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.diff))
                Y = convert(Vector{Float64}, temp.treat)
                β = nothing
                try
                    β = (X \ Y) 
                catch e
                    @warn "Direct solve failed, using pseudoinverse: $e"
                    β = pinv(X' * X) * X' * Y
                end
                att_state[i] = β[2]
            end
            ri_att[j] = mean(att_state)
        elseif agg == "unweighted"
            X = convert(Matrix{Float64}, hcat(ones(nrow(ri_df)), ri_df.diff))
            Y = convert(Vector{Float64}, ri_df.treat)
            β = nothing
            try
                β = (X \ Y) 
            catch e
                @warn "Direct solve failed, using pseudoinverse: $e"
                β = pinv(X' * X) * X' * Y
            end
            ri_att[j] = β[2]
        end 

        if verbose && j % 100 == 0
            println("Completed $(j) of $(nperm - 1) permutations")
        end
    end
      
    return (sum(abs.(ri_att) .> abs(original_att)) / length(ri_att))
end 

function get_sep_info(date::String)
    sep_positions = findall(c -> c in ['/', '-'], date)
    sep_types = unique(date[i] for i in sep_positions)
    return (sep_positions, sep_types)
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

function ri_common_adoption(diff_df::DataFrame, original_att::Number,
                            nperm::Integer, verbose::Bool)

    # Check total number of possible permutations
    n = nrow(diff_df)
    k = nrow(diff_df[diff_df.trt .== 1,:]) 
    n_unique_perms = binomial(n, k)
    if nperm > n_unique_perms
        @warn "'nperm' was set to $nperm but only $n_unique_perms unique permutations exist. \n 
Setting 'nperm' to $n_unique_perms."
        nperm = n_unique_perms
    end 
    if nperm < 500
        @warn "'nperm' is less than 500!"
    end 

    # create nperm unique permutations
    seen = Set{String}()
    key = join(diff_df.trt, "")
    push!(seen, key)
    i = 1
    while i < nperm
        new_perm = shuffle(diff_df.trt)
        key = join(new_perm, "")
        if !(key in seen)
            diff_df[!, Symbol("trt_randomized_$i")] = new_perm
            push!(seen, key)
            i += 1
        end
    end

    # Run all the regressions with the unique permutations of trt
    ri_att = Vector{Float64}(undef, nperm-1)
    for i in 1:nperm - 1
        X = convert(Matrix{Float64}, hcat(ones(nrow(diff_df)), diff_df[!, Symbol("trt_randomized_$i")]))
        Y = convert(Vector{Float64}, diff_df.diff)
        β = nothing
        try
            β = (X \ Y) 
        catch e
            @warn "Direct solve failed, using pseudoinverse: $e"
            β = pinv(X' * X) * X' * Y
        end
        ri_att[i] = β[2]
        if verbose && i % 100 == 0
            println("Completed $(i) of $(nperm - 1) permutations")
        end
    end
    return (sum(abs.(ri_att) .> abs(original_att)) / length(ri_att))
end