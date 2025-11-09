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

function compute_jknife_se(X::Matrix{<:Number}, Y::Vector{<:Number}, original_att::Number)
    
    n = length(Y)
    if n == 1 
        return missing
    end 
    jknife_beta = Vector{Float64}(undef, n)
    ncolx = size(X,2)
    treat_count = sum(X[:,ncolx] .!= 0)
    control_count = sum(X[:,ncolx] .== 0)
    if (treat_count < 2 || control_count < 2) & (ncolx > 1)
        return missing
    end
    for i in eachindex(Y)
        idx = [1:i-1; i+1:size(X, 1)]
        X_sub = X[idx, :]
        Y_sub = Y[idx]
        jknife_beta[i] = (X_sub \ Y_sub)[ncolx]
    end 
    jknife_se = sqrt(sum((jknife_beta .- original_att).^2) * ((n - 1) / n))
    return jknife_se
end 

function final_regression_results(X::Matrix{<:Number}, Y::Vector{<:Number};
                                  W::Vector{T} where T <: Union{Nothing, Number} = [nothing])
    beta_hat = nothing
    beta_hat_cov = nothing
    ncolx = size(X, 2)

    # Run OLS (normally if weights aren't provided, and scale (X,Y) -> (Xw,Yw) otherwise)
    if eltype(W) <: Nothing
        try
            beta_hat = (X \ Y) 
        catch e
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat = pinv(X' * X) * X' * Y
        end
        resid = Y - X * beta_hat
        omega = Diagonal(resid .^ 2)
        try
            beta_hat_cov = inv(X' * X) * (X' * omega * X) * inv(X' * X)
        catch e
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat_cov = pinv(X' * X) * (X' * omega * X) * pinv(X' * X)
        end 
        beta_hat_se_jknife = compute_jknife_se(X, Y, beta_hat[ncolx]) 
    elseif eltype(W) <: Number
        sw = sqrt.(W)           
        Xw = X .* sw            
        Yw = Y .* sw
        try
            beta_hat = (Xw) \ (Yw) 
        catch e
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat = pinv(Xw' * Xw) * Xw' * Yw
        end
        resid_w = Yw - Xw * beta_hat
        Ωw = Diagonal(resid_w.^2)
        try
            beta_hat_cov = inv(Xw'Xw) * (Xw' * Ωw * Xw) * inv(Xw'Xw)
        catch e 
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat_cov = pinv(Xw'Xw) * (Xw' * Ωw * Xw) * pinv(Xw'Xw)
        end 
        beta_hat_se_jknife = compute_jknife_se(Xw, Yw, beta_hat[ncolx]) 
    end 
    beta_hat_var = diag(beta_hat_cov)
    beta_hat_se = sqrt(beta_hat_var[ncolx]) 
    dof = length(Y) - ncolx
    pval_att = dof > 0 ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se))) : missing 
    pval_att_jknife = dof > 0 && !ismissing(beta_hat_se_jknife) ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se_jknife))) : missing
    result_dict = Dict("beta_hat" => beta_hat[ncolx], "beta_hat_se" => beta_hat_se, "pval_att" => pval_att,
                       "beta_hat_se_jknife" => beta_hat_se_jknife, "pval_att_jknife" => pval_att_jknife)
    return result_dict
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

function randomization_inference_v2(diff_df::DataFrame, nperm::Int, results::DataFrame,
                                    agg::AbstractString, verbose::Bool, seed::Number,
                                    data::DataFrame, weighting::AbstractString,
                                    use_pre_controls::Bool;
                                    dummy_cols::Union{Vector{Symbol}, Nothing} = nothing)
    
    # PART ONE: CREATE RANDOMIZED TREATMENT COLUMNS
    original_treated = unique(diff_df[diff_df.treat .== 1, [:state, :treated_time]])
    k = nrow(original_treated)  
    treatment_times = original_treated.treated_time
    treatment_states = original_treated.state
    all_states = unique(diff_df.state)

    n_unique_perms = compute_n_unique_assignments(treatment_times, length(all_states))
    if nperm > n_unique_perms
        @warn "'nperm' was set to $nperm but only $n_unique_perms unique permutations exist. \n 
                Setting 'nperm' to $(n_unique_perms - 1)."
        nperm = n_unique_perms
    end 
    if nperm < 500
        @warn "'nperm' is less than 500!"
    end 

    randomized_diff_df = diff_df
    Random.seed!(seed)
    
    i = 1
    seen = Set{String}()
    pairs = zip(original_treated.state, original_treated.treated_time)
    key = join(sort([string(s, "-", t) for (s, t) in pairs]), "")
    push!(seen, key)
    while i < nperm
        shuffled_states = shuffle(all_states)
        new_treated_states = shuffled_states[1:k]
        new_treated_times = treatment_times[randperm(k)]
        pairs = zip(new_treated_states, new_treated_times)
        key = join(sort([string(s, "-", t) for (s, t) in pairs]), "")
        if key in seen
            continue
        end 
        assigned_tt = Dict(new_treated_states[j] => new_treated_times[j] for j in 1:k)
        new_treat = Vector{Int}(undef, nrow(diff_df))
        for row_idx in 1:nrow(diff_df)
            state = diff_df[row_idx, :state]
            trt_time = diff_df[row_idx, :treated_time]
            if haskey(assigned_tt, state)
                if trt_time == assigned_tt[state]
                    new_treat[row_idx] = 1
                else 
                    if use_pre_controls
                        t = diff_df[row_idx, :t]
                        if t < assigned_tt[state]
                            new_treat[row_idx] = 0
                        else
                            new_treat[row_idx] = -1
                        end 
                    else
                        new_treat[row_idx] = -1
                    end
                end
            else
                new_treat[row_idx] = 0
            end
        end
        randomized_diff_df[!, Symbol("treat_random_", i)] = new_treat
        i += 1
    end

    # PART TWO: COMPUTE RI_ATT & RI_ATT_SUBGROUP
    # Force agg arguments for common adoption to either none or state
    if (length(unique(treatment_times)) == 1 && in(agg, ["cohort", "gt"])) || (length(unique(treatment_times)) == 1 && length(unique(treatment_states)) == 1)
        agg = "none"
    end 
    if length(unique(treatment_times)) == 1 && agg == "sgt"
        agg = "state"
    end 

    att_ri = Vector{Float64}(undef, nperm - 1)
    if agg == "cohort" 
        att_ri_cohort = Matrix{Float64}(undef, nperm - 1, length(treatment_times))
        for j in 1:nperm - 1 
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, length(treatment_times))
            for i in eachindex(treatment_times)
                # Compute sub aggregate ATT
                trt = treatment_times[i]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.treated_time .== trt), :]
                att_ri_cohort[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)
                
                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / length(W))
            end 
            att_ri[j] = dot(W, att_ri_cohort[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    elseif agg == "state"
        att_ri_state = Matrix{Float64}(undef, nperm - 1, length(treatment_times))
        for j in 1:nperm - 1 
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, length(treatment_states))
            for i in eachindex(treatment_states)
                # Compute sub aggregate ATT
                trt = treatment_times[i]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.treated_time .== trt), :]
                temp_treated_silos = unique(temp[temp[!, colname] .== 1, "state"])
                temp_treated_silo = shuffle(temp_treated_silos)[1]
                temp = temp[(temp.state .== temp_treated_silo) .|| (temp[!, colname] .== 0), :]
                att_ri_state[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)

                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / length(treatment_states))
            end 
            att_ri[j] = dot(W, att_ri_state[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end 
        end 
    elseif agg == "simple"
        unique_diffs = unique(select(diff_df[diff_df.treat .== 1,:], :t, :r1, :treated_time))
        att_ri_simple = Matrix{Float64}(undef, nperm - 1, nrow(unique_diffs))
        for j in 1:nperm - 1
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, nrow(unique_diffs))
            for i in 1:nrow(unique_diffs)
                # Compute sub aggregate ATT
                t = unique_diffs[i,"t"]
                r1 = unique_diffs[i,"r1"]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.t .== t) .&& (diff_df.r1 .== r1), :]
                att_ri_simple[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)

                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / nrow(unique_diffs))
            end 
            att_ri[j] = dot(W, att_ri_simple[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    elseif agg == "sgt"
        unique_diffs = unique(select(diff_df[diff_df.treat .== 1,:], :state, :t, :r1, :treated_time))
        att_ri_sgt = Matrix{Float64}(undef, nperm - 1, nrow(unique_diffs))
        for j in 1:nperm - 1
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, nrow(unique_diffs))
            for i in 1:nrow(unique_diffs)
                # Compute sub aggregate ATT
                t = unique_diffs[i,"t"]
                r1 = unique_diffs[i,"r1"]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.t .== t) .&& (diff_df.r1 .== r1), :]
                temp_treated_silos = unique(temp[temp[!, colname] .== 1, "state"])
                temp_treated_silo = shuffle(temp_treated_silos)[1]
                temp = temp[(temp.state .== temp_treated_silo) .|| (temp[!, colname] .== 0), :]
                att_ri_sgt[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)
                
                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / nrow(unique_diffs))
            end 
            att_ri[j] = dot(W, att_ri_sgt[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    elseif agg == "none"
        for j in 1:nperm - 1
            colname = Symbol("treat_random_$(j)")
            temp = diff_df[diff_df[!, colname] .!= -1,:]
            att_ri[j] = compute_ri_sub_agg_att(temp, weighting, colname)
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    elseif agg == "time"
        times = results.periods_post_treat
        att_ri_time = Matrix{Float64}(undef, nperm - 1, length(times))
        for j in 1:nperm - 1
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, length(times))
            for i in eachindex(times)
                # Compute sub aggregate ATT
                t = times[i]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.time_since_treatment .== t), :]
                Y = convert(Vector{Float64}, temp.diff)
                X = design_matrix_time_agg(temp, dummy_cols, colname)
                if in(weighting, ["diff", "both"]) 
                    W_diff = convert(Vector{Float64}, temp.n)
                    W_diff ./= sum(W_diff) 
                    sw = sqrt.(W_diff)           
                    X = X .* sw            
                    Y = Y .* sw
                end 
                att_ri_time[j,i] = (X \ Y)[end]
                
                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / length(times))
            end 
            att_ri[j] = dot(W, att_ri_time[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end
    end

    # PART THREE: COMPUTE P-VALS BASED ON RI_ATT & RI_ATT_SUBGROUP
    agg_att = results.agg_att[1]
    if agg == "cohort"
        for i in eachindex(treatment_times)
            sub_agg_att = results[results.treatment_time .== treatment_times[i], "att_cohort"][1]
            results[results.treatment_time .== treatment_times[i], "ri_pval_att_cohort"] .= (sum(abs.(att_ri_cohort[:,i]) .> abs(sub_agg_att))) / length(att_ri_cohort[:,i])
        end
    elseif agg == "state"
        for i in eachindex(treatment_states)
            sub_agg_att = results[results.state .== treatment_states[i], "att_s"][1]
            results[results.state .== treatment_states[i], "ri_pval_att_s"] .= (sum(abs.(att_ri_state[:,i]) .> abs(sub_agg_att))) / length(att_ri_state[:,i])
        end 
    elseif agg == "simple"
        for i in 1:nrow(unique_diffs)
            t = unique_diffs[i,"t"]
            r1 = unique_diffs[i,"r1"]
            gvar = unique_diffs[i, "treated_time"]
            sub_agg_att = results[(results.time .== t) .&& (results.r1 .== r1) .&& (results.gvar .== gvar), "att_gt"][1]
            results[(results.time .== t) .&& (results.r1 .== r1) .&& (results.gvar .== gvar), "ri_pval_att_gt"] .= (sum(abs.(att_ri_simple[:,i]) .> abs(sub_agg_att))) / length(att_ri_simple[:,i])
        end
    elseif agg == "sgt"
        for i in 1:nrow(unique_diffs)
            t = unique_diffs[i,"t"]
            gvar = unique_diffs[i, "treated_time"]
            state = unique_diffs[i, "state"]
            sub_agg_att = results[(results.t .== t) .&& (results.gvar .== gvar) .&& (results.state .== state), "att_sgt"][1]
            results[(results.t .== t) .&& (results.gvar .== gvar) .&& (results.state .== state), "ri_pval_att_sgt"] .= (sum(abs.(att_ri_sgt[:,i]) .> abs(sub_agg_att))) / length(att_ri_sgt[:,i])
        end 
    elseif agg == "time"
        for i in eachindex(times)
            sub_agg_att = results[results.periods_post_treat .== times[i], "att_t"][1]
            results[results.periods_post_treat .== times[i], "ri_pval_att_t"] .= (sum(abs.(att_ri_time[:,i]) .> abs(sub_agg_att))) / length(att_ri_time[:,i])
        end 
    end
    results.ri_pval_agg_att[1] = ((sum(abs.(att_ri) .> abs(agg_att))) / length(att_ri))
    results.nperm[1] = nperm - 1
    return results
end

function compute_n_unique_assignments(treatment_times::Vector, total_n_states::Number)
    # This computes the combinations formula * the multinomial coefficient
    n_assignments = length(treatment_times)
    unique_assignments = unique(treatment_times)
    num = factorial(big(total_n_states))           
    den = factorial(big(total_n_states - n_assignments))
    for m in unique_assignments
        n_m = sum(treatment_times .== m)
        den *= factorial(big(n_m))
    end
    return num ÷ den                  
end

function custom_sort_order(s)
    parsed = tryparse(Int, s)
    if isnothing(parsed)
        (0, lowercase(s))
    else
        (1, parsed)
    end 
end 

function compute_ri_sub_agg_att(temp::DataFrame, weighting::AbstractString, colname::Symbol)

    X = convert(Vector{Float64}, temp[!, colname])
    Y = convert(Vector{Float64}, temp.diff)
    if in(weighting, ["both", "diff"])
        W_diff  = convert(Vector{Float64}, temp.n)
        W_diff ./= sum(W_diff)
        sub_agg_att = (dot(W_diff[X .== 1], Y[X .== 1]) / sum(W_diff[X .== 1])) - 
                             (dot(W_diff[X .== 0], Y[X .== 0]) / sum(W_diff[X .== 0]))
    elseif in(weighting, ["att", "none"])
        sub_agg_att = mean(Y[X .== 1]) - mean(Y[X .== 0])
    end 
    return sub_agg_att
end 

function scale_weights_final(results::DataFrame, weighting::AbstractString)
    
    # Scale weights and convert to either Float64 or Nothing
    if weighting in ["att", "both"]
        results.weights ./= sum(results.weights)
        results.weights = convert(Vector{Float64}, results.weights)
    elseif weighting in ["diff", "none"]
        results.weights = convert(Vector{Nothing}, results.weights)
    end
    return results
end 

function design_matrix_time_agg(temp::DataFrame, dummy_cols::Vector{Symbol}, treat_col::Symbol)
    
    # Only keep the dummy columns where at least one row is 1
    active = [c for c in dummy_cols if any(temp[!, c] .== 1)]

    # Drop the first of the kept dummies to avoid perfect multicollinearity with intercept
    dummies_keep = length(active) ≤ 1 ? Symbol[] : active[2:end]

    # If no kept dummy, just omit 
    if isempty(dummies_keep)
        X = hcat(ones(nrow(temp)), temp[:, treat_col])               
    else 
        X = hcat(ones(nrow(temp)), Matrix(temp[:, dummies_keep]), temp[:, treat_col])            
    end
    return X
end

function check_dates_by_state(df::DataFrame; state::Symbol = :state_71X9yTx,
                              time::Symbol = :time_71X9yTx)

    states  = unique(df[!, state])
    ref  = Set(df[df[!, state] .== states[1], time])   
    for s in states[2:end]                                
        if Set(df[df[!, state] .== s, time]) ≠ ref
            @warn "State $s has a different set of time values than $(states[1]). Consider setting options for \n'freq', 'start_time', and 'end_time'."
            return false
        end
    end
    return true
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