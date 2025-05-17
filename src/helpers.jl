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
        error("Er01: Please specify a date_format listed here: $possible_formats.")
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
        error("Er02: Unsupported period type: $period_type, try day(s), week(s), month(s), or year(s).")
    end
end

function compute_jknife_se(X::Matrix{<:Number}, Y::Vector{<:Number}, original_att::Number)
    
    n = length(Y)
    if n == 1 
        return missing
    end 
    jknife_beta = Vector{Float64}(undef, n)
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
            beta_hat = (Xw'Xw) \ (Xw'Yw) 
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
    pval_att_jknife = dof > 0 ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se_jknife))) : missing
    result_dict = Dict("beta_hat" => beta_hat[ncolx], "beta_hat_se" => beta_hat_se, "pval_att" => pval_att,
                       "beta_hat_se_jknife" => beta_hat_se_jknife, "pval_att_jknife" => pval_att_jknife)
    return result_dict
end 

function randomization_inference_didint(diff_df::DataFrame, agg::AbstractString,
                                        original_att::Number, nperm::Integer, control_states::AbstractVector,
                                        treated_states::AbstractVector, verbose::Bool; warnings::Bool = true)

    # Check total number of possible permutations
    n = length(control_states) + length(treated_states)
    k = length(treated_states)
    n_unique_perms = binomial(n, k)
    if nperm > n_unique_perms
        if warnings 
            @warn "'nperm' was set to $nperm but only $n_unique_perms unique permutations exist. \n 
Setting 'nperm' to $n_unique_perms."
        end 
        nperm = n_unique_perms
    end 
    if nperm < 500 && warnings
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

    # Make sure unique_diffs is defined, no need to repeatedly redefine during the loop
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
                X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
                Y = convert(Vector{Float64}, temp.diff)
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
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
                Y = convert(Vector{Float64}, temp.diff)
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
        elseif agg == "none"
            X = convert(Matrix{Float64}, hcat(ones(nrow(ri_df)), ri_df.treat))
            Y = convert(Vector{Float64}, ri_df.diff)
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
    pval = (sum(abs.(ri_att) .> abs(original_att)) / length(ri_att))
    result_dict = Dict("ri_pval" => pval, "ri_nperm" => nperm)
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

function randomization_inference_v2(diff_df::DataFrame, nperm::Int, results::DataFrame,
                                    agg::AbstractString, verbose::Bool, seed::Number,
                                    data::DataFrame, weighting::AbstractString)
    
    # PART ONE: CREATE RANDOMIZED TREATMENT COLUMNS
    original_treated = unique(diff_df[diff_df.treat .== 1, [:state, :treated_time]])
    k = nrow(original_treated)  
    treatment_times = original_treated.treated_time
    treatment_states = original_treated.state
    all_states = unique(diff_df.state)

    n_unique_perms = compute_n_unique_assignments(treatment_times, length(all_states))
    if nperm > n_unique_perms
        @warn "'nperm' was set to $nperm but only $n_unique_perms unique permutations exist. \n 
                Setting 'nperm' to $n_unique_perms."
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
                    new_treat[row_idx] = -1
                end
            else
                new_treat[row_idx] = 0
            end
        end
        randomized_diff_df[!, Symbol("treat_random_", i)] = new_treat
        i += 1
    end

    # PART TWO: COMPUTE RI_ATT & RI_ATT_SUBGROUP
    if length(unique(treatment_times)) == 1
        if agg in ["cohort", "simple"]
            agg = "none"
        end
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
                X = convert(Vector{Float64}, temp[!, colname])
                Y = convert(Vector{Float64}, temp.diff)
                att_ri_cohort[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])
                
                # Compute weights
                if weighting == "default"
                    matched_states = Set(temp[(temp[!, colname] .== 1) .&& (temp.treated_time .== trt), :].state)
                    count = sum((data.time_dmG5fpM .>= trt) .&& in.(data.state_71X9yTx, Ref(matched_states)))
                    W[i] = count
                end 
            end

            # Compute aggregate ATT
            if weighting == "default"
                W ./= sum(W)
            elseif weighting == "equal"
                W .= (1 / length(treatment_times))
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
                X = convert(Vector{Float64}, temp[!, colname])
                Y = convert(Vector{Float64}, temp.diff)
                att_ri_state[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])

                # Compute weights
                if weighting == "default"
                    count = sum((data.time_dmG5fpM .>= trt) .&& (data.state_71X9yTx .== temp_treated_silo))
                    W[i] = count
                end 
            end

            # Compute aggregate ATT
            if weighting == "default"
                W ./= sum(W)
            elseif weighting == "equal"
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
                X = convert(Vector{Float64}, temp[!, colname])
                Y = convert(Vector{Float64}, temp.diff)
                att_ri_simple[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])

                # Compute weights
                if weighting == "default"
                    matched_states = Set(temp[(temp[!, colname] .== 1) .&& (temp.t .== t) .&& (temp.r1 .== r1), :].state)
                    count = sum((data.time_dmG5fpM .== t) .&& in.(data.state_71X9yTx, Ref(matched_states)))
                    W[i] = count
                end 
            end

            # Compute aggregate ATT
            if weighting == "default"
                W ./= sum(W)
            elseif weighting == "equal"
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
                X = convert(Vector{Float64}, temp[!, colname])
                Y = convert(Vector{Float64}, temp.diff)
                att_ri_sgt[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])
                
                # Compute weights
                if weighting == "default"
                    count = sum((data.time_dmG5fpM .== t) .&& (data.state_71X9yTx .== temp_treated_silo))
                    W[i] = count
                end 
            end

            # Compute aggregate ATT
            if weighting == "default"
                W ./= sum(W)
            elseif weighting == "equal"
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
            X = convert(Vector{Float64}, temp[!, colname])
            Y = convert(Vector{Float64}, temp.diff)
            att_ri[j] = mean(Y[X .== 1]) - mean(Y[X .== 0])
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

function compute_weights(results::DataFrame, data::DataFrame,
                         agg::AbstractString, weighting::AbstractString,
                         treated_states::Vector{<:AbstractString},
                         treatment_times::Vector{Date})
    
    if weighting == "equal"
        results.weights .= nothing
    elseif weighting == "default"
        results = compute_default_weights(results, data, agg, treated_states, treatment_times)
    end 

    return results
end 

function compute_default_weights(results::DataFrame, data::DataFrame, agg::AbstractString,
                                 treated_states::Vector{<:AbstractString},
                                 treatment_times::Vector{Date})

    if agg == "sgt"
        states = results.state
        times = results.t
        sgt_weights = Vector{Float64}(undef, length(states))
        for (i, s) in enumerate(states)
            t = times[i]
            count = sum((data.time_dmG5fpM .== t) .&& (data.state_71X9yTx .== s))
            sgt_weights[i] = count
        end 
        sgt_weights ./= sum(sgt_weights)
        results.weights = sgt_weights
    elseif agg == "state"
        states = results.state
        state_weights = Vector{Float64}(undef, length(states))
        for (i, s) in enumerate(states)
            idx = findall(x -> x == s, treated_states)
            matched_time = treatment_times[idx]
            count = sum((data.time_dmG5fpM .>= matched_time) .&& (data.state_71X9yTx .== s))
            state_weights[i] = count
        end 
        state_weights ./= sum(state_weights)
        results.weights = state_weights
    elseif agg == "cohort"
        unique_treatment_times = results.treatment_time
        unique_treatment_weights = Vector{Float64}(undef, length(unique_treatment_times))
        for (i, t) in enumerate(unique_treatment_times)
            idxs = findall(x -> x == t, treatment_times)
            matched_states = Set(treated_states[idxs])
            count = sum((data.time_dmG5fpM .>= t) .&& in.(data.state_71X9yTx, Ref(matched_states)))
            unique_treatment_weights[i] = count
        end
        unique_treatment_weights ./= sum(unique_treatment_weights)
        results.weights = unique_treatment_weights
    elseif agg == "simple"
        gt_weights = Vector{Float64}(undef, nrow(results))
        gvars = results.gvar
        times = results.time
        for (i, g) in enumerate(gvars)
            t = times[i]
            idxs = findall(x -> x == g, treatment_times)
            matched_states = Set(treated_states[idxs])
            count = sum((data.time_dmG5fpM .== t) .&& in.(data.state_71X9yTx, Ref(matched_states)))
            gt_weights[i] = count
        end 
        gt_weights ./= sum(gt_weights)
        results.weights = gt_weights
    end 

    return results
end 