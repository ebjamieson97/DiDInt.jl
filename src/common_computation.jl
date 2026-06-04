
## The following functions have to do with actual computation common to both didint_plot and didint_estimate

function safe_solve(X::Matrix, Y::Vector)
    try
        return X \ Y
    catch e
        @warn "$e"
        return missing
    end
end

function compute_hc_covariance(X::Matrix, resid::Vector, hc::AbstractString)
    
    n, k = size(X)
    XXinv = inv(X' * X)
    
    # Compute hat matrix diagonal if needed for HC2/HC3/HC4
    if hc in ["hc2", "hc3", "hc4"]
        # H = X(X'X)⁻¹X' but we only need diagonal
        # hᵢᵢ = xᵢ'(X'X)⁻¹xᵢ
        h = [X[i,:]' * XXinv * X[i,:] for i in 1:n]
    end
    
    # Construct Ω based on HC type
    if hc == "hc0"
        omega_diag = resid .^ 2
    elseif hc == "hc1"
        omega_diag = (n / (n - k)) .* (resid .^ 2)
    elseif hc == "hc2"
        omega_diag = (resid .^ 2) ./ (1 .- h)
    elseif hc == "hc3"
        omega_diag = (resid .^ 2) ./ ((1 .- h) .^ 2)
    elseif hc == "hc4"
        h_bar = mean(h) 
        delta = min.(Ref(4), h ./ h_bar)
        omega_diag = (resid .^ 2) ./ ((1 .- h) .^ delta)
    end
    
    # Sandwich estimator: (X'X)⁻¹ X'ΩX (X'X)⁻¹
    return XXinv * (X' * (omega_diag .* X)) * XXinv
end

function iterative_demean(data_working, ccc, covariates_to_include, staggered_adoption, hc, edgecase; 
                          agg = nothing, time_to_index = nothing, treatment_times = nothing, match_to_these_dates = nothing,
                          treated_states = nothing, unique_states = nothing, use_pre_controls = nothing)

    affected_cells   = nothing
    do_full_edgecase = false

    if edgecase 
        # Case 1: Common adoption, only 1 control state, everything gets computed as edgecase
        if !staggered_adoption
            do_full_edgecase = true
        else
            # Case 2: Staggered adoption, check all gts
            states = eltype(treated_states)[]
            cohorts = eltype(treatment_times)[]
            ts = eltype(match_to_these_dates)[]
            r1s = eltype(match_to_these_dates)[]
            for (i, trt_state) in enumerate(treated_states)
                g = treatment_times[i]
                r1 = match_to_these_dates[time_to_index[g] - 1]
            
                for t in match_to_these_dates[match_to_these_dates .>= g]
                    push!(states, trt_state)
                    push!(cohorts, g)
                    push!(ts, t)
                    push!(r1s, r1)
                end
            end
            gts = DataFrame(state = states, cohort = cohorts, t = ts, r1 = r1s)
            treated_count = Vector{Int}(undef, nrow(gts))
            observed = Set(zip(data_working.state_71X9yTx, data_working.time_dmG5fpM))
            for (i, row) in enumerate(eachrow(gts))
                treated_count[i] = Int(((row.state, row.t) in observed) && ((row.state, row.r1) in observed))
            end
            gts.treated_count = treated_count
            gts = gts[treated_count .== 1, :]
            if agg == "time"
                gts.periods_since_treatment = [time_to_index[t] - time_to_index[g] for (g,t) in zip(gts.cohort, gts.t)]
            end
            gts_individual = copy(gts)

            treated_time_map = Dict(state => treatment_times[i] for (i, state) in enumerate(treated_states))
            if agg in ["simple", "cohort", "time"]
                    gts = combine(groupby(gts, [:cohort, :t, :r1]), :treated_count => sum => :treated_count)
            end
            control_states_names = Vector{Vector{eltype(unique_states)}}(undef, nrow(gts))
            control_count = Vector{Int}(undef, nrow(gts))
            for (i, row) in enumerate(eachrow(gts))
                cs_candidates = filter(s -> get(treated_time_map, s, nothing) != row.cohort, unique_states)
                if !use_pre_controls
                    cs_candidates = cs_candidates[.!in.(cs_candidates, Ref(treated_states))]
                end
                valid = eltype(unique_states)[]
                for cs in cs_candidates
                    ok = if haskey(treated_time_map, cs)
                        ((cs, row.t) in observed) && ((cs, row.r1) in observed) && (row.t < treated_time_map[cs])
                    else
                        ((cs, row.t) in observed) && ((cs, row.r1) in observed)
                    end
                    ok && push!(valid, cs)
                end
                control_count[i] = length(valid)
                control_states_names[i] = valid
            end
            gts.control_count = control_count
            gts.control_states = control_states_names
            gts = gts[control_count .>= 1, :]

            # The check for time agg actually works out since the edgecase is also the only time where the cohort dummy vars arent used in the diff regression
            # that is, when there is only one control long diff and one treated long diff at a specific (g,t)
            if agg == "time"
                gts.periods_since_treatment = [time_to_index[t] - time_to_index[g] for (g,t) in zip(gts.cohort, gts.t)]
                gts_check = combine(groupby(gts, [:periods_since_treatment]), :treated_count => sum => :treated_count, :control_count => sum => :control_count)
            elseif agg == "cohort"
                gts_check = combine(groupby(gts, [:cohort]), :treated_count => sum => :treated_count, :control_count => sum => :control_count)
            elseif agg == "state"
                gts_check = combine(groupby(gts, [:state, :cohort]), :treated_count => sum => :treated_count, :control_count => sum => :control_count)
            else
                gts_check = gts
            end
            long_diff_count = gts_check.control_count .+ gts_check.treated_count
            gts_check = gts_check[long_diff_count .== 2, :]

            join_key = if agg == "cohort"
                [:cohort]
            elseif agg == "state"
                [:state, :cohort]
            elseif agg in ["simple"]
                [:cohort, :t, :r1]
            elseif agg == "time"
                [:periods_since_treatment]
            else
                [:state, :cohort, :t, :r1]
            end
            gts_saturated = innerjoin(gts_individual, select(gts_check, join_key), on = join_key)
            

            gts_join_key = agg in ["cohort", "simple", "time"] ? [:cohort, :t, :r1] : [:state, :cohort, :t, :r1]
            gts_saturated = innerjoin(gts_saturated, select(gts, [gts_join_key..., :control_states]),
                          on = gts_join_key)

            # Build final (state, time) output
            out_states = eltype(unique_states)[]
            out_times  = eltype(match_to_these_dates)[]
                
            for row in eachrow(gts_saturated)
                push!(out_states, row.state, row.state)
                push!(out_times,  row.r1,    row.t)
                for cs in row.control_states
                    push!(out_states, cs,      cs)
                    push!(out_times,  row.r1,  row.t)
                end
            end

            if !isempty(out_states)
                affected_cells = unique(DataFrame(state = out_states, time = out_times))
                if ccc in ["hom"]
                    do_full_edgecase = true
                end
            end
        end
    end 

    data_working = copy(data_working)
    y = Float64.(data_working.outcome_71X9yTx)

    # Divide data: edgecase cells get lambda+vcov from compute_lambda_edgecase;
    # the remainder get lambda only from the normal procedure below.
    has_edgecase = edgecase && (do_full_edgecase || affected_cells !== nothing)
    skeleton = nothing; vcov_lambda_global = nothing

    if has_edgecase
        if do_full_edgecase
            data_ec     = data_working
            data_normal = data_working[1:0, :]
        elseif ccc == "state"
            m = in.(data_working.state_71X9yTx, Ref(unique(affected_cells.state)))
            data_ec = data_working[m, :]; data_normal = data_working[.!m, :]
        elseif ccc == "time"
            m = in.(data_working.time_dmG5fpM, Ref(unique(affected_cells.time)))
            data_ec = data_working[m, :]; data_normal = data_working[.!m, :]
        elseif ccc == "int"
            ec_pairs = Set(zip(affected_cells.state, affected_cells.time))
            m = [(row.state_71X9yTx, row.time_dmG5fpM) in ec_pairs for row in eachrow(data_working)]
            data_ec = data_working[m, :]; data_normal = data_working[.!m, :]
        end

        # Global skeleton uses time_71X9yTx (String) as the time key so that
        # cell_id_map inside compute_lambda_edgecase matches sub_df.time_71X9yTx
        skeleton = sort!(unique(select(data_working, [:state_71X9yTx, :time_71X9yTx])),
                         [:state_71X9yTx, :time_71X9yTx])
        n_sk = nrow(skeleton)
        skeleton.lambda       = fill(NaN, n_sk)
        skeleton.lambda_index = 1:n_sk
        rename!(skeleton, :state_71X9yTx => :state, :time_71X9yTx => :time)
        vcov_lambda_global = fill(0.0, n_sk, n_sk)

        skeleton, vcov_lambda_global = compute_lambda_edgecase(data_ec, ccc, covariates_to_include,
                                                               skeleton, vcov_lambda_global, hc)
        if isempty(data_normal)
            staggered_adoption && (skeleton.time = Date.(skeleton.time))
            select!(skeleton, [:state, :time, :lambda, :lambda_index])
            skeleton.ccc .= ccc
            return skeleton, vcov_lambda_global
        end
        data_working = data_normal   # normal procedure runs only on non-affected data
    end

    # Within-cell demean y and covariates
    data_working.y_demeaned = Float64.(data_working.outcome_71X9yTx)
    for cov in covariates_to_include
        data_working[!, Symbol(cov * "_demeaned")] = Float64.(data_working[!, Symbol(cov)])
    end

    data_working = transform(groupby(data_working, [:state_71X9yTx, :time_71X9yTx])) do df
        df = copy(df)
        df.y_demeaned = df.y_demeaned .- mean(df.y_demeaned)
        for cov in covariates_to_include
            col = Symbol(cov * "_demeaned")
            df[!, col] = df[!, col] .- mean(skipmissing(df[!, col]))
        end
        df
    end

    demeaned_covs = [cov * "_demeaned" for cov in covariates_to_include]

    cell_means = combine(groupby(data_working, [:state_71X9yTx, :time_71X9yTx])) do df
        row = DataFrame(cell_mean_y = mean(df.outcome_71X9yTx))
        for cov in covariates_to_include
            row[!, Symbol("cell_mean_" * cov)] = [mean(skipmissing(df[!, Symbol(cov)]))]
        end
        row
    end

    # Helper: estimate beta from within-cell-demeaned data, within a grouping
    # No intercept because data is already cell-demeaned
    function estimate_beta(df, demeaned_covs)
        isempty(demeaned_covs) && return Float64[]
        y = Float64.(collect(df.y_demeaned))
        X = Matrix(Float64.(df[:, Symbol.(demeaned_covs)]))
        return X \ y
    end

    if ccc == "int"
        lambda_df = combine(groupby(data_working, [:state_71X9yTx, :time_71X9yTx])) do df
            active_covs = filter(covariates_to_include) do c
                length(unique(skipmissing(df[!, Symbol(c)]))) > 1
            end
            y = Float64.(collect(df.outcome_71X9yTx))
            if isempty(active_covs)
                X = ones(nrow(df), 1)
            else
                X = hcat(ones(nrow(df)), prune_covariates(Matrix(Float64.(df[:, Symbol.(active_covs)]))))
            end
            beta = X \ y
            DataFrame(lambda = [beta[1]])
        end

    elseif ccc == "none"
        lambda_df = combine(groupby(data_working, [:state_71X9yTx, :time_71X9yTx])) do df
            y = Float64.(collect(df.outcome_71X9yTx))
            DataFrame(lambda = [mean(y)])
        end

    elseif ccc == "time"
        # Beta is time-specific - estimate within each time period using cell-demeaned data
        betas = combine(groupby(data_working, :time_71X9yTx)) do df
            b = estimate_beta(df, demeaned_covs)
            DataFrame(beta = isempty(b) ? [zeros(0)] : [b])
        end

    elseif ccc == "state"
        # Beta is state-specific
        betas = combine(groupby(data_working, :state_71X9yTx)) do df
            b = estimate_beta(df, demeaned_covs)
            DataFrame(beta = isempty(b) ? [zeros(0)] : [b])
        end

  elseif ccc == "add"
    states = sort(unique(data_working.state_71X9yTx))
    times = sort(unique(data_working.time_71X9yTx))
    n = nrow(data_working)
    ncovs = length(demeaned_covs)
    
    state_idx = Dict(s => i for (i,s) in enumerate(states))
    time_idx = Dict(t => i for (i,t) in enumerate(times))

    I_rows = Int[]
    J_cols = Int[]
    V_vals = Float64[]

    for (row_i, row) in enumerate(eachrow(data_working))
        s = state_idx[row.state_71X9yTx]
        t = time_idx[row.time_71X9yTx]
        for (j, cov) in enumerate(demeaned_covs)
            x_val = Float64(row[Symbol(cov)])
            push!(I_rows, row_i); push!(J_cols, (s-1)*ncovs + j);                        push!(V_vals, x_val)
            push!(I_rows, row_i); push!(J_cols, length(states)*ncovs + (t-1)*ncovs + j); push!(V_vals, x_val)
        end
    end

    Z_sparse = sparse(I_rows, J_cols, V_vals, n, (length(states) + length(times)) * ncovs)
    y_dem = Float64.(data_working.y_demeaned)
    β_all = Z_sparse \ y_dem

    β_s_all = reshape(β_all[1:length(states)*ncovs], ncovs, length(states))
    β_t_all = reshape(β_all[length(states)*ncovs+1:end], ncovs, length(times))

    lambda_df = cell_means
    lambda_df.lambda = map(eachrow(lambda_df)) do row
        isempty(covariates_to_include) && return row.cell_mean_y
        xmeans = [row[Symbol("cell_mean_" * cov)] for cov in covariates_to_include]
        s_i = state_idx[row.state_71X9yTx]
        t_i = time_idx[row.time_71X9yTx]
        row.cell_mean_y - dot(β_s_all[:, s_i], xmeans) - dot(β_t_all[:, t_i], xmeans)
    end

    elseif ccc == "hom"
        # Single global beta
        b = isempty(demeaned_covs) ? Float64[] : begin
            y2 = Float64.(data_working.y_demeaned)
            X2 = Matrix(Float64.(data_working[:, Symbol.(demeaned_covs)]))
            X2 \ y2
        end
    end    

    # Apply lambda = cell_mean_y - beta' * cell_mean_X
    if ccc == "time"
        lambda_df = leftjoin(cell_means, betas, on = :time_71X9yTx)
        lambda_df.lambda = map(eachrow(lambda_df)) do row
            isempty(covariates_to_include) && return row.cell_mean_y
            b = row.beta
            xmeans = [row[Symbol("cell_mean_" * cov)] for cov in covariates_to_include]
            row.cell_mean_y - dot(b, xmeans)
        end

    elseif ccc == "state"
        lambda_df = leftjoin(cell_means, betas, on = :state_71X9yTx)
        lambda_df.lambda = map(eachrow(lambda_df)) do row
            isempty(covariates_to_include) && return row.cell_mean_y
            b = row.beta
            xmeans = [row[Symbol("cell_mean_" * cov)] for cov in covariates_to_include]
            row.cell_mean_y - dot(b, xmeans)
        end

    elseif ccc == "hom"
        lambda_df = cell_means
        lambda_df.lambda = map(eachrow(lambda_df)) do row
            isempty(covariates_to_include) && return row.cell_mean_y
            xmeans = [row[Symbol("cell_mean_" * cov)] for cov in covariates_to_include]
            row.cell_mean_y - dot(b, xmeans)
        end
    end

    sort!(lambda_df, [:state_71X9yTx, :time_71X9yTx])
    n = nrow(lambda_df)
    lambda_df.lambda_index = 1:nrow(lambda_df)
    rename!(lambda_df, :state_71X9yTx => :state, :time_71X9yTx => :time)
    select!(lambda_df, [:state, :time, :lambda, :lambda_index])

    # Merge normal-procedure lambdas into global skeleton (both still have String time here)
    if has_edgecase
        cell_map = Dict((skeleton.state[i], skeleton.time[i]) => i for i in 1:nrow(skeleton))
        for row in eachrow(lambda_df)
            skeleton.lambda[cell_map[(row.state, row.time)]] = row.lambda
        end
        staggered_adoption && (skeleton.time = Date.(skeleton.time))
        skeleton.ccc .= ccc
        return skeleton, vcov_lambda_global
    end

    staggered_adoption && (lambda_df.time = Date.(lambda_df.time))
    lambda_df.ccc .= ccc
    vcov_lambda = fill(NaN, n, n)
    return lambda_df, vcov_lambda
end

function compute_lambda_edgecase(data_working, ccc, covariates_to_include, lambda_df, vcov_lambda, hc)
    n_lambda = nrow(lambda_df)

    # lambda_df.time is String matches data_working.time_71X9yTx
    cell_id_map = Dict((lambda_df.state[i], lambda_df.time[i]) =>
                       lambda_df.lambda_index[i] for i in 1:n_lambda)

    # Returns (V_block, global_idx, lambda_vals) where lambda_vals are the (s,t) cell intercepts
    function block_regression(sub_df; int = false)
        n_sub = nrow(sub_df)

        cells_seen = sort(unique([(sub_df.state_71X9yTx[i], sub_df.time_71X9yTx[i])
                                  for i in 1:n_sub]))
        local_idx  = Dict(c => i for (i, c) in enumerate(cells_seen))
        n_cells    = length(cells_seen)
        global_idx = [cell_id_map[c] for c in cells_seen]

        # Dummy block
        D_block = zeros(n_sub, n_cells)
        for i in 1:n_sub
            ci = local_idx[(sub_df.state_71X9yTx[i], sub_df.time_71X9yTx[i])]
            D_block[i, ci] = 1.0
        end

        # Filter to covariates with variation
        active = filter(covariates_to_include) do c
            stds = combine(groupby(sub_df, [:state_71X9yTx, :time_71X9yTx])) do df
                vals = collect(skipmissing(df[!, Symbol(c)]))
                DataFrame(has_var = length(unique(vals)) > 1)
            end
            any(stds.has_var)
        end
        
        if isempty(active)
            Z = D_block
        elseif int
            Z = hcat(D_block, prune_covariates(Matrix(Float64.(sub_df[:, Symbol.(active)]))))
        else
            Z = hcat(D_block, Matrix(Float64.(sub_df[:, Symbol.(active)])))
        end


        y     = Float64.(sub_df.outcome_71X9yTx)
        β = Z \ y
        # If Z is not full rank or not enough rows then in this edgecase of an edgecase we cannot compute the se(ATT)
        if (size(Z, 1) > size(Z, 2)) && (rank(Z) == size(Z, 2))
            resid = y .- Z * β
            V     = compute_hc_covariance(Z, resid, hc)
        else
            V = fill(NaN, n_cells, n_cells)
        end
        return V[1:n_cells, 1:n_cells], global_idx, β[1:n_cells]
    end

    function place_block!(V_block, idx, lambda_vals)
        for j in 1:length(idx), i in 1:length(idx)
            vcov_lambda[idx[i], idx[j]] = V_block[i, j]
        end
        for i in 1:length(idx)
            lambda_df.lambda[idx[i]] = lambda_vals[i]
        end
    end

    if ccc == "int"
        for grp in groupby(data_working, [:state_71X9yTx, :time_71X9yTx])
            place_block!(block_regression(grp, int = true)...)
        end

    elseif ccc == "state"
        for grp in groupby(data_working, :state_71X9yTx)
            place_block!(block_regression(grp)...)
        end

    elseif ccc == "time"
        for grp in groupby(data_working, :time_71X9yTx)
            place_block!(block_regression(grp)...)
        end

    elseif ccc == "hom"
        place_block!(block_regression(data_working)...)

    end

    return lambda_df, vcov_lambda
end

function prune_covariates(W)
    # This function is to be used when calculating lambda with the two-way intersection (int) DID-INT model
    # its necessary since we cant allow for non-uniquely identified lambda values, whereas we can with the
    # other ccc options as the ATTs calculated from the those DID-INT variations ultimately yield a uniquely
    # identified ATT

    Wcol = size(W, 2)
    Wrow = size(W, 1)

    Wqr = qr(W, ColumnNorm())
    m = min(Wcol, Wrow)
    tol = abs(Wqr.factors[1,1]) * eps(Float64) * m
    # Search for first column below tolerance
    rank = something(findfirst(i -> abs(Wqr.factors[i,i]) <= tol, 1:m), m+1) - 1

    # Reorder from most to least informative columns
    W = W[:, Wqr.p]

    if Wcol != rank
        # Keep only the first rank columns in pivoted order
        W = W[:, 1:rank]
    end

    # Account for fact that we are adding an intercept column afterwards
    Wcol = size(W, 2)
    row_surplus = Wrow - Wcol # Needs to be >= 1

    # Prune the least important covariates according to RRQR if theres still a dimensionality issue
    if row_surplus < 1
        nprune = iszero(row_surplus) ? 1 : abs(row_surplus) + 1
        # Drop the last nprune columns in pivoted order (least informative)
        W = W[:, 1:Wcol - nprune]
    end

    return W
end 