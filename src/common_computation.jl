
## The following functions have to do with actual computation common to both didint_plot and didint_estimate
function construct_formula(ccc, covariates_to_include; forplot = false)

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
    elseif ccc == "none" && forplot
            formula_str *= ""
    else 
        error("'ccc' must be set to one of \"int\", \"time\", \"state\", \"add\", or \"hom\".")
    end 
    formula_str *= ")"
    formula_expr = Meta.parse(formula_str)
    formula = eval(formula_expr)

    return formula

end

function run_fixed_effects_model(data_copy, formula, ccc, covariates, covariates_to_include, hc; common_adoption = false, staggered_adoption = true, recover = true,
                                 iterative = true, fem = false)

    # If iterative is true and fem is false
    if iterative == true && fem == false
        return iterative_demean(data_copy, ccc, covariates_to_include, staggered_adoption, hc)
    end

    # For the particular case of ccc == "hom" and recover == true, we can check for and remove collinearity before running the FE model
    if ccc == "hom" && recover == true
        covariates_to_include = drop_invariant_covariates(data_copy, covariates_to_include)
        formula = construct_formula(ccc, covariates_to_include)
    end

    # Call GC.gc() before running big FixedEffectsModels regression
    GC.gc()

    # Run FixedEffectsModels regression (stage 1), capturing any collinearity warnings
    collinear_set = Set{String}()
    log_buffer = IOBuffer()
    logger = SimpleLogger(log_buffer, Logging.Info)
    stage1 = with_logger(logger) do
        reg(data_copy, formula, Vcov.robust(), contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
            save = false)
    end

    # Parse captured log for collinear variable names
    log_output = String(take!(log_buffer))
    for line in split(log_output, "\n")
        if occursin("is collinear", line)
            m = match(r"RHS-variable (.+) is collinear", line)
            if !isnothing(m)
                push!(collinear_set, strip(m.captures[1]))
            end
        end
    end

    # Recover lambdas
    state_time = split.(replace.(coefnames(stage1), "state_time: " => ""), "0IQR7q6Wei7Ejp4e")
    coefs = coef(stage1)
    if ccc == "hom" && !isnothing(covariates)
        keep_mask = map(x -> !(x[1] in covariates), state_time)
        state_time = state_time[keep_mask]
        coefs = coefs[keep_mask]
        keep_idx = findall(keep_mask)
        vcov_lambda = vcov(stage1)[keep_idx, keep_idx]
    else
        vcov_lambda = vcov(stage1)
    end 
    
    lambda_df = DataFrame(state = first.(state_time), time = last.(state_time))
    lambda_df.lambda = coefs
    # positional index into vcov(stage1), note that any recovered will not have an index value
    lambda_df.lambda_index = 1:nrow(lambda_df)   

    expected_df = unique(data_copy[:, [:state_71X9yTx, :time_71X9yTx]])
    expected_df = rename(expected_df, :state_71X9yTx => :state, :time_71X9yTx => :time)
    lambda_df = leftjoin(expected_df, lambda_df, on = [:state, :time])

    # Mark estimated: must be non-missing AND not flagged as collinear by FixedEffectModels
    # Also flag any lambdas that get zero'd out by FE
    # Collinear coefficients are zeroed by FixedEffectModels and logged as Info messages
    lambda_df.estimated = map(row -> begin
        ismissing(row.lambda) && return false
        row.lambda == 0 && recover && return false # Unfortunately not a great way to distinguish between "zero'd out" 0 lambdas, and true 0 lambdas
        key = row.state * "0IQR7q6Wei7Ejp4e" * row.time
        collinear_label = "state_time: " * key
        !(key in collinear_set) && !(collinear_label in collinear_set)
    end, eachrow(lambda_df))

    if recover
        any_missing = any(i -> i == false, lambda_df.estimated)
        lambda_df_missings = lambda_df[lambda_df.estimated .== false, :]
        if any_missing && ccc == "int"
            for row in eachrow(lambda_df_missings)
                s, t = row.state, row.time
                subset_mask = (data_copy.state_71X9yTx .== s) .& (data_copy.time_71X9yTx .== t)
                lambda_df, vcov_lambda = recover_lambdas_subset(lambda_df, data_copy, ccc, covariates_to_include, subset_mask, vcov_lambda)
            end
        elseif any_missing && ccc == "state"
            for s in unique(lambda_df_missings.state)
                subset_mask = data_copy.state_71X9yTx .== s
                lambda_df, vcov_lambda = recover_lambdas_subset(lambda_df, data_copy, ccc, covariates_to_include, subset_mask, vcov_lambda)
            end
        elseif any_missing && ccc == "time"
            for t in unique(lambda_df_missings.time)
                subset_mask = data_copy.time_71X9yTx .== t
                lambda_df, vcov_lambda = recover_lambdas_subset(lambda_df, data_copy, ccc, covariates_to_include, subset_mask, vcov_lambda)
            end
        end
    end

    # Parse the treatment_times and the time column to dates for staggered adoption
    if staggered_adoption
        lambda_df.time = Date.(lambda_df.time)
    end

    # Note the ccc for this loop
    lambda_df.ccc .= ccc

    return lambda_df, vcov_lambda

end

function safe_solve(X::Matrix, Y::Vector)
    try
        return X \ Y
    catch e
        @warn "Unable to compute an ATT value, see $e"
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

function drop_invariant_covariates(temp, covariates)
    # This filters covariates down to just those covariates
    # that have variation within at least one state-time cell
    filter(covariates) do c
        any(
            length(unique(skipmissing(g[:, Symbol(c)]))) > 1
            for g in groupby(temp, [:state_71X9yTx, :time_71X9yTx])
        )
    end
end

function recover_lambdas_subset(lambda_df, data_copy, ccc, covariates, subset_mask, vcov_lambda)
    temp = data_copy[subset_mask, vcat([:state_71X9yTx, :time_71X9yTx, :outcome_71X9yTx, :state_time], Symbol.(covariates))]
    covariates_temp = drop_invariant_covariates(temp, covariates)

    if ccc == "int"
        # Single cell, since data is already filtered down, don't need to put it thru the FE model
        X = isempty(covariates_temp) ? ones(nrow(temp), 1) : hcat(ones(nrow(temp)), Matrix(Float64.(temp[:, Symbol.(covariates_temp)])))
        y = Float64.(temp.outcome_71X9yTx)
        beta = X \ y
        recovered_lambda = beta[1]
        resid = y - X*beta
        # Just need to recover the variance of state-time lambda coefficient (all cov terms are zero for ccc int)
        var_lambda_temp = compute_hc_covariance(X, resid, "hc1")[1, 1]
        mask = (lambda_df.state .== temp.state_71X9yTx[1]) .& (lambda_df.time .== temp.time_71X9yTx[1])
        row_idx = findfirst(mask)
        replace_idx = lambda_df[row_idx, :lambda_index][1]
        vcov_lambda[replace_idx, replace_idx] = var_lambda_temp
        lambda_df[row_idx, :lambda] = recovered_lambda
        lambda_df[row_idx, :estimated] = true
    else
        formula_temp = construct_formula(ccc, covariates_temp)
        GC.gc()
        stage1_temp = reg(temp, formula_temp, Vcov.robust(),
                          contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
                          save = false)
        state_time_temp = split.(replace.(coefnames(stage1_temp), "state_time: " => ""), "0IQR7q6Wei7Ejp4e")
        coefs_temp = coef(stage1_temp)
        lambda_df_temp = DataFrame(state = first.(state_time_temp), time = last.(state_time_temp))
        lambda_df_temp.lambda = coefs_temp

        # Join the lambda_df to the temp_lambda_df in order to get the matching indexes
        vcov_lambda_temp = vcov(stage1_temp)
        lambda_df_temp.lambda_index_temp = 1:nrow(lambda_df_temp)
        replacement_lambda = leftjoin(lambda_df_temp[:, [:state, :time, :lambda_index_temp]],
                                      lambda_df[:, [:state, :time, :lambda_index]], on = [:state, :time])
        valid = .!ismissing.(replacement_lambda.lambda_index)
        orig_idx = Int.(replacement_lambda.lambda_index[valid])
        temp_idx = Int.(replacement_lambda.lambda_index_temp[valid])
        vcov_lambda[orig_idx, orig_idx] .= vcov_lambda_temp[temp_idx, temp_idx]

        for row in eachrow(lambda_df_temp)
            mask = (lambda_df.state .== row.state) .& (lambda_df.time .== row.time)
            row_idx = findfirst(mask)
            lambda_df[row_idx, :lambda]
            lambda_df[row_idx, :lambda] = row.lambda
            lambda_df[row_idx, :estimated] = true
        end
   end

   
        return lambda_df, vcov_lambda

end

function iterative_demean(data_working, ccc, covariates_to_include, staggered_adoption, hc)
    data_working = copy(data_working)
    y = Float64.(data_working.outcome_71X9yTx)
    
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
            X = isempty(active_covs) ? ones(nrow(df), 1) :
                hcat(ones(nrow(df)), Matrix(Float64.(df[:, Symbol.(active_covs)])))
            beta = X \ y
            DataFrame(lambda = [beta[1]])
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

    elseif ccc == "none"
        lambda_df = cell_means
        lambda_df.lambda = lambda_df.cell_mean_y

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
    if hc != "skip"
        vcov_lambda = compute_vcov_lambda(data_working, ccc, covariates_to_include,
                                          lambda_df, hc)
    else
        vcov_lambda = fill(NaN, n, n)
    end


    # Only keep needed columns
    select!(lambda_df, [:state, :time, :lambda, :lambda_index])

    if staggered_adoption
        lambda_df.time = Date.(lambda_df.time)
    end

    lambda_df.ccc .= ccc

    return lambda_df, vcov_lambda
end

function compute_vcov_lambda(data_working, ccc, covariates_to_include, lambda_df, hc)
    n_lambda = nrow(lambda_df)
    ncovs    = length(covariates_to_include)
    vcov_lambda = zeros(n_lambda, n_lambda)

    cell_id_map = Dict((lambda_df.state[i], lambda_df.time[i]) =>
                       lambda_df.lambda_index[i] for i in 1:n_lambda)

    # ------------------------------------------------------------------
    # Drop rows with missing covariate values to mirror what skipmissing
    # gives the demean step.
    # ------------------------------------------------------------------
    if ncovs > 0
        cov_syms = Symbol.(covariates_to_include)
        keep = trues(nrow(data_working))
        for s in cov_syms
            keep .&= .!ismissing.(data_working[!, s])
        end
        data_working = data_working[keep, :]
    end

    # ------------------------------------------------------------------
    # Per-block OLS of y on [D_block | W_block] with within-cell-constant
    # covariates dropped from W_block. Returns the cell-dummy block of the
    # HC-robust vcov plus global lambda_index for each row/col.
    # ------------------------------------------------------------------
    function block_regression(sub_df)
        n_sub = nrow(sub_df)

        cells_seen = sort(unique([(sub_df.state_71X9yTx[i], sub_df.time_71X9yTx[i])
                                  for i in 1:n_sub]))
        local_idx  = Dict(c => i for (i, c) in enumerate(cells_seen))
        n_cells    = length(cells_seen)
        global_idx = [cell_id_map[c] for c in cells_seen]

        D_block = zeros(n_sub, n_cells)
        for i in 1:n_sub
            ci = local_idx[(sub_df.state_71X9yTx[i], sub_df.time_71X9yTx[i])]
            D_block[i, ci] = 1.0
        end

        active = filter(covariates_to_include) do c
            stds = combine(groupby(sub_df, [:state_71X9yTx, :time_71X9yTx])) do df
                vals = collect(skipmissing(df[!, Symbol(c)]))
                DataFrame(has_var = length(unique(vals)) > 1)
            end
            any(stds.has_var)
        end

        W_block = isempty(active) ? zeros(n_sub, 0) :
                  Matrix(Float64.(sub_df[:, Symbol.(active)]))
        Z = hcat(D_block, W_block)

        if size(Z, 1) <= size(Z, 2)
            return fill(NaN, n_cells, n_cells), global_idx
        end

        y     = Float64.(sub_df.outcome_71X9yTx)
        β̂    = Z \ y
        resid = y .- Z * β̂
        V     = compute_hc_covariance(Z, resid, hc)
        return V[1:n_cells, 1:n_cells], global_idx
    end

    function place_block!(V_block, idx)
        for j in 1:length(idx), i in 1:length(idx)
            vcov_lambda[idx[i], idx[j]] = V_block[i, j]
        end
    end

    if ccc == "int"
        for grp in groupby(data_working, [:state_71X9yTx, :time_71X9yTx])
            place_block!(block_regression(grp)...)
        end

    elseif ccc == "state"
        for grp in groupby(data_working, :state_71X9yTx)
            place_block!(block_regression(grp)...)
        end

    elseif ccc == "time"
        for grp in groupby(data_working, :time_71X9yTx)
            place_block!(block_regression(grp)...)
        end

    elseif ccc == "hom" || ccc == "none"
        place_block!(block_regression(data_working)...)

    elseif ccc == "add"
        # ------------------------------------------------------------------
        # state and time effects on β are jointly estimated -- can't decouple.
        # Build Z = [D | state-interacted X (drop state 1) | time-interacted X]
        # and prune any remaining collinear W columns via rank-revealing QR.
        # All D columns are always retained (they're independent by construction).
        # ------------------------------------------------------------------
        n_obs    = nrow(data_working)
        states   = sort(unique(data_working.state_71X9yTx))
        times    = sort(unique(data_working.time_71X9yTx))
        S, T     = length(states), length(times)
        s_idx    = Dict(s => i for (i, s) in enumerate(states))
        t_idx    = Dict(t => i for (i, t) in enumerate(times))

        # D: cell dummies aligned to lambda_df.lambda_index
        D = zeros(n_obs, n_lambda)
        for i in 1:n_obs
            ci = cell_id_map[(data_working.state_71X9yTx[i],
                              data_working.time_71X9yTx[i])]
            D[i, ci] = 1.0
        end

        if ncovs == 0
            Z = D
        else
            X_raw = Matrix(Float64.(data_working[:, Symbol.(covariates_to_include)]))
            # state 1 dropped as reference -> (S - 1 + T) * ncovs columns
            W = zeros(n_obs, (S - 1 + T) * ncovs)
            for i in 1:n_obs
                s = s_idx[data_working.state_71X9yTx[i]]
                t = t_idx[data_working.time_71X9yTx[i]]
                for j in 1:ncovs
                    if s > 1
                        W[i, (s - 2) * ncovs + j] = X_raw[i, j]
                    end
                    W[i, (S - 1) * ncovs + (t - 1) * ncovs + j] = X_raw[i, j]
                end
            end
            Z = hcat(D, W)
        end

        # Rank-revealing QR: prune redundant W columns, keep all D columns.
        if size(Z, 2) > n_lambda
            F   = qr(Z, ColumnNorm())
            tol = size(Z, 1) * eps(Float64) * maximum(abs.(diag(F.R)))
            r   = count(>(tol), abs.(diag(F.R)))
            if r < size(Z, 2)
                keep = sort(union(1:n_lambda, F.p[1:r]))
                Z = Z[:, keep]
            end
        end

        if size(Z, 1) <= size(Z, 2)
            vcov_lambda .= NaN
        else
            y      = Float64.(data_working.outcome_71X9yTx)
            α̂     = Z \ y
            resid  = y .- Z * α̂
            V_full = compute_hc_covariance(Z, resid, hc)
            vcov_lambda = V_full[1:n_lambda, 1:n_lambda]
        end
    end

    return vcov_lambda
end