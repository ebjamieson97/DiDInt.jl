
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

function iterative_demean(data_working, ccc, covariates_to_include, staggered_adoption, hc, edgecase)
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
    if edgecase
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