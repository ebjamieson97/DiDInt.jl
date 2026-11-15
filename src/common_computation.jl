
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

    # A row with leverage 1 (such as the single treated row of an ATT(s,t) regression) leaves
    # HC2-HC4 undefined. Rounding can put h a hair below 1 and turn 0/0 into an arbitrary
    # finite number, so such rows are set to NaN and the SE is NaN (update spec D17).
    if hc in ["hc2", "hc3", "hc4"]
        omega_diag[h .> 1 - 1e-10] .= NaN
    end

    # Sandwich estimator: (X'X)⁻¹ X'ΩX (X'X)⁻¹
    return XXinv * (X' * (omega_diag .* X)) * XXinv
end

function hc_mean_covariance(Z::Matrix, hc::AbstractString)

    # HC covariance matrix of the column means of Z: the multivariate version of the HC variance
    # of the intercept in a regression on a constant (every row has leverage 1/n), so for one
    # column it equals compute_hc_covariance(ones(n, 1), z .- mean(z), hc)
    n = size(Z, 1)
    centered = Z .- mean(Z, dims = 1)
    h = 1 / n
    scale = hc == "hc0" ? 1.0 : hc == "hc1" ? n / (n - 1) : hc == "hc3" ? 1 / (1 - h)^2 : 1 / (1 - h)
    return scale .* (centered' * centered) ./ n^2
end

# First-stage estimates beyond the cell intercepts lambda (update spec §3.1), keyed by cell
# (state, period). `slopes` holds the covariate slopes gamma_{s,t} implied by the ccc, in the
# order of `covariates` (0 for a covariate dropped in that cell), and `xbar` the cell's covariate
# means. For the cells of an edge case under ccc = "int", `cell_vcov` holds the HC covariance of
# (lambda, gamma) from the cell's regression and `cell_obs` the cell's outcome vector and
# covariate matrix; edge_case_se_integrated needs both (spec Appendix B.4).
struct CovariateFit
    covariates::Vector{String}
    slopes::Dict{Tuple{String, Date}, Vector{Float64}}
    xbar::Dict{Tuple{String, Date}, Vector{Float64}}
    cell_vcov::Dict{Tuple{String, Date}, Matrix{Float64}}
    cell_obs::Dict{Tuple{String, Date}, Tuple{Vector{Float64}, Matrix{Float64}}}
end

function CovariateFit(covariates::Vector{String})
    return CovariateFit(covariates, Dict{Tuple{String, Date}, Vector{Float64}}(), Dict{Tuple{String, Date}, Vector{Float64}}(),
                        Dict{Tuple{String, Date}, Matrix{Float64}}(),
                        Dict{Tuple{String, Date}, Tuple{Vector{Float64}, Matrix{Float64}}}())
end

function iterative_demean(data_working, ccc, covariates_to_include, hc, edgecase;
                          time_to_index = nothing, treatment_times = nothing, match_to_these_dates = nothing,
                          treated_states = nothing, use_pre_controls = nothing)

    affected_cells   = nothing
    do_full_edgecase = false

    if edgecase
        # Edge case (update spec §3.6): an ATT(s,t) whose regression has exactly two rows, the
        # treated state s and a single valid control c, whatever the agg. Its SE comes from the
        # covariances of the lambdas of the cells (s, r1), (s, t), (c, r1) and (c, t), so those
        # cells are collected here and estimated by compute_lambda_edgecase below.
        observed = Set(zip(data_working.state_71X9yTx, data_working.time_dmG5fpM))
        states = unique(data_working.state_71X9yTx)
        treated_time_map = Dict(state => treatment_times[i] for (i, state) in enumerate(treated_states))
        out_states = eltype(states)[]
        out_times  = eltype(match_to_these_dates)[]
        for (i, trt_state) in enumerate(treated_states)
            g = treatment_times[i]
            r1 = match_to_these_dates[time_to_index[g] - 1]
            for t in match_to_these_dates[match_to_these_dates .>= g]
                ((trt_state, t) in observed && (trt_state, r1) in observed) || continue

                # Valid controls: never-treated states and, with use_pre_controls, states not yet
                # treated at t, observed in both t and r1
                valid_controls = eltype(states)[]
                for cs in states
                    cs_time = get(treated_time_map, cs, nothing)
                    is_control = isnothing(cs_time) || (use_pre_controls && t < cs_time)
                    if is_control && ((cs, t) in observed) && ((cs, r1) in observed)
                        push!(valid_controls, cs)
                    end
                end
                if length(valid_controls) == 1
                    cs = valid_controls[1]
                    push!(out_states, trt_state, trt_state, cs, cs)
                    push!(out_times,  r1,        t,         r1, t)
                end
            end
        end

        if !isempty(out_states)
            affected_cells = unique(DataFrame(state = out_states, time = out_times))
            if ccc in ["hom"]
                do_full_edgecase = true
            end
        end
    end

    data_working = copy(data_working)
    y = Float64.(data_working.outcome_71X9yTx)

    # Covariate slopes and cell means of the covariates (update spec §3.1). The means come from
    # the same rows as the first stage; the slopes are recorded where each ccc estimates them.
    fit = CovariateFit(String[string(cov) for cov in covariates_to_include])
    for cell in groupby(data_working, [:state_71X9yTx, :time_dmG5fpM])
        fit.xbar[(cell.state_71X9yTx[1], cell.time_dmG5fpM[1])] = Float64[mean(cell[!, Symbol(cov)]) for cov in covariates_to_include]
    end

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
                                                               skeleton, vcov_lambda_global, hc, fit)
        if isempty(data_normal)
            skeleton.time = Date.(skeleton.time)
            select!(skeleton, [:state, :time, :lambda, :lambda_index])
            skeleton.ccc .= ccc
            return skeleton, vcov_lambda_global, fit
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
            gamma = zeros(length(covariates_to_include))
            if isempty(active_covs)
                X = ones(nrow(df), 1)
                beta = X \ y
            else
                W = Matrix(Float64.(df[:, Symbol.(active_covs)]))
                kept = prune_covariate_columns(W)
                X = hcat(ones(nrow(df)), W[:, kept])
                beta = X \ y
                gamma[Int.(indexin(active_covs[kept], covariates_to_include))] = beta[2:end]
            end
            # The cell's own slopes gamma_{s,t}; a covariate dropped in this cell keeps slope 0
            fit.slopes[(df.state_71X9yTx[1], df.time_dmG5fpM[1])] = gamma
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
    # Without covariates there is nothing to solve for (and the sparse solver rejects an empty design)
    β_all = ncovs == 0 ? Float64[] : Z_sparse \ y_dem

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

    # Slopes gamma_{s,t} implied by the ccc for every cell (update spec §3.1): gamma for "hom",
    # gamma_s for "state", gamma_t for "time", gamma_s + gamma_t for "add". ccc = "int" recorded
    # its cell-specific slopes above.
    if ccc != "int"
        no_slopes = isempty(covariates_to_include) || ccc == "none"
        if !no_slopes && ccc in ["time", "state"]
            group_col = ccc == "time" ? :time_71X9yTx : :state_71X9yTx
            beta_by_group = Dict(zip(betas[!, group_col], betas.beta))
        end
        for cell in eachrow(unique(select(data_working, [:state_71X9yTx, :time_71X9yTx, :time_dmG5fpM])))
            fit.slopes[(cell.state_71X9yTx, cell.time_dmG5fpM)] =
                if no_slopes
                    zeros(length(covariates_to_include))
                elseif ccc == "hom"
                    b
                elseif ccc == "add"
                    β_s_all[:, state_idx[cell.state_71X9yTx]] .+ β_t_all[:, time_idx[cell.time_71X9yTx]]
                else
                    beta_by_group[cell[group_col]]
                end
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
        skeleton.time = Date.(skeleton.time)
        skeleton.ccc .= ccc
        return skeleton, vcov_lambda_global, fit
    end

    lambda_df.time = Date.(lambda_df.time)
    lambda_df.ccc .= ccc
    vcov_lambda = fill(NaN, n, n)
    return lambda_df, vcov_lambda, fit
end

function compute_lambda_edgecase(data_working, ccc, covariates_to_include, lambda_df, vcov_lambda, hc, fit::CovariateFit)
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
        
        # Covariates in the regression: those with variation, pruned further under "int"
        kept_covs = active
        if isempty(active)
            Z = D_block
        elseif int
            W = Matrix(Float64.(sub_df[:, Symbol.(active)]))
            kept = prune_covariate_columns(W)
            kept_covs = active[kept]
            Z = hcat(D_block, W[:, kept])
        else
            Z = hcat(D_block, Matrix(Float64.(sub_df[:, Symbol.(active)])))
        end


        y     = Float64.(sub_df.outcome_71X9yTx)
        β = Z \ y
        # If Z is not full rank or not enough rows then in this edgecase of an edgecase we cannot compute the se(ATT)
        estimable = (size(Z, 1) > size(Z, 2)) && (rank(Z) == size(Z, 2))
        if estimable
            resid = y .- Z * β
            V     = compute_hc_covariance(Z, resid, hc)
        else
            V = fill(NaN, n_cells, n_cells)
        end

        # Record the slopes of every cell in the block (update spec §3.1): common to the block
        # for "hom", "state" and "time", the cell's own for "int"
        kept_positions = Int.(indexin(kept_covs, covariates_to_include))
        gamma = zeros(length(covariates_to_include))
        gamma[kept_positions] = β[n_cells+1:end]
        for cell in cells_seen
            fit.slopes[(cell[1], Date(cell[2]))] = gamma
        end

        # Under "int" with covariates the block is one cell. Keep the HC covariance of its
        # (lambda, gamma), with zero rows and columns for dropped covariates, and its
        # observations: edge_case_se_integrated needs them (update spec Appendix B.4).
        if int && !isempty(covariates_to_include)
            n_coef = 1 + length(covariates_to_include)
            positions = vcat(1, 1 .+ kept_positions)
            cell_vcov = estimable ? zeros(n_coef, n_coef) : fill(NaN, n_coef, n_coef)
            if estimable
                cell_vcov[positions, positions] = V
            end
            cell_key = (cells_seen[1][1], Date(cells_seen[1][2]))
            fit.cell_vcov[cell_key] = cell_vcov
            fit.cell_obs[cell_key] = (y, Matrix(Float64.(sub_df[:, Symbol.(covariates_to_include)])))
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

function prune_covariate_columns(W)
    # This function is to be used when calculating lambda with the two-way intersection (int) DID-INT model
    # its necessary since we cant allow for non-uniquely identified lambda values, whereas we can with the
    # other ccc options as the ATTs calculated from the those DID-INT variations ultimately yield a uniquely
    # identified ATT
    # Returns the indices of the columns of W to keep, from most to least informative
    # Collinear covariates are dropped (rank-revealing QR), and a dropped covariate gets slope 0
    # in that cell. This is a deliberate choice (confirmed by the package author): didintpy instead
    # keeps every covariate and uses the minimum-norm least-squares solution. The two agree when a
    # covariate is constant within the cell, but can give different slopes, and so different
    # integrated ATT(s,t)'s, for other collinearities.

    Wcol = size(W, 2)
    Wrow = size(W, 1)

    Wqr = qr(W, ColumnNorm())
    m = min(Wcol, Wrow)
    tol = abs(Wqr.factors[1,1]) * eps(Float64) * m
    # Search for first column below tolerance
    rank = something(findfirst(i -> abs(Wqr.factors[i,i]) <= tol, 1:m), m+1) - 1

    # Keep only the first rank columns in pivoted order (most to least informative)
    kept = Wqr.p[1:rank]

    # Account for fact that we are adding an intercept column afterwards
    row_surplus = Wrow - length(kept) # Needs to be >= 1

    # Prune the least important covariates according to RRQR if theres still a dimensionality issue
    if row_surplus < 1
        nprune = iszero(row_surplus) ? 1 : abs(row_surplus) + 1
        # Drop the last nprune columns in pivoted order (least informative)
        kept = kept[1:length(kept) - nprune]
    end

    return kept
end

