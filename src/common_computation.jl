###########################################################################################
# common_computation.jl: the first stage, shared by didint() and didint_plot()            #
#                                                                                         #
#   Section 1: Least squares and HC covariance estimators                                 #
#   Section 2: First-stage results beyond the lambdas (CovariateFit)                      #
#   Section 3: The first stage: lambda_{s,t} for every state and period                   #
#   Section 4: The first stage with covariances, for edge-case standard errors            #
#                                                                                         #
# The first stage of DID-INT (Karim & Webb 2025, eq. 45) is the regression                #
#   Y_i = sum_{s,t} lambda_{s,t} I(s,t) + f(X_i) + e_i,                                   #
# where I(s,t) marks the observations of state s in period t (a "cell") and f(X) is       #
# linear in the covariates with slopes that are common to all cells ("hom"), vary by      #
# state ("state"), by period ("time"), by state plus period ("add") or by cell ("int").   #
# lambda_{s,t} is the covariate-adjusted mean outcome of the cell.                        #
###########################################################################################


###########################################################################################
# Section 1: Least squares and HC covariance estimators                                   #
###########################################################################################

function safe_solve(X::Matrix, Y::Vector)

    # Least-squares coefficients of Y on X, or missing (with a warning) if they cannot be
    # computed. Only failures of the linear algebra are turned into missing; anything else
    # (such as running out of memory) is an error.
    try
        return X \ Y
    catch e
        e isa Union{SingularException, LAPACKException, PosDefException, RankDeficientException} || rethrow()
        @warn "$e"
        return missing
    end
end

function compute_hc_covariance(X::Matrix, resid::Vector, hc::AbstractString)

    # Heteroskedasticity-consistent covariance matrix of least-squares coefficients,
    # (X'X)^-1 X' Omega X (X'X)^-1 with Omega diagonal: HC0 (White 1980), HC1-HC3 (MacKinnon &
    # White 1985) or HC4 (Cribari-Neto 2004)
    n, k = size(X)
    XXinv = inv(X' * X)

    # Leverages (the diagonal of the hat matrix X(X'X)^-1 X'), needed for HC2/HC3/HC4
    if hc in ["hc2", "hc3", "hc4"]
        h = [X[i,:]' * XXinv * X[i,:] for i in 1:n]
    end

    # Diagonal of Omega
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
        # Rounding can put a leverage of 1 slightly above 1; without the max, the negative base
        # raised to a non-integer power would throw an error
        omega_diag = (resid .^ 2) ./ (max.(1 .- h, 0.0) .^ delta)
    end

    # A row with leverage 1 (such as the single treated row of an ATT(s,t) regression) leaves
    # HC2-HC4 undefined. Rounding can put h a hair below 1 and turn 0/0 into an arbitrary
    # finite number, so such rows are set to NaN and the SE is NaN.
    if hc in ["hc2", "hc3", "hc4"]
        omega_diag[h .> 1 - 1e-10] .= NaN
    end

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


###########################################################################################
# Section 2: First-stage results beyond the lambdas (CovariateFit)                        #
###########################################################################################

# First-stage estimates other than the cell intercepts lambda, keyed by cell (state, period).
# `slopes` holds the covariate slopes gamma_{s,t} implied by the ccc, in the order of
# `covariates` (0 for a covariate dropped in that cell), and `xbar` the cell's covariate
# means; the integrated long differences need both (Karim & Webb 2025, eqs. 46-47). For the
# cells of an edge case under ccc = "int", `cell_vcov` holds the HC covariance of (lambda,
# gamma) from the cell's regression and `cell_obs` the cell's outcome vector and covariate
# matrix, which edge_case_se_integrated needs.
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


###########################################################################################
# Section 3: The first stage: lambda_{s,t} for every state and period                     #
###########################################################################################

# Columns that identify a cell: the state and the period (as a string)
const CELL_COLUMNS = [:state_71X9yTx, :time_71X9yTx]

function first_stage(data_working, ccc, covariates_to_include, hc, edgecase;
                     time_to_index = nothing, treatment_times = nothing, match_to_these_dates = nothing,
                     treated_states = nothing, use_pre_controls = nothing)

    # Estimate lambda_{s,t} for every cell under the chosen ccc (Karim & Webb 2025, eq. 45).
    # Returns
    #   lambda_df    one row per cell: state, time, lambda, lambda_index, ccc
    #   vcov_lambda  the covariance matrix of the lambdas, indexed by lambda_index. It is only
    #                computed with `edgecase`, for the cells that an edge-case standard error
    #                needs (zero elsewhere); otherwise it is `nothing`
    #   fit          the covariate slopes and covariate means of every cell (a CovariateFit)
    # The lambdas are not estimated from the one large regression. Under "int" each cell has
    # its own regression. Under "hom", "state", "time" and "add" the slopes are estimated from
    # outcomes and covariates demeaned within cells (the Frisch-Waugh-Lovell theorem) and
    # lambda_{s,t} = ybar_{s,t} - gamma' xbar_{s,t}. Both give the lambdas of the large regression.

    # Covariate means of every cell
    fit = CovariateFit(String[string(cov) for cov in covariates_to_include])
    for cell in groupby(data_working, [:state_71X9yTx, :time_dmG5fpM])
        fit.xbar[(cell.state_71X9yTx[1], cell.time_dmG5fpM[1])] = Float64[mean(cell[!, Symbol(cov)]) for cov in covariates_to_include]
    end

    # With `edgecase`, the cells behind an ATT(s,t) with a single control (if there are any) get
    # their lambdas and the covariances of their lambdas from first_stage_edgecase, together
    # with every cell that shares slopes with them. The other cells go through the regular
    # procedure below, which is cheaper.
    edgecase_cells = edgecase ? find_edgecase_cells(data_working, treated_states, treatment_times, time_to_index,
                                                    match_to_these_dates, use_pre_controls) : nothing
    has_edgecase = !isnothing(edgecase_cells)
    if has_edgecase
        # One row per cell of the whole data, to be filled by both procedures
        skeleton = sort!(unique(select(data_working, CELL_COLUMNS)), CELL_COLUMNS)
        n_cells = nrow(skeleton)
        skeleton.lambda = fill(NaN, n_cells)
        skeleton.lambda_index = 1:n_cells
        rename!(skeleton, :state_71X9yTx => :state, :time_71X9yTx => :time)
        vcov_lambda = fill(0.0, n_cells, n_cells)

        in_edgecase = edgecase_rows(data_working, ccc, edgecase_cells)
        skeleton, vcov_lambda = first_stage_edgecase(data_working[in_edgecase, :], ccc, covariates_to_include,
                                                     skeleton, vcov_lambda, hc, fit)
        data_working = data_working[.!in_edgecase, :]
    end

    # Regular procedure
    if !has_edgecase || !isempty(data_working)
        lambda_df = cell_lambdas!(fit, data_working, ccc, covariates_to_include)
        sort!(lambda_df, CELL_COLUMNS)
        lambda_df.lambda_index = 1:nrow(lambda_df)
        rename!(lambda_df, :state_71X9yTx => :state, :time_71X9yTx => :time)
        select!(lambda_df, [:state, :time, :lambda, :lambda_index])
    end

    # With an edge case, put the regular procedure's lambdas into the table of all cells
    if has_edgecase
        if !isempty(data_working)
            position = Dict((skeleton.state[i], skeleton.time[i]) => i for i in 1:nrow(skeleton))
            for row in eachrow(lambda_df)
                skeleton.lambda[position[(row.state, row.time)]] = row.lambda
            end
        end
        lambda_df = skeleton
    else
        vcov_lambda = nothing
    end

    lambda_df.time = Date.(lambda_df.time)
    lambda_df.ccc .= ccc
    return lambda_df, vcov_lambda, fit
end

function cell_lambdas!(fit::CovariateFit, data, ccc, covariates_to_include)

    # lambda of every cell of `data` by the regular procedure, as a DataFrame with the cell
    # columns and `lambda`. Also records every cell's covariate slopes in `fit`.
    if ccc == "int"
        return cell_lambdas_int!(fit, data, covariates_to_include)
    elseif ccc == "none"
        # No covariate adjustment (used by didint_plot): lambda is the cell mean of the outcome
        return combine(groupby(data, CELL_COLUMNS)) do df
            fit.slopes[(df.state_71X9yTx[1], df.time_dmG5fpM[1])] = zeros(length(covariates_to_include))
            DataFrame(lambda = [mean(Float64.(collect(df.outcome_71X9yTx)))])
        end
    else
        return cell_lambdas_fwl!(fit, data, ccc, covariates_to_include)
    end
end

function cell_lambdas_int!(fit::CovariateFit, data, covariates_to_include)

    # ccc = "int": the slopes are specific to the cell, so each cell has its own regression of
    # the outcome on an intercept and the covariates, and lambda is the intercept. Covariates
    # that do not vary in the cell, or that are collinear, are dropped from that cell's
    # regression and keep slope 0 there (see prune_covariate_columns).
    return combine(groupby(data, CELL_COLUMNS)) do df
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
        fit.slopes[(df.state_71X9yTx[1], df.time_dmG5fpM[1])] = gamma
        DataFrame(lambda = [beta[1]])
    end
end

function cell_lambdas_fwl!(fit::CovariateFit, data, ccc, covariates_to_include)

    # ccc = "hom", "state", "time" or "add": the slopes are shared across cells. By the
    # Frisch-Waugh-Lovell theorem they can be estimated from the outcome and the covariates
    # demeaned within cells (which removes the cell dummies), and then
    #   lambda_{s,t} = ybar_{s,t} - gamma_{s,t}' xbar_{s,t}
    # with gamma_{s,t} the slopes that apply to the cell: gamma ("hom"), gamma_s ("state"),
    # gamma_t ("time") or gamma_s + gamma_t ("add").
    # A covariate that does not vary within a cell (for example one measured at the state-period
    # level) is absorbed by the cell dummies. Its demeaned values are set to exactly zero:
    # subtracting the mean of identical values can leave rounding noise of about 1e-17, and a
    # regression on a column of pure rounding noise returns an arbitrary slope.
    demean(x) = allequal(x) ? zeros(length(x)) : (x = Float64.(x); x .- mean(skipmissing(x)))
    demeaned_covs = [cov * "_demeaned" for cov in covariates_to_include]
    demeaned = transform(groupby(data, CELL_COLUMNS),
                         :outcome_71X9yTx => (y -> y .- mean(y)) => :y_demeaned,
                         [Symbol(cov) => demean => Symbol(cov * "_demeaned") for cov in covariates_to_include]...)
    cell_means = combine(groupby(data, CELL_COLUMNS),
                         :time_dmG5fpM => first => :date,
                         :outcome_71X9yTx => (y -> mean(y)) => :cell_mean_y,
                         [Symbol(cov) => (x -> mean(skipmissing(x))) => Symbol("cell_mean_" * cov) for cov in covariates_to_include]...)

    # Slopes from a regression of the demeaned outcome on the demeaned covariates, without an
    # intercept as both are already demeaned. A covariate without any within-cell variation in
    # the data of the regression is left out and has slope 0. The least-squares solution comes
    # from a column-pivoted QR decomposition, which also handles collinear covariates.
    function estimate_slopes(df)
        X = Matrix(Float64.(df[:, Symbol.(demeaned_covs)]))
        varies = [any(!iszero, view(X, :, j)) for j in axes(X, 2)]
        slopes = zeros(size(X, 2))
        if any(varies)
            slopes[varies] = qr(X[:, varies], ColumnNorm()) \ Float64.(collect(df.y_demeaned))
        end
        return slopes
    end

    # cell_slopes(state, time) gives the slope vectors that apply to a cell (two for "add")
    if ccc == "hom"
        gamma = estimate_slopes(demeaned)
        cell_slopes = (state, time) -> (gamma,)
    elseif ccc == "state"
        gamma_by_state = Dict(df.state_71X9yTx[1] => estimate_slopes(df) for df in groupby(demeaned, :state_71X9yTx))
        cell_slopes = (state, time) -> (gamma_by_state[state],)
    elseif ccc == "time"
        gamma_by_time = Dict(df.time_71X9yTx[1] => estimate_slopes(df) for df in groupby(demeaned, :time_71X9yTx))
        cell_slopes = (state, time) -> (gamma_by_time[time],)
    elseif ccc == "add"
        gamma_states, gamma_times, state_idx, time_idx = estimate_additive_slopes(demeaned, demeaned_covs)
        cell_slopes = (state, time) -> (gamma_states[:, state_idx[state]], gamma_times[:, time_idx[time]])
    end

    lambdas = Vector{Float64}(undef, nrow(cell_means))
    for (i, row) in enumerate(eachrow(cell_means))
        xmeans = Float64[row[Symbol("cell_mean_" * cov)] for cov in covariates_to_include]
        slopes = cell_slopes(row.state_71X9yTx, row.time_71X9yTx)
        lambdas[i] = foldl((lambda, gamma) -> lambda - dot(gamma, xmeans), slopes; init = row.cell_mean_y)
        fit.slopes[(row.state_71X9yTx, row.date)] = length(slopes) == 1 ? slopes[1] : slopes[1] .+ slopes[2]
    end
    lambda_df = select(cell_means, CELL_COLUMNS)
    lambda_df.lambda = lambdas
    return lambda_df
end

function estimate_additive_slopes(demeaned, demeaned_covs)

    # ccc = "add": state slopes gamma_s and period slopes gamma_t are estimated jointly, from the
    # regression of the demeaned outcome on a sparse design with one block of covariate columns
    # per state followed by one per period. An observation of state s in period t has its
    # demeaned covariates in the block of s and in the block of t. Returns the slopes as
    # matrices with one column per state and per period, and the column of each state and period.
    states = sort(unique(demeaned.state_71X9yTx))
    times = sort(unique(demeaned.time_71X9yTx))
    state_idx = Dict(s => i for (i, s) in enumerate(states))
    time_idx = Dict(t => i for (i, t) in enumerate(times))
    n, ncovs = nrow(demeaned), length(demeaned_covs)

    row_state = [state_idx[s] for s in demeaned.state_71X9yTx]
    row_time = [time_idx[t] for t in demeaned.time_71X9yTx]
    rows = Int[]; cols = Int[]; vals = Float64[]
    for (j, cov) in enumerate(demeaned_covs)
        x = Float64.(demeaned[!, Symbol(cov)])
        append!(rows, 1:n); append!(cols, (row_state .- 1) .* ncovs .+ j);                          append!(vals, x)
        append!(rows, 1:n); append!(cols, length(states) * ncovs .+ (row_time .- 1) .* ncovs .+ j); append!(vals, x)
    end
    Z = dropzeros!(sparse(rows, cols, vals, n, (length(states) + length(times)) * ncovs))

    # Without covariates, or without any within-cell variation in them, there is nothing to
    # solve for (and the sparse solver rejects an empty design): every slope is zero
    gamma = nnz(Z) == 0 ? zeros(size(Z, 2)) : Z \ Float64.(demeaned.y_demeaned)
    gamma_states = reshape(gamma[1:length(states)*ncovs], ncovs, length(states))
    gamma_times = reshape(gamma[length(states)*ncovs+1:end], ncovs, length(times))
    return gamma_states, gamma_times, state_idx, time_idx
end

function prune_covariate_columns(W)

    # Columns of the covariate matrix W to keep in a cell's regression under ccc = "int", in the
    # order of a rank-revealing (column-pivoted) QR decomposition, from most to least
    # informative. Collinear covariates are dropped, and a dropped covariate has slope 0 in that
    # cell. A cell's lambda and slopes must be uniquely identified under "int" (the other ccc
    # options tolerate rank deficiency, as their ATTs are identified anyway).
    # The decomposition is of the covariates centered at their cell means, because the
    # regression has an intercept: covariates can be collinear with the intercept without being
    # collinear with each other (for example the dummies of a factor whose reference category
    # does not occur in the cell add up to one).
    # This is a choice: the Python implementation didintpy keeps every covariate and uses the
    # minimum-norm least-squares solution instead. The two agree when a covariate is constant
    # within the cell, but can give different slopes, and so different ATT(s,t)'s, for other
    # collinearities.

    Wcol = size(W, 2)
    Wrow = size(W, 1)

    Wqr = qr(W .- mean(W, dims = 1), ColumnNorm())
    m = min(Wcol, Wrow)
    # A column counts as collinear when what is left of it after removing the earlier columns
    # is within rounding error of zero. Centering a column leaves rounding errors in proportion
    # to its uncentered values, so the tolerance is relative to the largest uncentered column
    # norm, with the usual factor of the larger matrix dimension. (A tighter tolerance lets the
    # decision depend on the processor, as the rounding errors differ between them.)
    tol = eps(Float64) * max(Wcol, Wrow) * maximum(norm, eachcol(W))
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


###########################################################################################
# Section 4: The first stage with covariances, for edge-case standard errors              #
#                                                                                         #
# An ATT(s,t) with a single control c is a combination of the lambdas of four cells,      #
# (s, r1), (s, t), (c, r1) and (c, t), and its edge-case standard error needs their       #
# covariances. Those come from regressions of the outcome on cell dummies and the         #
# covariates, in blocks that follow the ccc: one regression per cell ("int"), per state   #
# ("state"), per period ("time") or a single one ("hom").                                 #
###########################################################################################

function find_edgecase_cells(data_working, treated_states, treatment_times, time_to_index, match_to_these_dates, use_pre_controls)

    # The cells (state, period) behind every ATT(s,t) that has exactly one valid control,
    # whatever the agg: for treated state s, period t and control c, the cells (s, r1), (s, t),
    # (c, r1) and (c, t). Returns nothing if there is no such ATT(s,t).
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
    return isempty(out_states) ? nothing : unique(DataFrame(state = out_states, time = out_times))
end

function edgecase_rows(data_working, ccc, edgecase_cells)

    # Observations that go through first_stage_edgecase: those of the edge-case cells, and of
    # every cell that shares slopes with one of them (all cells under "hom", the cells of the
    # same state under "state", of the same period under "time")
    if ccc == "hom"
        return trues(nrow(data_working))
    elseif ccc == "state"
        return in.(data_working.state_71X9yTx, Ref(unique(edgecase_cells.state)))
    elseif ccc == "time"
        return in.(data_working.time_dmG5fpM, Ref(unique(edgecase_cells.time)))
    elseif ccc == "int"
        cells = Set(zip(edgecase_cells.state, edgecase_cells.time))
        return [cell in cells for cell in zip(data_working.state_71X9yTx, data_working.time_dmG5fpM)]
    end
end

function first_stage_edgecase(data_working, ccc, covariates_to_include, lambda_df, vcov_lambda, hc, fit::CovariateFit)

    # Lambdas of the cells in `data_working` with their HC covariance matrix, from regressions
    # of the outcome on cell dummies and the covariates. Fills the rows of `lambda_df` and the
    # blocks of `vcov_lambda` that belong to those cells, and records the cells' slopes in `fit`.
    # Cells in different blocks are estimated separately, so their covariances stay zero.
    n_lambda = nrow(lambda_df)

    # Position of each cell in lambda_df and vcov_lambda (lambda_df.time is still a string here,
    # as is data_working.time_71X9yTx)
    cell_id_map = Dict((lambda_df.state[i], lambda_df.time[i]) =>
                       lambda_df.lambda_index[i] for i in 1:n_lambda)

    # Regression of one block. Returns the covariance matrix of its lambdas, their positions
    # and their values.
    function block_regression(sub_df; int = false)
        n_sub = nrow(sub_df)

        cells_seen = sort(unique([(sub_df.state_71X9yTx[i], sub_df.time_71X9yTx[i])
                                  for i in 1:n_sub]))
        local_idx  = Dict(c => i for (i, c) in enumerate(cells_seen))
        n_cells    = length(cells_seen)
        global_idx = [cell_id_map[c] for c in cells_seen]

        # Cell dummies
        D_block = zeros(n_sub, n_cells)
        for i in 1:n_sub
            ci = local_idx[(sub_df.state_71X9yTx[i], sub_df.time_71X9yTx[i])]
            D_block[i, ci] = 1.0
        end

        # Covariates that vary within at least one cell of the block (the others are absorbed
        # by the cell dummies)
        active = filter(covariates_to_include) do c
            stds = combine(groupby(sub_df, CELL_COLUMNS)) do df
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
        # Without residual degrees of freedom, or without full rank, the covariance matrix
        # (and so the edge-case SE) cannot be computed
        estimable = (size(Z, 1) > size(Z, 2)) && (rank(Z) == size(Z, 2))
        if estimable
            resid = y .- Z * β
            V     = compute_hc_covariance(Z, resid, hc)
        else
            V = fill(NaN, n_cells, n_cells)
        end

        # Record the slopes of every cell in the block: common to the block for "hom", "state"
        # and "time", the cell's own for "int"
        kept_positions = Int.(indexin(kept_covs, covariates_to_include))
        gamma = zeros(length(covariates_to_include))
        gamma[kept_positions] = β[n_cells+1:end]
        for cell in cells_seen
            fit.slopes[(cell[1], Date(cell[2]))] = gamma
        end

        # Under "int" with covariates the block is one cell. Keep the HC covariance of its
        # (lambda, gamma), with zero rows and columns for dropped covariates, and its
        # observations: edge_case_se_integrated needs them.
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

    # Write a block's lambdas and their covariance matrix into lambda_df and vcov_lambda
    function place_block!(V_block, idx, lambda_vals)
        for j in 1:length(idx), i in 1:length(idx)
            vcov_lambda[idx[i], idx[j]] = V_block[i, j]
        end
        for i in 1:length(idx)
            lambda_df.lambda[idx[i]] = lambda_vals[i]
        end
    end

    if ccc == "int"
        for grp in groupby(data_working, CELL_COLUMNS)
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
