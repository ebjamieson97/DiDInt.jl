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


###########################################################################################
# Section 2: First-stage results beyond the lambdas (CovariateFit)                        #
###########################################################################################

# First-stage estimates other than the cell intercepts lambda, keyed by cell (state, period).
# `slopes` holds the covariate slopes gamma_{s,t} implied by the ccc, in the order of
# `covariates` (0 for a covariate dropped in that cell), and `xbar` the cell's covariate
# means; the integrated long differences need both (Karim & Webb 2025, eqs. 46-47). For the
# cells of an edge case under ccc = "int", `cell_vcov` holds the HC covariance of (lambda,
# gamma) from the cell's regression, which edge_case_se_integrated needs.
# The last two fields describe what the data leave undetermined. Under ccc = "int",
# `null_space` holds, for every cell whose regression drops covariates, the directions in
# which that cell's slopes can move without changing its fit (one column per dropped
# covariate; see null_space_basis), and `covariate_scale` the root mean square of every
# covariate over the data, the unit in which differences between covariate means are judged
# (see dropped_covariate_effects). Under ccc = "add", `lambda_shifts` holds, for every
# direction in which the additive slopes are undetermined and the lambdas change, the change
# in every cell's lambda and the covariates involved, and `additive_rank` the rank of the
# first-stage design and the rank expected (see additive_lambda_shifts!); both are empty
# when the lambdas are uniquely determined.
const LambdaShift = NamedTuple{(:shift, :covariates), Tuple{Dict{Tuple{String, Date}, Float64}, Vector{String}}}

struct CovariateFit
    covariates::Vector{String}
    slopes::Dict{Tuple{String, Date}, Vector{Float64}}
    xbar::Dict{Tuple{String, Date}, Vector{Float64}}
    cell_vcov::Dict{Tuple{String, Date}, Matrix{Float64}}
    null_space::Dict{Tuple{String, Date}, Matrix{Float64}}
    covariate_scale::Vector{Float64}
    lambda_shifts::Vector{LambdaShift}
    additive_rank::Vector{Int}
end

function CovariateFit(covariates::Vector{String})
    return CovariateFit(covariates, Dict{Tuple{String, Date}, Vector{Float64}}(), Dict{Tuple{String, Date}, Vector{Float64}}(),
                        Dict{Tuple{String, Date}, Matrix{Float64}}(), Dict{Tuple{String, Date}, Matrix{Float64}}(),
                        Float64[], LambdaShift[], Int[])
end

# Relative tolerance of the checks of what the data leave undetermined (dropped_covariate_effects
# under ccc = "int", additive_lambda_shifts! and additive_att_shifts under ccc = "add"):
# a change smaller than this share of the size of the quantities involved counts as zero
const IDENTIFICATION_RTOL = 1e-8


###########################################################################################
# Section 3: The first stage: lambda_{s,t} for every state and period                     #
###########################################################################################

# Columns that identify a cell: the state and the period (as a string)
const CELL_COLUMNS = [:state_71X9yTx, :time_71X9yTx]

function first_stage(data_working, ccc, covariates_to_include, hc, edgecase;
                     time_to_index = nothing, treatment_times = nothing, match_to_these_dates = nothing,
                     treated_states = nothing, use_pre_controls = nothing, warn = true)

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
    # `warn = false` silences the warnings about collinear covariates (used when the first
    # stage is re-estimated for the jackknife, after the full sample has been reported on).

    # Covariate means of every cell
    fit = CovariateFit(String[string(cov) for cov in covariates_to_include])
    for cell in groupby(data_working, [:state_71X9yTx, :time_dmG5fpM])
        fit.xbar[(cell.state_71X9yTx[1], cell.time_dmG5fpM[1])] = Float64[mean(cell[!, Symbol(cov)]) for cov in covariates_to_include]
    end

    # With covariates: the outcome and the covariates demeaned within cells, and the covariates
    # that each slope block keeps. The kept covariates are determined once, here, and used by
    # the regular and the edge-case procedure alike, so that both estimate the same model.
    kept = Dict{Tuple{String, String}, Vector{Int}}()
    if ccc != "none" && !isempty(covariates_to_include)
        data_working = with_demeaned_columns(data_working, fit.covariates)
        kept, null_space = kept_covariates(data_working, ccc, fit.covariates; warn = warn)
        # Under "int": what the cells that drop covariates leave undetermined, for the check
        # of which ATT(s,t)'s that affects (warn_dropped_covariates_int)
        for ((state, time), basis) in null_space
            fit.null_space[(state, Date(time))] = basis
        end
        append!(fit.covariate_scale, Float64[sqrt(mean(abs2, data_working[!, Symbol(cov)])) for cov in covariates_to_include])
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
                                                     skeleton, vcov_lambda, hc, fit, kept)
        data_working = data_working[.!in_edgecase, :]
    end

    # Regular procedure
    if !has_edgecase || !isempty(data_working)
        lambda_df = cell_lambdas!(fit, data_working, ccc, covariates_to_include, kept; warn = warn)
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

function cell_lambdas!(fit::CovariateFit, data, ccc, covariates_to_include, kept; warn = true)

    # lambda of every cell of `data` by the regular procedure, as a DataFrame with the cell
    # columns and `lambda`. Also records every cell's covariate slopes in `fit`. `kept` gives
    # the covariates that each slope block keeps (see kept_covariates).
    if ccc == "int"
        return cell_lambdas_int!(fit, data, covariates_to_include, kept)
    elseif ccc == "none" || isempty(covariates_to_include)
        # No covariate adjustment (ccc = "none" is used by didint_plot): lambda is the cell
        # mean of the outcome
        return combine(groupby(data, CELL_COLUMNS)) do df
            fit.slopes[(df.state_71X9yTx[1], df.time_dmG5fpM[1])] = zeros(length(covariates_to_include))
            DataFrame(lambda = [mean(Float64.(collect(df.outcome_71X9yTx)))])
        end
    else
        return cell_lambdas_fwl!(fit, data, ccc, covariates_to_include, kept; warn = warn)
    end
end

function cell_lambdas_int!(fit::CovariateFit, data, covariates_to_include, kept)

    # ccc = "int": the slopes are specific to the cell, so each cell has its own regression of
    # the outcome on an intercept and the covariates that the cell keeps, and lambda is the
    # intercept. A covariate that does not vary in the cell, or that is collinear with
    # covariates listed before it, is not kept and has slope 0 there (see kept_covariates).
    return combine(groupby(data, CELL_COLUMNS)) do df
        y = Float64.(collect(df.outcome_71X9yTx))
        gamma = zeros(length(covariates_to_include))
        cell_kept = get(kept, (df.state_71X9yTx[1], df.time_71X9yTx[1]), Int[])
        if isempty(cell_kept)
            X = ones(nrow(df), 1)
            beta = X \ y
        else
            X = hcat(ones(nrow(df)), Matrix{Float64}(df[:, covariates_to_include[cell_kept]]))
            beta = X \ y
            gamma[cell_kept] = beta[2:end]
        end
        fit.slopes[(df.state_71X9yTx[1], df.time_dmG5fpM[1])] = gamma
        DataFrame(lambda = [beta[1]])
    end
end

function cell_lambdas_fwl!(fit::CovariateFit, data, ccc, covariates_to_include, kept; warn = true)

    # ccc = "hom", "state", "time" or "add": the slopes are shared across cells. By the
    # Frisch-Waugh-Lovell theorem they can be estimated from the outcome and the covariates
    # demeaned within cells (which removes the cell dummies), and then
    #   lambda_{s,t} = ybar_{s,t} - gamma_{s,t}' xbar_{s,t}
    # with gamma_{s,t} the slopes that apply to the cell: gamma ("hom"), gamma_s ("state"),
    # gamma_t ("time") or gamma_s + gamma_t ("add"). `data` holds the demeaned outcome and
    # covariates (see with_demeaned_columns).
    demeaned_covs = demeaned_name.(covariates_to_include)
    cell_means = combine(groupby(data, CELL_COLUMNS),
                         :time_dmG5fpM => first => :date,
                         :outcome_71X9yTx => (y -> mean(y)) => :cell_mean_y,
                         [Symbol(cov) => (x -> mean(skipmissing(x))) => Symbol("cell_mean_" * cov) for cov in covariates_to_include]...)

    # Slopes of one slope block ("hom", "state", "time"), from the regression of the demeaned
    # outcome on the block's kept covariates, demeaned; there is no intercept, as both sides
    # are demeaned. The kept covariates are linearly independent, so the least-squares
    # solution is unique. A covariate that is not kept has slope 0.
    function estimate_slopes(df)
        slopes = zeros(length(covariates_to_include))
        block_kept = kept[slope_block(ccc, df.state_71X9yTx[1], df.time_71X9yTx[1])]
        if !isempty(block_kept)
            slopes[block_kept] = qr(Matrix{Float64}(df[:, demeaned_covs[block_kept]])) \ df.outcome_demeaned_71X9yTx
        end
        return slopes
    end

    # cell_slopes(state, time) gives the slope vectors that apply to a cell (two for "add")
    if ccc == "hom"
        gamma = estimate_slopes(data)
        cell_slopes = (state, time) -> (gamma,)
    elseif ccc == "state"
        gamma_by_state = Dict(df.state_71X9yTx[1] => estimate_slopes(df) for df in groupby(data, :state_71X9yTx))
        cell_slopes = (state, time) -> (gamma_by_state[state],)
    elseif ccc == "time"
        gamma_by_time = Dict(df.time_71X9yTx[1] => estimate_slopes(df) for df in groupby(data, :time_71X9yTx))
        cell_slopes = (state, time) -> (gamma_by_time[time],)
    elseif ccc == "add"
        additive = estimate_additive_slopes(data, demeaned_covs; null_space = warn)
        cell_slopes = (state, time) -> (additive.gamma_states[:, additive.state_idx[state]],
                                        additive.gamma_times[:, additive.time_idx[time]])
        # What the data leave undetermined beyond the split of the slopes between states and periods
        warn && additive_lambda_shifts!(fit, additive, cell_means, covariates_to_include)
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

function estimate_additive_slopes(demeaned, demeaned_covs; null_space = true)

    # ccc = "add": state slopes gamma_s and period slopes gamma_t are estimated jointly, from the
    # regression of the demeaned outcome on a sparse design with one block of covariate columns
    # per state followed by one per period. An observation of state s in period t has its
    # demeaned covariates in the block of s and in the block of t.
    # This design is always rank deficient: adding a constant to every state's slope of a
    # covariate and subtracting it from every period's slope leaves the slope of every cell,
    # gamma_s + gamma_t, unchanged. Every least-squares solution therefore gives the same
    # lambdas, and no covariate is dropped. The system is solved by a sparse QR decomposition
    # (SuiteSparseQR, which is what `\` uses for a rectangular sparse matrix). For a rank
    # deficient design it returns a basic solution: the coefficients of the columns that it
    # finds linearly dependent on others are set to zero. It is not the minimum-norm solution.
    # Returns
    #   gamma_states, gamma_times  the slopes, as matrices with one column per state and period
    #   state_idx, time_idx        the column of each state and period
    #   null_space                 a basis of the directions in which the solution can move
    #                              without changing the fit, one per column (only computed with
    #                              `null_space = true`; additive_lambda_shifts! uses it)
    #   rank, expected_rank        the rank of the design and the rank that the split between
    #                              states and periods alone would leave
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
    n_columns = (length(states) + length(times)) * ncovs
    Z = dropzeros!(sparse(rows, cols, vals, n, n_columns))

    # Without any within-cell variation in the covariates there is nothing to solve for (and
    # the sparse solver rejects an empty design): every slope is zero, and every direction
    # leaves the fit unchanged
    if nnz(Z) == 0
        gamma = zeros(n_columns)
        rank_found = 0
        directions = null_space ? Matrix{Float64}(I, n_columns, n_columns) : zeros(n_columns, 0)
    else
        Zqr = qr(Z)
        gamma = Zqr \ Float64.(demeaned.outcome_demeaned_71X9yTx)
        rank_found = rank(Zqr)
        directions = null_space ? sparse_qr_null_space(Zqr, n_columns) : zeros(n_columns, 0)
    end

    # The split between states and periods costs one rank per covariate and per connected
    # component of the graph that links each state to the periods it is observed in (one
    # component, unless the panel falls apart into groups of states observed in disjoint sets
    # of periods)
    components = null_space ? count_panel_components(demeaned.state_71X9yTx, demeaned.time_71X9yTx) : 1
    expected_rank = ncovs * (length(states) + length(times) - components)

    return (gamma_states = reshape(gamma[1:length(states)*ncovs], ncovs, length(states)),
            gamma_times = reshape(gamma[length(states)*ncovs+1:end], ncovs, length(times)),
            state_idx = state_idx, time_idx = time_idx, null_space = directions,
            rank = rank_found, expected_rank = expected_rank)
end

function sparse_qr_null_space(Zqr, n_columns)

    # A basis of the null space of the matrix Z behind the sparse QR decomposition Zqr, one
    # vector per column of the matrix returned. With r the rank that the decomposition found,
    # Z[:, pcol] = Q R with R = [R11 R12], R11 upper triangular of size r: the columns pcol[1:r]
    # are the ones the decomposition treats as independent (it calls a column dependent when
    # its norm, after removing the earlier columns, is at most 20 * (rows + columns) * eps *
    # (largest column norm), the default of SuiteSparseQR), and every other column is a
    # combination of them. For each dependent column pcol[r + j] the basis vector has a one in
    # that position and -inv(R11) R12[:, j] in the positions pcol[1:r]. This is the null space
    # that goes with the solution that `\` returns, which sets the dependent columns'
    # coefficients to zero.
    r = rank(Zqr)
    basis = zeros(n_columns, n_columns - r)
    r == n_columns && return basis
    if r > 0
        basis[Zqr.pcol[1:r], :] = -(UpperTriangular(Matrix(Zqr.R[1:r, 1:r])) \ Matrix(Zqr.R[1:r, r+1:end]))
    end
    for (j, column) in enumerate(Zqr.pcol[r+1:end])
        basis[column, j] = 1.0
    end
    return basis
end

function additive_lambda_shifts!(fit::CovariateFit, additive, cell_means, covariates_to_include)

    # ccc = "add": find the directions in which the least-squares solution is undetermined and
    # that change the lambdas, and record them in `fit`.
    # Moving the solution along a null-space direction d changes the slopes of cell (s,t) by
    # d_s + d_t and its lambda by -(d_s + d_t)' xbar_{s,t}. For the split between states and
    # periods d_s + d_t = 0 in every cell, so the lambdas do not change: that direction is
    # harmless and is skipped. Any other direction changes some lambdas; whether it changes
    # an ATT is checked once the ATT(s,t)'s are known (warn_additive_identification).
    # A lambda counts as changed when the change exceeds IDENTIFICATION_RTOL times
    # sum_i w_i scale_i, with w_i the largest entry of d for covariate i and scale_i the root
    # mean square of covariate i over the data. The covariates recorded as involved are those
    # whose own term of the change, (d_s + d_t)_i xbar_i, exceeds the same threshold in some cell.
    n_states, p = length(additive.state_idx), length(covariates_to_include)
    # Covariate means of every cell (one column per cell), and each cell's state and period
    xbar = permutedims(Matrix{Float64}(cell_means[:, [Symbol("cell_mean_" * cov) for cov in covariates_to_include]]))
    state_of_cell = [additive.state_idx[s] for s in cell_means.state_71X9yTx]
    time_of_cell = [additive.time_idx[t] for t in cell_means.time_71X9yTx]
    for direction in eachcol(additive.null_space)
        d_states = reshape(direction[1:n_states*p], p, n_states)
        d_times = reshape(direction[n_states*p+1:end], p, :)
        weight = vec(max.(maximum(abs, d_states, dims = 2), maximum(abs, d_times, dims = 2)))

        # Each covariate's term of the change in every cell's lambda, and the changes themselves
        terms = (d_states[:, state_of_cell] .+ d_times[:, time_of_cell]) .* xbar
        shifts = -vec(sum(terms, dims = 1))
        threshold = IDENTIFICATION_RTOL * dot(weight, fit.covariate_scale)
        if maximum(abs, shifts) > threshold
            involved = covariates_to_include[vec(maximum(abs, terms, dims = 2)) .> threshold]
            shift = Dict(zip(zip(cell_means.state_71X9yTx, cell_means.date), shifts))
            push!(fit.lambda_shifts, (shift = shift, covariates = String.(involved)))
        end
    end
    if !isempty(fit.lambda_shifts)
        append!(fit.additive_rank, [additive.rank, additive.expected_rank])
    end
    return nothing
end

function count_panel_components(states, times)

    # Number of connected components of the graph whose nodes are the states and the periods,
    # with a link between a state and every period it is observed in (union-find)
    parent = Dict{Any, Any}()
    function root(node)
        while parent[node] != node
            node = parent[node]
        end
        return node
    end
    for (state, time) in unique(zip(states, times))
        a = root(get!(parent, (:state, state), (:state, state)))
        b = root(get!(parent, (:time, time), (:time, time)))
        parent[a] = b
    end
    return count(node -> root(node) == node, collect(keys(parent)))
end

# Relative tolerance of the collinearity rule in prune_covariate_columns: a covariate is
# dropped when less than this share of it (in norm) is left after removing what the
# covariates before it explain. 1e-7 is also the tolerance of R's lm().
const COLLINEARITY_RTOL = 1e-7

# Name of the column that holds a covariate demeaned within cells (the suffix keeps it from
# clashing with a column of the user's data)
demeaned_name(covariate) = covariate * "_demeaned_71X9yTx"

function with_demeaned_columns(data, covariates)

    # `data` with two kinds of added columns: the outcome demeaned within cells (outcome_demeaned_71X9yTx)
    # and every covariate demeaned within cells (named by demeaned_name). The existing columns
    # are shared with `data`, not copied.
    # A covariate that does not vary within a cell (for example one measured at the state-period
    # level) is absorbed by the cell dummies. Its demeaned values are set to exactly zero:
    # subtracting the mean of identical values can leave rounding noise of about 1e-17, and a
    # regression on a column of pure rounding noise returns an arbitrary slope.
    cell_rows = [parentindices(cell)[1] for cell in groupby(data, CELL_COLUMNS)]
    function demean(x, exact_zero)
        demeaned = zeros(length(x))
        for rows in cell_rows
            values = view(x, rows)
            exact_zero && allequal(values) && continue
            demeaned[rows] .= values .- mean(values)
        end
        return demeaned
    end
    demeaned_data = copy(data, copycols = false)
    demeaned_data.outcome_demeaned_71X9yTx = demean(data.outcome_71X9yTx, false)
    for cov in covariates
        demeaned_data[!, demeaned_name(cov)] = demean(data[!, cov], true)
    end
    return demeaned_data
end

# The slope block of a cell: the observations that share one set of covariate slopes are all
# the data ("hom"), a state ("state"), a period ("time") or a cell ("int")
function slope_block(ccc, state, time)
    return ccc == "hom" ? ("", "") : ccc == "state" ? (state, "") : ccc == "time" ? ("", time) : (state, time)
end

function kept_covariates(data, ccc, covariates; warn = true)

    # The covariates that each slope block keeps under ccc = "hom", "state", "time" and "int",
    # as positions in `covariates`, keyed by slope_block. `data` holds the covariates demeaned
    # within cells (see with_demeaned_columns).
    # When the demeaned covariates of a block are linearly dependent, its slopes, and with them
    # the lambdas, are not determined by the data. The ATTs would only be unaffected if the
    # lambdas could change by amounts that are additive in state and period, which cancel in
    # the long differences; collinear covariates do not have that form in general. Redundant
    # covariates are therefore dropped, by the rule of prune_covariate_columns, and have slope
    # 0 in that block.
    # Under "hom", "state" and "time" a warning names every covariate that was dropped although
    # it varies within cells.
    # Under "int" the block is a cell and demeaning within the cell is centering at the cell
    # means, as needed for a regression with an intercept: covariates can be collinear with the
    # intercept without being collinear with each other (the dummies of a factor whose
    # reference category does not occur in the cell add up to one). No warning is given here:
    # whether a dropped covariate affects an ATT(s,t) depends on where the cell's slopes are
    # evaluated, which is checked once the ATT(s,t)'s are known (warn_dropped_covariates_int).
    # For that check the second value returned holds, for every cell that drops covariates,
    # a basis of the directions in which its slopes are undetermined (see null_space_basis).
    # ccc = "add" keeps every covariate (see estimate_additive_slopes).
    kept = Dict{Tuple{String, String}, Vector{Int}}()
    null_space = Dict{Tuple{String, String}, Matrix{Float64}}()
    ccc in ("hom", "state", "time", "int") || return kept, null_space

    X = Matrix{Float64}(data[:, demeaned_name.(covariates)])
    block_columns = ccc == "state" ? [:state_71X9yTx] : ccc == "time" ? [:time_71X9yTx] : CELL_COLUMNS
    blocks = ccc == "hom" ? [collect(1:nrow(data))] : [parentindices(block)[1] for block in groupby(data, block_columns)]
    collinear_blocks = [Tuple{String, String}[] for _ in covariates]
    for rows in blocks
        block = slope_block(ccc, data.state_71X9yTx[rows[1]], data.time_71X9yTx[rows[1]])
        kept[block], collinear = prune_covariate_columns(view(X, rows, :))
        for j in collinear
            push!(collinear_blocks[j], block)
        end
        if ccc == "int" && length(kept[block]) < length(covariates)
            null_space[block] = null_space_basis(view(X, rows, :), kept[block])
        end
    end
    warn && ccc != "int" && warn_collinear_covariates(covariates, collinear_blocks, ccc)
    return kept, null_space
end

function prune_covariate_columns(X)

    # Which columns of X (the covariates of one slope block, demeaned within cells) to keep,
    # going through them in order, which is the order in which the user listed the covariates
    # (with the dummies of a factor in its place). A column is kept if it is not, up to the
    # relative tolerance COLLINEARITY_RTOL, a linear combination of the columns kept before it:
    # the norm of what is left of it after projecting on those columns must exceed the
    # tolerance times its own norm. This is a Gram-Schmidt orthogonalization without column
    # pivoting, so the decision does not depend on the units of the covariates, only on their
    # order. (Each column is orthogonalized twice, which keeps the result accurate to rounding
    # error when columns are nearly dependent.)
    # Returns the kept columns and the columns dropped as collinear. A column of zeros (a
    # covariate without any within-cell variation in the block, which the cell dummies absorb)
    # is in neither.
    kept = Int[]
    collinear = Int[]
    basis = Vector{Float64}[]
    for j in axes(X, 2)
        column = X[:, j]
        column_norm = norm(column)
        column_norm == 0 && continue
        for _ in 1:2, q in basis
            column .-= dot(q, column) .* q
        end
        remainder = norm(column)
        if remainder > COLLINEARITY_RTOL * column_norm
            push!(kept, j)
            push!(basis, column ./ remainder)
        else
            push!(collinear, j)
        end
    end
    return kept, collinear
end

function null_space_basis(X, kept)

    # A basis of the null space of X (the covariates of one cell, centered at the cell means),
    # consistent with the columns that prune_covariate_columns keeps: one basis vector per
    # dropped column j, as the columns of the matrix returned. If column j is zero (the
    # covariate does not vary in the cell), the vector is the unit vector e_j. Otherwise
    # column j is a linear combination of the kept columns before it, x_j = sum_i c_i x_i, and
    # the vector is e_j - sum_i c_i e_i. Adding any multiple of a basis vector to the cell's
    # slopes leaves its fitted values on its own data unchanged.
    p = size(X, 2)
    dropped = setdiff(1:p, kept)
    basis = zeros(p, length(dropped))
    for (column, j) in enumerate(dropped)
        basis[j, column] = 1.0
        kept_before = filter(<(j), kept)
        if !isempty(kept_before) && !iszero(view(X, :, j))
            basis[kept_before, column] = -(qr(X[:, kept_before]) \ X[:, j])
        end
    end
    return basis
end

function warn_collinear_covariates(covariates, collinear_blocks, ccc)

    # One warning that names every covariate dropped as collinear and where it was dropped
    # (ccc = "hom", "state" and "time")
    any(!isempty, collinear_blocks) || return nothing
    where_label(block) = ccc == "state" ? "state $(block[1])" : "period $(block[2])"
    descriptions = String[]
    for (cov, blocks) in zip(covariates, collinear_blocks)
        isempty(blocks) && continue
        if ccc == "hom"
            push!(descriptions, cov)
        else
            shown = join(where_label.(first(blocks, 3)), ", ")
            more = length(blocks) > 3 ? " and $(length(blocks) - 3) more" : ""
            push!(descriptions, "$cov ($shown$more)")
        end
    end
    @warn "Under ccc = \"$ccc\", covariates that are collinear with covariates listed before them were dropped" *
          (ccc == "hom" ? "" : " where their slopes are estimated") * ": " * join(descriptions, "; ") *
          ". The estimates depend on which of the collinear covariates is dropped, which is decided by the order of 'covariates'."
    return nothing
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

function first_stage_edgecase(data_working, ccc, covariates_to_include, lambda_df, vcov_lambda, hc, fit::CovariateFit, kept)

    # Lambdas of the cells in `data_working` with their HC covariance matrix, from regressions
    # of the outcome on cell dummies and the covariates. Fills the rows of `lambda_df` and the
    # blocks of `vcov_lambda` that belong to those cells, and records the cells' slopes in `fit`.
    # Cells in different blocks are estimated separately, so their covariances stay zero.
    # The blocks are the slope blocks of the ccc, and each regression includes the covariates
    # that its block keeps (`kept`, see kept_covariates): the same covariates as the regular
    # procedure, so both give the same lambdas (the Frisch-Waugh-Lovell theorem).
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

        # Covariates in the regression: the ones that the block keeps. The others either do not
        # vary within any cell of the block (the cell dummies absorb them) or are collinear
        # with covariates listed before them.
        kept_positions = get(kept, slope_block(ccc, sub_df.state_71X9yTx[1], sub_df.time_71X9yTx[1]), Int[])
        Z = isempty(kept_positions) ? D_block :
            hcat(D_block, Matrix{Float64}(sub_df[:, covariates_to_include[kept_positions]]))

        y     = Float64.(sub_df.outcome_71X9yTx)
        β = Z \ y
        # Without residual degrees of freedom the covariance matrix (and so the edge-case SE)
        # cannot be computed. The design has full column rank, as the kept covariates are
        # linearly independent after removing the cell dummies; the rank is still checked, as
        # a guard against inverting a singular matrix.
        estimable = (size(Z, 1) > size(Z, 2)) && (rank(Z) == size(Z, 2))
        if estimable
            resid = y .- Z * β
            V     = compute_hc_covariance(Z, resid, hc)
        else
            V = fill(NaN, n_cells, n_cells)
        end

        # Record the slopes of every cell in the block: common to the block for "hom", "state"
        # and "time", the cell's own for "int"
        gamma = zeros(length(covariates_to_include))
        gamma[kept_positions] = β[n_cells+1:end]
        for cell in cells_seen
            fit.slopes[(cell[1], Date(cell[2]))] = gamma
        end

        # Under "int" with covariates the block is one cell. Keep the HC covariance of its
        # (lambda, gamma), with zero rows and columns for dropped covariates:
        # edge_case_se_integrated needs it.
        if int && !isempty(covariates_to_include)
            n_coef = 1 + length(covariates_to_include)
            positions = vcat(1, 1 .+ kept_positions)
            cell_vcov = estimable ? zeros(n_coef, n_coef) : fill(NaN, n_coef, n_coef)
            if estimable
                cell_vcov[positions, positions] = V
            end
            cell_key = (cells_seen[1][1], Date(cells_seen[1][2]))
            fit.cell_vcov[cell_key] = cell_vcov
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
