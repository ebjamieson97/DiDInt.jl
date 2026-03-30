
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

function run_fixed_effects_model(data_copy, formula, ccc, covariates, covariates_to_include; common_adoption = false, staggered_adoption = true, recover = true)

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