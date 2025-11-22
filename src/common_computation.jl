
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

function run_fixed_effects_model(data_copy, formula, ccc, covariates; common_adoption = false, staggered_adoption = true)

    # Call GC.gc() before running big FixedEffectsModels regression
    GC.gc()

    # Determine if corner case (i.e. only two states common adoption)
    cornercase = false
    if common_adoption && length(unique(data_copy.state_71X9yTx)) == 2
        cornercase = true
    end 

    # Run FixedEffectsModels regression (stage 1)
    if cornercase
        stage1 = reg(data_copy, formula, Vcov.robust(), contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
                     save = false)
    else
        stage1 = reg(data_copy, formula, contrasts = Dict(:state_71X9yTx => DummyCoding(), :time_71X9yTx => DummyCoding()),
                     save = false)
    end

    # Recover lambdas
    state_time = split.(replace.(coefnames(stage1), "state_time: " => ""), "0IQR7q6Wei7Ejp4e")
    coefs = coef(stage1)
    if ccc == "hom" && !isnothing(covariates)
        keep_mask = map(x -> !(x[1] in covariates), state_time)
        state_time = state_time[keep_mask]
        coefs = coefs[keep_mask]
    end 
    
    lambda_df = DataFrame(state = first.(state_time), time = last.(state_time))
    lambda_df.lambda = coefs

    # Handle edge case where FE structure drops some state-time coefficients,
    # that is, effect is entirely absorbed by fixed effects thus the estimate from
    # the state_time dummy should be exactly zero
    expected_df = unique(data_copy[:, [:state_71X9yTx, :time_71X9yTx]])
    expected_df = rename(expected_df, :state_71X9yTx => :state, :time_71X9yTx => :time)
    lambda_df = leftjoin(expected_df, lambda_df, on = [:state, :time])
    lambda_df.lambda = coalesce.(lambda_df.lambda, 0.0)

    # Parse the treatment_times and the time column to dates for staggered adoption (left as true for didint_plot)
    if staggered_adoption
        lambda_df.time = Date.(lambda_df.time)
    end    

    # Note the ccc for this loop
    lambda_df.ccc .= ccc

    # Common adoption and only 2 states can still compue SE manually
    if cornercase
        vcov_matrix = vcov(stage1)[1:4,1:4]

        treat_post_var = vcov_matrix[1,1]
        treat_pre_var = vcov_matrix[2,2]
        treat_cov = vcov_matrix[2,1]
        control_post_var = vcov_matrix[3,3]
        control_pre_var = vcov_matrix[4,4]
        control_cov = vcov_matrix[4,3]

        cornercase_se = sqrt(treat_post_var + treat_pre_var + control_post_var + control_pre_var - 2*treat_cov - 2*control_cov)
    else
        cornercase_se = nothing
    end 

    return lambda_df, cornercase, cornercase_se

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
    
    Ω = Diagonal(omega_diag)
    
    # Sandwich estimator: (X'X)⁻¹ X'ΩX (X'X)⁻¹
    return XXinv * (X' * Ω * X) * XXinv
end