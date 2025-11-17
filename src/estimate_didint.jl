"""
    didint(outcome::AbstractString, state::AbstractString,
           time::AbstractString, data::DataFrame,
           treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number},
           treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date};
           date_format::Union{AbstractString, Nothing} = nothing,
           covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
           ccc::AbstractString = "int", agg::AbstractString = "state",
           ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
           freq::Union{AbstractString, Nothing} = nothing, freq_multiplier::Number = 1,
           autoadjust::Bool = false, nperm::Number = 1000, verbose::Bool = true,
           seed::Number = rand(1:1000000))

The `didint()` function estimates the average effect of treatment on the treated (ATT)
while accounting for covariates that may vary by state, time, or by both state and time simultaneously.

# Details
The arguments `treated_states` and `treated_times` should be entered such that the first element in
`treated_states` refers to the state treated at the date entered as the first element in `treated_times`,
the second element in `treated_states` refers to the state treated at the date entered as the second element
in `treated_times`, and so on. 

# Parameters
- `outcome::AbstractString` 
    Input the name of the column which identifies the outcome of interest.
- `state::AbstractString` 
    Input the name of the column which identifies the state membership of the observation.
- `time::AbstractString`
    Input the name of the column which identifies the date of the observation.
- `data::DataFrame`
    The DataFrame to be used for the analysis.
- `treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number}`
    A vector of strings (or a single string) noting the treated state(s).
- `treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date}`
    A vector (or single entry) denoting the associated treatment times of the 'treated_states'.
- `covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing}` 
    A vector of covariates entered as strings (or a single covariate string),
    or, `nothing` (default).
- `ccc::AbstractString = "int"`
    Specify which version of DID-INT should be used.
    Options are: `"hom"`, `"time"`, `"state"`, "`add`", and `"int"` (default). 
- `agg::AbstractString = "cohort"` 
    Enter the weighting method as a string.
    Options are: `"cohort"` (default), `"simple"`, `"state"`, `"sgt"`, `"none"`.
- `ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing`
    A dictionary specifying which category in a categorical variable should be used
    as the reference (baseline) category.
- `freq::Union{AbstractString, Nothing} = nothing`
    A string indicating the desired timeframe of a period for the analysis for staggered adoption scenarios.
    Options are: `"year"`, `"month"`, `"week"`, `"day"`.
- `freq_multiplier::Number = 1`
    An integer by which the 'freq' argument should be multiplied in a staggered adoption scenario, e.g. if a two-year
    period is desired, set `freq = "year"` and `freq_multiplier = 2`.
- `autoadjust::Bool = false`
    A boolean option for the length of a period to be considered in staggered adoption scenario to be automatically
    determined by DiDInt.
- `nperm::Number = 1000`
    The number of unique permutations to be considered when performing the randomization inference.
- `verbose::Bool = true`
    A boolean option for displaying progress of the randomization procedure. 
- `seed::Number = rand(1:1000000)`
    An integer to set the random seed for the randomization inference procedure.

# Returns
A DataFrame of results including the estimate of the ATT as well as standard errors and p-values.

"""
function didint(outcome::Union{AbstractString, Symbol},
                state::Union{AbstractString, Symbol},
                time::Union{AbstractString, Symbol},
                data::DataFrame;
                gvar::Union{AbstractString, Symbol, Nothing} = nothing,
                treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
                treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
                date_format::Union{AbstractString, Nothing} = nothing,
                covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
                ccc::AbstractString = "int",
                agg::AbstractString = "cohort",
                weighting::AbstractString = "both",
                ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
                freq::Union{AbstractString, Nothing} = nothing,
                freq_multiplier::Number = 1,
                start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
                nperm::Number = 1001,
                verbose::Bool = true,
                seed::Number = rand(1:1000000),
                use_pre_controls::Bool = false,
                notyet::Union{Nothing, Bool} = nothing,
                hc::Union{AbstractString, Number} = "hc3")

    # Check hc args
    hc = hc_checks(hc)

    # Let notyet override use_pre_controls
    use_pre_controls = isnothing(notyet) ? use_pre_controls : notyet

    # Turn outcome, state, time, and gvar to strings if they were inputted as symbols
    outcome, state, time, gvar = symbol_to_string(outcome, state, time, gvar)

    # Determine if the gvar method or the treatment_times & treated_states method is being used
    treated_states, treatment_times = gvar_or_treatment_times(gvar, treatment_times, treated_states, data, state;
                                                              event = false, check = "estimate_didint")

    # Do checks for start_date and end_date
    start_date = start_end_date_checks(start_date, "start_date", date_format)
    end_date = start_end_date_checks(end_date, "end_date", date_format)           

    # Check that seed is set correctly
    seed = Int(round(seed))
    if seed <= 0 
        error("'seed' must be a positive integer.")
    end 

    # Check that agg args are passed correctly
    agg = lowercase(agg)
    agg_options = ["cohort", "state", "simple", "none", "sgt", "time"]
    if !(agg in agg_options)
        error("'agg' must be one of $(agg_options)")
    end 

    # Check that weighting arg is passed correctly
    weighting = lowercase(weighting)
    weighting_options = ["both", "none", "att", "diff"]
    if !(weighting in weighting_options) 
        error("'weighting' must be one of $(weighting_options)")
    end

    # Round nperm
    nperm = Int(round(nperm))
    if nperm < 1
        error("'nperm' must be a positive integer > 0.")
    end 

    # Create data copy
    data_copy = DataFrame(data; copycols = false)    

    # Do validation checks to see if specified columns exist, and that they are able to be processed correctly
    data_copy, treatment_times = validate_data(data_copy, outcome, state, time, covariates, treatment_times, date_format;
                                               warn_missing_covariates = false)
    treatment_times = !(treatment_times isa AbstractVector) ? [treatment_times] : treatment_times

    # Classify the requested estimation as either a staggered adoption or common adoption scenario
    if length(unique(treatment_times)) == 1
        common_adoption = true
        staggered_adoption = false
    elseif length(unique(treatment_times)) > 1 
        common_adoption = false
        staggered_adoption = true
    else
        error("'treatment_times' must have at least one entry.")
    end
    
    # Make sure treated_states and control_states are vectors, make sure control_states exist
    treated_states = !(treated_states isa AbstractVector) ? [treated_states] : treated_states 
    unique_states = unique(data_copy[!, state])
    control_states = setdiff(unique_states, treated_states)
    if !(control_states isa AbstractVector)
        control_states = [control_states]
    end 
    missing_control_states = isempty(control_states)
    if missing_control_states && !use_pre_controls
        error("No control states were found.")
    elseif missing_control_states && use_pre_controls
        # If there are at least two unique treatment times, then should still be able to compute
        # some ATTs if using the notyet treated cells as controls 
        if !staggered_adoption
            error("No valid control states were found.")
        end 
    end

    # Do intial check before filtering that treatment_times length is equal treated_states length
    if staggered_adoption && length(treatment_times) != length(treated_states)
        error("'treatment_times' should be the same length as the 'treated_states'.")
    end 

    # Check missing values
    data_copy = check_missing_vals(data_copy, state, time, covariates)

    # Ensure the state column is a string or number and that the nonmissingtype(treated_states) == nonmissingtype(state column)
    validate_state_types(treated_states, data_copy, state)
    data_copy.state_71X9yTx = data_copy[!, state]

    # Do some checks after (potentially) dropping rows to make sure that the analysis
    # is still valid
    if (length(unique(data_copy.state_71X9yTx)) == 1)
        error("Not enough non-missing observations.\nFound", length(unique(data_copy.state_71X9yTx)), "states after dropping missing observations.")
    end 

    # Validate and convert dates to Date objects, filter data by date range
    data_copy, treatment_times, start_date, end_date, all_times, freq = validate_and_convert_dates(data_copy, time, treatment_times, date_format,
                                                                                                   freq, start_date, end_date) 

    # In the case of staggered adoption, do date matching procedure should be done
    if staggered_adoption
        data_copy, treatment_times, treated_states, match_to_these_dates, time_to_index, period = perform_date_matching(data_copy, all_times, freq, freq_multiplier,
                                                                                                                        start_date, end_date, treatment_times,
                                                                                                                        treated_states)
    end

    # Check that treated states actually exist
    treated_states, treatment_times = validate_treated_states(treated_states, treatment_times, data_copy)

    # Do some checks for treatment_times vector length
    if common_adoption && length(treatment_times) != length(treated_states)
        if length(treatment_times) != 1
            error("'treatment_times' should either be the same length as the 'treated_states' vector or of length 1.")
        else
            treatment_times = fill(treatment_times[1], length(treated_states))
        end 
    elseif staggered_adoption && length(treatment_times) != length(treated_states)
        error("'treatment_times' should be the same length as the 'treated_states'.")
    end 

    # Also need to make sure that start_times < treatment_times < end_times is true for each state
    for i in eachindex(treated_states)
        s = treated_states[i]
        treat_time = treatment_times[i]
        state_dates = data_copy[data_copy.state_71X9yTx .== s, :time_71X9yTx]
        earliest = minimum(state_dates)
        latest = maximum(state_dates)
        if !(earliest < treat_time)
            earliest = string(earliest)
            treat_time = string(treat_time)
            error("For state $s, the earliest date for which there is non-missing value ($earliest) is not strictly less than the treatment time ($treat_time).")
        end
        if !(treat_time <= latest)
            treat_time = string(treat_time)
            latest = string(latest)
            error("For state $s, the treatment time ($treat_time) is greater than the last date ($latest) for which there is a non-missing observation.")
        end
    end

    # Ensure state column is a string (as are treated_states)
    data_copy, treated_states = validate_string_treated_states(data_copy; treated_states = treated_states)

    # Once any date matching procedures are done, convert time_71X9yTx back to a string for processing in `categorical()`
    # Also keep a column vector copy as a date
    data_copy.time_dmG5fpM = data_copy.time_71X9yTx
    if common_adoption
        data_copy.time_71X9yTx = ifelse.(data_copy.time_71X9yTx .>= treatment_times[1], "post", "pre")
        data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)
        data_copy = check_states_for_common_adoption(data_copy, treated_states)
    elseif staggered_adoption
        data_copy.time_71X9yTx = string.(data_copy.time_71X9yTx)
    else
        error("A non-common or non-staggered adoption scenario was discovered!?")
    end 

    # Create dummies for each time and state interaction 
    data_copy.state_time = categorical(data_copy.state_71X9yTx .* "0IQR7q6Wei7Ejp4e" .* data_copy.time_71X9yTx)

    # Convert factor covariates into multiple numeric dummy variable columns
    data_copy, covariates_to_include = process_covariates(covariates, data_copy, ref)

    # Force outcome to float64 to speed up regression (runs faster if <:Number rather than <:Union{Number, Missing})
    data_copy.outcome_71X9yTx = convert(Vector{Float64}, data_copy.outcome_71X9yTx)

    # Construct formula, depending on DID-INT variation
    formula = construct_formula(ccc, covariates_to_include)

    # Run the fixed effects model and get back the dataframe of means (or means residualized by covariates) for each period at each state
    lambda_df, cornercase, cornercase_se = run_fixed_effects_model(data_copy, formula, ccc, covariates,
                                                                   common_adoption = common_adoption, staggered_adoption = staggered_adoption) 

    # Compute diff for each treated state
    unique_states = unique(lambda_df.state)
    diff_df = DataFrame(state = String[], treated_time = Date[],
                        t = Date[], r1 = Date[], diff = Float64[],
                        treat = Int[], n = Int[], n_t = Int[])

    # Define a function that initializes the results dataframe columns
    init_column() = Vector{Union{Missing, Float64}}(missing, nrows)

    # Compute ATT
    if staggered_adoption
        # Compute diffs for each treated_state
        for i in eachindex(treated_states)
            idx = time_to_index[treatment_times[i]]
            one_period_prior_treatment = match_to_these_dates[idx - 1]
            temp = lambda_df[(lambda_df.state .== treated_states[i]) .& (lambda_df.time .>= one_period_prior_treatment), :]
            lambda_r1 = temp[temp.time .== one_period_prior_treatment, "lambda"]
            lambda_r1 = isempty(lambda_r1) ? missing : lambda_r1[1]
            temp = temp[temp.time .> one_period_prior_treatment, :]
            sort!(temp, :time)
            years = temp.time
            diffs = Vector{Union{Float64, Missing}}(undef, nrow(temp))
            n = Vector{Int}(undef, nrow(temp))
            n_t = Vector{Int}(undef, nrow(temp))
            for j in eachindex(diffs)
                lambda_j = temp[j, "lambda"]
                lambda_j = isempty(lambda_j) ? missing : lambda_j[1]
                diffs[j] = lambda_j - lambda_r1
                post = temp[j, "time"]
                pre = one_period_prior_treatment
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([post, pre])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
                n_t[j] = sum((data_copy.time_dmG5fpM .== post) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            end 
            temp_df = DataFrame(state = treated_states[i], treated_time = treatment_times[i],
                                t = years, r1 = one_period_prior_treatment, diff = diffs, treat = 1,
                                n = n, n_t = n_t)
            diff_df = vcat(diff_df, temp_df) 
        end
        
        # Compute diffs for control states
        unique_diffs = unique(select(diff_df, :t, :r1, :treated_time))
        control_states = setdiff(unique_states, treated_states)     
        for i in eachindex(control_states)
            temp = lambda_df[(lambda_df.state .== control_states[i]), :]
            diffs = Vector{Union{Float64, Missing}}(undef, nrow(unique_diffs))
            n = Vector{Int}(undef, nrow(unique_diffs))
            n_t = Vector{Int}(undef, nrow(unique_diffs))
            for j in 1:nrow(unique_diffs)
                t = unique_diffs[j,"t"]
                r1 = unique_diffs[j,"r1"]
                lambda_t = temp[temp.time .== t, "lambda"]
                lambda_t = isempty(lambda_t) ? missing : lambda_t[1]
                lambda_r1 = temp[temp.time .== r1, "lambda"]
                lambda_r1 = isempty(lambda_r1) ? missing : lambda_r1[1]
                diffs[j] = lambda_t - lambda_r1
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== control_states[i]))
                n_t[j] = sum((data_copy.time_dmG5fpM .== t) .&& (data_copy.state_71X9yTx .== control_states[i]))
            end 
            trtd_time = [match_to_these_dates[findfirst(==(t), match_to_these_dates) + 1] for t in unique_diffs.r1]
            temp_df = DataFrame(state = control_states[i], treated_time = trtd_time,
                                t = unique_diffs.t, r1 = unique_diffs.r1, diff = diffs, treat = 0,
                                n = n, n_t = n_t)
            diff_df = vcat(diff_df, temp_df)
        end

        # Compute diffs for randomization inference
        ri_diff_df = DataFrame(state = String[], treated_time = Date[],
                               t = Date[], r1 = Date[], diff = Float64[],
                               treat = Int[], n = Int[], n_t = Int[])
        for i in eachindex(treated_states)
            ri_diffs = unique(select(diff_df[(diff_df.state .!= treated_states[i]) .& (diff_df.treated_time .!= treatment_times[i]), :], :t, :r1))
            temp = lambda_df[(lambda_df.state .== treated_states[i]), :]
            diffs = Vector{Union{Float64, Missing}}(undef, nrow(ri_diffs))
            n = Vector{Int}(undef, nrow(ri_diffs))
            n_t = Vector{Int}(undef, nrow(ri_diffs))
            for j in 1:nrow(ri_diffs)
                t = unique_diffs[j,"t"]
                r1 = unique_diffs[j,"r1"]
                lambda_t = temp[temp.time .== t, "lambda"]
                lambda_t = isempty(lambda_t) ? missing : lambda_t[1]
                lambda_r1 = temp[temp.time .== r1, "lambda"]
                lambda_r1 = isempty(lambda_r1) ? missing : lambda_r1[1]
                diffs[j] = lambda_t - lambda_r1
                n[j] = sum(in.(data_copy.time_dmG5fpM, Ref([t, r1])) .&& (data_copy.state_71X9yTx .== treated_states[i]))
                n_t[j] = sum((data_copy.time_dmG5fpM .== t) .&& (data_copy.state_71X9yTx .== treated_states[i]))
            end
            trtd_time = [match_to_these_dates[findfirst(==(t), match_to_these_dates) + 1] for t in ri_diffs.r1]
            temp_df = DataFrame(state = treated_states[i], treated_time = trtd_time,
                                t = ri_diffs.t, r1 = ri_diffs.r1, diff = diffs, treat = -1,
                                n = n, n_t = n_t)
            ri_diff_df = vcat(ri_diff_df, temp_df)
        end

        # Add the periods since treatment column to ri_diff_df and diff_df
        ordered_t = match_to_these_dates   
        t_period = indexin(diff_df.t, ordered_t)             
        treated_time_period = indexin(diff_df.treated_time, ordered_t) 
        diff_df.time_since_treatment = t_period .- treated_time_period 
        t_period = indexin(ri_diff_df.t, ordered_t)             
        treated_time_period = indexin(ri_diff_df.treated_time, ordered_t) 
        ri_diff_df.time_since_treatment = t_period .- treated_time_period 

        # If use_pre_controls is true (default) re-structure the ri_diff_df and diff_df
        # to add the pre-treatment periods from treated silos as control contrasts
        if use_pre_controls
            for (t,s) in zip(treatment_times, treated_states)
                ri_diff_df[ri_diff_df.state .== s .&& ri_diff_df.t .< t, "treat"] .= 0
            end
            diff_df = vcat(diff_df, ri_diff_df[ri_diff_df.treat .== 0,:])
            ri_diff_df = ri_diff_df[ri_diff_df.treat .== -1,:]
            if missing_control_states
                # If there are no pure control states, do extra processing to remove cells
                # where there are no treat-control pairs
                diff_df = subset(
                            groupby(diff_df, [:t, :r1]),
                            :treat => x -> any(x .== 1) && any(x .== 0)
                          )
                ri_diff_df = semijoin(ri_diff_df, diff_df, on=[:t, :r1])
            end 
        end

        # In the case where the diff_df has missing values for diffs, drop those rows
        if any(ismissing.(diff_df.diff))
            diff_df = dropmissing(diff_df, :diff)
            diff_df = subset(
                            groupby(diff_df, [:t, :r1]),
                            :treat => x -> any(x .== 1) && any(x .== 0)
                          )
            if nrow(diff_df) == 0
                error("Could not find a single pair of non-missing valued treated & untreated differences for the same (g,t) group.\nUnable to compute any ATTs.")
            end
            ri_diff_df = semijoin(ri_diff_df, diff_df, on=[:t, :r1])
            ri_diff_df = dropmissing(ri_diff_df, :diff)
        end

        # Do check to make sure that each state has each (g,t) group for RI 
        randomize = check_prior_to_ri(diff_df, ri_diff_df)

        # These are used in computing the sub-agg ATTs so, overwrite them if we dropped missing vals
        unique_diffs = unique(select(diff_df, :t, :r1, :treated_time))
        original_treated = sort(unique(diff_df[diff_df.treat .== 1, [:state, :treated_time]]), :treated_time)
        treatment_times = original_treated.treated_time
        treated_states = original_treated.state

        # Show period length in results
        if isnothing(date_format) || date_format == "yyyy"
            date_format = assume_date_format(period)
        end
        period = string(period)
        start_date = parse_date_to_string_didint(start_date, date_format)
        end_date = parse_date_to_string_didint(maximum(data_copy.time_dmG5fpM), date_format)

        # Run final regression to compute ATT based on weighting/aggregation method
        if agg == "cohort"

            # Define nrows and vector to iterate through for sub aggregate ATTs
            unique_treatment_times = sort(unique(diff_df.treated_time))
            nrows = length(unique_treatment_times)

            # Create results df
            results = DataFrame(treatment_time = Vector{Date}(undef, nrows),
                                att_cohort = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_cohort = init_column(),
                                pval_att_cohort = init_column(),
                                jknifese_att_cohort = init_column(),
                                jknifepval_att_cohort = init_column(),
                                ri_pval_att_cohort = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Compute the sub aggregate ATTs
            for i in eachindex(unique_treatment_times)
                trt = unique_treatment_times[i]
                temp = diff_df[diff_df.treated_time .== trt, :]
                X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end            
                results.treatment_time[i] = trt
                result_dict = final_regression_results(X, Y, W = W, hc = hc)
                results.att_cohort[i] = result_dict["beta_hat"]
                results.se_att_cohort[i] = result_dict["beta_hat_se"]
                results.pval_att_cohort[i] = result_dict["pval_att"] 
                results.jknifese_att_cohort[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_cohort[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Compute aggregate ATT and return results
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_cohort)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = !randomize ? results : randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "cohort",
                                                      verbose, seed, data_copy, weighting, use_pre_controls)

        elseif agg == "simple"

            # Define nrows and object to iterate through for sub aggregate ATTs (unique_diffs)
            nrows = nrow(unique_diffs)

            # Create results df
            results = DataFrame(r1 = Vector{Date}(undef, nrows),
                                time = Vector{Date}(undef, nrows),
                                gvar = Vector{Date}(undef, nrows),
                                att_gt = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_gt = init_column(),
                                pval_att_gt = init_column(),
                                jknifese_att_gt = init_column(),
                                jknifepval_att_gt = init_column(),
                                ri_pval_att_gt = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))
            
            # Compute the sub aggregate ATTs
            for i in 1:nrow(unique_diffs)
                t = unique_diffs[i,"t"]
                r1 = unique_diffs[i,"r1"]
                gvar = unique_diffs[i, "treated_time"]
                temp = diff_df[(diff_df.t .== t) .& (diff_df.r1 .== r1), :]
                X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.time[i] = t
                results.r1[i] = r1
                results.gvar[i] = gvar
                result_dict = final_regression_results(X, Y, W = W, hc = hc)
                results.att_gt[i] = result_dict["beta_hat"]
                results.se_att_gt[i] = result_dict["beta_hat_se"]
                results.pval_att_gt[i] = result_dict["pval_att"] 
                results.jknifese_att_gt[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_gt[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Sort results dataframe and compute aggregate ATT, return results
            sort!(results, [:time])
            sort!(results, [:r1])
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_gt)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = !randomize ? results : randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "simple",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)

        elseif agg == "state"

            # Define nrows and identify object to iterate thru for the sub aggregate ATTs
            nrows = length(treated_states)

            # Create results df
            results = DataFrame(state = Vector{String}(undef, nrows),
                                att_s = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_s = init_column(),
                                pval_att_s = init_column(),
                                jknifese_att_s = init_column(),
                                jknifepval_att_s = init_column(),
                                ri_pval_att_s = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Compute sub aggregate ATTs
            for i in eachindex(treated_states)
                state = treated_states[i]
                trt = treatment_times[i]
                temp_treated = diff_df[diff_df.state .== state, :]
                temp_control = diff_df[(diff_df.treat .== 0) .&& (diff_df.treated_time .== trt), :]
                temp = vcat(temp_control, temp_treated)
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.state[i] = state
                result_dict = final_regression_results(X, Y, W = W, hc = hc)
                results.att_s[i] = result_dict["beta_hat"]
                results.se_att_s[i] = result_dict["beta_hat_se"]
                results.pval_att_s[i] = result_dict["pval_att"] 
                results.jknifese_att_s[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_s[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Sort results dataframe and compute aggregate ATT, return results
            results.tuple_state = custom_sort_order.(results.state)
            sort!(results,[order(:tuple_state)])
            select!(results, Not([:tuple_state]))
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_s)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = !randomize ? results : randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "state",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)

        elseif agg == "none"  
            
            # Create results df
            results = DataFrame(agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                nperm = Vector{Union{Missing, Float64}}(missing, 1))
            
            # Compute aggregate ATT, return results
            X = convert(Matrix{Float64}, hcat(ones(nrow(diff_df)), diff_df.treat))
            Y = convert(Vector{Float64}, diff_df.diff)
            if in(weighting, ["diff", "both"])
                W = convert(Vector{Float64}, diff_df.n)
                W ./= sum(W)
            else
                W = fill(nothing, nrow(diff_df))
            end 
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = !randomize ? results : randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "none",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)

        elseif agg == "sgt"
            
            # Compute nrows and define object to iterate thru for the sub aggregate ATTs
            unique_sgt = unique(select(diff_df[diff_df.treat .== 1,:], :state, :t, :treated_time))
            nrows = nrow(unique_sgt)

            # Create results df
            results = DataFrame(state = Vector{String}(undef, nrows),
                                gvar = Vector{Date}(undef, nrows),
                                t = Vector{Date}(undef, nrows),
                                att_sgt = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_sgt = init_column(),
                                pval_att_sgt = init_column(),
                                jknifese_att_sgt = init_column(),
                                jknifepval_att_sgt = init_column(),
                                ri_pval_att_sgt = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))
            
            # Compute the sub aggregate ATTs 
            for i in 1:nrow(unique_sgt)
                state = unique_sgt[i, "state"]
                gvar = unique_sgt[i, "treated_time"]
                t = unique_sgt[i, "t"]
                temp_treated = diff_df[(diff_df.state .== state) .&& (diff_df.treated_time .== gvar) .&& (diff_df.t .== t),:]
                temp_control = diff_df[(diff_df.treat .== 0) .&& (diff_df.treated_time .== gvar) .&& (diff_df.t .== t),:]
                temp = vcat(temp_treated, temp_control)
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.state[i] = state
                results.gvar[i] = gvar
                results.t[i] = t
                result_dict = final_regression_results(X, Y, W = W, hc = hc)
                results.att_sgt[i] = result_dict["beta_hat"]
                results.se_att_sgt[i] = result_dict["beta_hat_se"]
                results.pval_att_sgt[i] = result_dict["pval_att"] 
                results.jknifese_att_sgt[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_sgt[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end 

            # Compute aggregate ATT, return results
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_sgt)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = !randomize ? results : randomization_inference_v2(vcat(diff_df, ri_diff_df), nperm, results, "sgt",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)

        elseif agg == "time"

            # Compute nrows, define times to iterate thru to get sub aggregate ATTs
            nrows = length(unique(diff_df.time_since_treatment))
            times = sort(unique(diff_df.time_since_treatment))

            # Create results df
            results = DataFrame(periods_post_treat = init_column(),
                                att_t = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_t = init_column(),
                                pval_att_t = init_column(),
                                jknifese_att_t = init_column(),
                                jknifepval_att_t = init_column(),
                                ri_pval_att_t = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Special procedure for adding time or gvar dummies (yields equivalent results)
            diff = vcat(diff_df, ri_diff_df)
            treated_category = categorical(string.(diff.treated_time); ordered = true)
            levels_tt   = levels(treated_category)         
            dummy_cols  = Symbol[]                         
            for i in eachindex(levels_tt) 
                tt = levels_tt[i]                    
                colname = Symbol("tt_", i)
                diff[!, colname] = Int.(treated_category .== tt)
                push!(dummy_cols, colname)
            end

            # Compute sub aggregate ATTs
            for i in eachindex(times)
                t = times[i]
                temp = diff[(diff.time_since_treatment .== t) .&& (diff.treat .!= -1), :]
                X = design_matrix_time_agg(temp, dummy_cols, Symbol("treat"))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.periods_post_treat[i] = t
                result_dict = final_regression_results(X, Y, W = W, hc = hc)
                results.att_t[i] = result_dict["beta_hat"]
                results.se_att_t[i] = result_dict["beta_hat_se"]
                results.pval_att_t[i] = result_dict["pval_att"] 
                results.jknifese_att_t[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_t[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Compute aggregate ATT, return results
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_t)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = !randomize ? results : randomization_inference_v2(diff, nperm, results, "time",
                                                 verbose, seed, data_copy, weighting, use_pre_controls,
                                                 dummy_cols = dummy_cols)
        end

        # Return staggered adoption results
        results.period .= period
        results.start_date .= start_date
        results.end_date .= end_date
        return results

    elseif common_adoption

        # Create diff matrix
        states = unique(lambda_df.state)
        diff_df = DataFrame(state = states,
                            diff = Vector{Float64}(undef, length(states)),
                            treat = Vector{Float64}(undef, length(states)),
                            treated_time = fill(unique(treatment_times)[1], length(states)),
                            n = Vector{Float64}(undef, length(states)),
                            n_t = Vector{Float64}(undef, length(states)))
        for i in eachindex(states)
            state = states[i]
            if state in treated_states
                trt = 1
            else 
                trt = 0
            end 
            diff = lambda_df[(lambda_df.state .== state) .&& (lambda_df.time .== "post"), "lambda"][1] - lambda_df[(lambda_df.state .== state) .&& (lambda_df.time .== "pre"), "lambda"][1]
            diff_df.diff[i] = diff
            diff_df.treat[i] = trt
            diff_df.n[i] = sum(data_copy.state_71X9yTx .== state)
            diff_df.n_t[i] = sum((data_copy.state_71X9yTx .== state) .&& (data_copy.time_71X9yTx .== "post"))
        end 

        # Compute results
        if in(agg, ["cohort", "gt", "none"]) || length(treated_states) == 1
            results = DataFrame(agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                se_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifese_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                jknifepval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                                nperm = Vector{Union{Missing, Float64}}(missing, 1))
            X = convert(Matrix{Float64}, hcat(ones(nrow(diff_df)), diff_df.treat))
            Y = convert(Vector{Float64}, diff_df.diff)
            if in(weighting, ["diff", "both"])
                W = convert(Vector{Float64}, diff_df.n)
                W ./= sum(W)
            else
                W = fill(nothing, nrow(diff_df))
            end 
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            if cornercase
                results.se_agg_att[1] = cornercase_se
            else
                results.se_agg_att[1] = result_dict["beta_hat_se"]
            end 
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(diff_df, nperm, results, agg,
                                                 verbose, seed, data_copy, weighting,
                                                 use_pre_controls)

            return results
        elseif in(agg, ["state", "sgt"])

            # Define nrows and identify object to iterate thru for the sub aggregate ATTs
            nrows = length(treated_states)

            # Create results df
            results = DataFrame(state = Vector{String}(undef, nrows),
                                att_s = init_column(),
                                agg_att = init_column(),
                                se_agg_att = init_column(),
                                pval_agg_att = init_column(),
                                jknifese_agg_att = init_column(),
                                jknifepval_agg_att = init_column(),
                                ri_pval_agg_att = init_column(),
                                nperm = init_column(),
                                se_att_s = init_column(),
                                pval_att_s = init_column(),
                                jknifese_att_s = init_column(),
                                jknifepval_att_s = init_column(),
                                ri_pval_att_s = init_column(),
                                weights = Vector{Union{Nothing, Float64}}(nothing, nrows))

            # Compute sub aggregate ATTs
            for i in eachindex(treated_states)
                state = treated_states[i]
                trt = treatment_times[i]
                temp_treated = diff_df[diff_df.state .== state, :]
                temp_control = diff_df[(diff_df.treat .== 0) .&& (diff_df.treated_time .== trt), :]
                temp = vcat(temp_control, temp_treated)
                X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
                Y = convert(Vector{Float64}, temp.diff)
                if in(weighting, ["diff", "both"])
                    W = convert(Vector{Float64}, temp.n)
                    W ./= sum(W)
                else
                    W = fill(nothing, nrow(temp))
                end 
                results.state[i] = state
                result_dict = final_regression_results(X, Y, W = W, hc = hc)
                results.att_s[i] = result_dict["beta_hat"]
                results.se_att_s[i] = result_dict["beta_hat_se"]
                results.pval_att_s[i] = result_dict["pval_att"] 
                results.jknifese_att_s[i] = result_dict["beta_hat_se_jknife"]
                results.jknifepval_att_s[i] = result_dict["pval_att_jknife"]
                if weighting in ["att", "both"]
                    results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
                end
            end

            # Sort results dataframe and compute aggregate ATT, return results
            results.tuple_state = custom_sort_order.(results.state)
            sort!(results,[order(:tuple_state)])
            select!(results, Not([:tuple_state]))
            results = scale_weights_final(results, weighting)
            X = ones(nrow(results), 1)
            Y = convert(Vector{Float64}, results.att_s)
            W = results.weights
            result_dict = final_regression_results(X, Y, W = W, hc = hc)
            results.agg_att[1] = result_dict["beta_hat"]
            results.se_agg_att[1] = result_dict["beta_hat_se"]
            results.pval_agg_att[1] = result_dict["pval_att"]
            results.jknifese_agg_att[1] = result_dict["beta_hat_se_jknife"]
            results.jknifepval_agg_att[1] = result_dict["pval_att_jknife"]
            results = randomization_inference_v2(diff_df, nperm, results, "state",
                                                 verbose, seed, data_copy, weighting, use_pre_controls)
            return results

        end 
    end 

end

### The following functions are parent functions within the didint() function:
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

function check_states_for_common_adoption(data_copy, treated_states)

    states = unique(data_copy.state_71X9yTx)
    states_to_drop = []
    for state in states
        periods = unique(data_copy[data_copy.state_71X9yTx .== state, :].time_71X9yTx)
        if !("pre" in periods && "post" in periods)
            @warn "State $state does not have at least one post-treatment and one pre-treatment observation.\nDropping $state."
            push!(states_to_drop, state)
        end
    end

    data_copy = isempty(states_to_drop) ? data_copy : filter(row -> row.state_71X9yTx ∉ states_to_drop, data_copy)

    remaining_states = unique(data_copy.state_71X9yTx)
    n_states = length(remaining_states)
    if n_states < 2
        error("Only found $n_states with observations both post-treatment and pre-treatment!")
    end
    have_treated = any(state -> state in remaining_states, treated_states)
    have_control = any(state -> state ∉ treated_states, remaining_states)

    if have_treated && have_control
        return data_copy
    else
        error("After dropping states that did not have at least one post-treatment and one pre-treatment observation,\nthere was not a single pair of a treated and untreated state required in order to compute ATTs.")
    end

end

function design_matrix_time_agg(temp::DataFrame, dummy_cols::Vector{Symbol}, treat_col::Symbol)

    # This function adds dummy variables corresponding to the different gvars
    # e.g. in the 4 periods since treatment group you could have multiple different 
    # treatment times, so this function accounts for those differences

    # This function is called both from within didint() and during the randomization inference
    
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

function final_regression_results(X::Matrix, Y::Vector; W::Vector = [nothing], hc::AbstractString = "hc3")

    n = length(Y)

    # Check for and remove missing values
    # Build mask for each component
    x_valid = [!any(ismissing, X[i, :]) for i in 1:n]
    y_valid = .!ismissing.(Y)

    if eltype(W) <: Nothing
        valid_mask = x_valid .& y_valid
    else
        w_valid = .!ismissing.(W)
        valid_mask = x_valid .& y_valid .& w_valid
    end

    # Filter out missing observations
    X = X[valid_mask, :]
    Y = Y[valid_mask]
    if eltype(W) <: Number
        W = W[valid_mask]
    end
    
    # Check if we have enough observations left
    ncolx = size(X, 2)
    n = length(Y)
    if n <= ncolx
        @warn "Insufficient observations after dropping missing values (n = $n, k = $ncolx)"
        return Dict("beta_hat" => missing, "beta_hat_se" => missing, 
                   "pval_att" => missing, "beta_hat_se_jknife" => missing, 
                   "pval_att_jknife" => missing)
    end

    beta_hat = nothing
    beta_hat_cov = nothing

    # Run OLS (normally if weights aren't provided, and scale (X,Y) -> (Xw,Yw) otherwise)
    if eltype(W) <: Nothing
        beta_hat = safe_solve(X, Y)
        if ismissing(beta_hat)
            return Dict("beta_hat" => missing, "beta_hat_se" => missing, "pval_att" => missing,
                       "beta_hat_se_jknife" => missing, "pval_att_jknife" => missing)
        end
        resid = Y - X * beta_hat
        beta_hat_cov = compute_hc_covariance(X, resid, hc)
        beta_hat_se_jknife = compute_jknife_se(X, Y, beta_hat[ncolx]) 
    elseif eltype(W) <: Number
        sw = sqrt.(W)           
        Xw = X .* sw            
        Yw = Y .* sw
        beta_hat = safe_solve(Xw, Yw)
        if ismissing(beta_hat)
            return Dict("beta_hat" => missing, "beta_hat_se" => missing, "pval_att" => missing,
                       "beta_hat_se_jknife" => missing, "pval_att_jknife" => missing)
        end
        resid_w = Yw - Xw * beta_hat
        beta_hat_cov = compute_hc_covariance(Xw, resid_w, hc)
        beta_hat_se_jknife = compute_jknife_se(Xw, Yw, beta_hat[ncolx]) 
    end 

    beta_hat_var = diag(beta_hat_cov)
    beta_hat_se = sqrt(beta_hat_var[ncolx]) 
    dof = n - ncolx
    pval_att = dof > 0 ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se))) : missing 
    pval_att_jknife = dof > 0 && !ismissing(beta_hat_se_jknife) ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se_jknife))) : missing
    result_dict = Dict("beta_hat" => beta_hat[ncolx], "beta_hat_se" => beta_hat_se, "pval_att" => pval_att,
                       "beta_hat_se_jknife" => beta_hat_se_jknife, "pval_att_jknife" => pval_att_jknife)
    return result_dict
end

function custom_sort_order(s)
    parsed = tryparse(Int, s)
    if isnothing(parsed)
        (0, lowercase(s))
    else
        (1, parsed)
    end 
end

function check_prior_to_ri(diff_df, ri_diff_df)
    
    randomize = true

    # Get all required (t, r1) combinations and check for each state, each pair exists in the combined data
    combined_df = vcat(diff_df, ri_diff_df)
    all_states = unique(combined_df.state)
    required_pairs = unique(select(combined_df, :t, :r1))
    for state in all_states
        state_data = combined_df[combined_df.state .== state, :]
        state_pairs = unique(select(state_data, :t, :r1))
        missing_pairs = antijoin(required_pairs, state_pairs, on=[:t, :r1])
        if nrow(missing_pairs) > 0
            randomize = false
            break
        end
    end
    if !randomize
        @warn "Randomization inference cannot be performed.\nDid not find all the non-missing values for differences required to implement randomization inference procedure."
    end

    return randomize

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

## compute_jknife_se() is used within final_regression_results()
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

## The following functions are used within randomization_inference_v2():
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



