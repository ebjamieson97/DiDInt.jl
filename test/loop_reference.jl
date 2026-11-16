# Slow reference for the batched engine (update spec §12.2 items 7, 11, 12): the per-draw loops
# that the jackknife and randomization inference used before the engine replaced them. They
# compute the estimates for one assignment at a time, in closed form, from a lookup of long
# differences. Kept only for the tests.

using DiDInt: CovariateFit, DEFAULT_POLICY, group_columns, treated_unit_weights, slope_change_term, policy_pool_weights
using Dates
using LinearAlgebra
using Random
using Statistics

function att_st_point_estimates(cells::DataFrame, lookup::Dict, labels::Dict, states::AbstractVector,
                                weighting::AbstractString, use_pre_controls::Bool;
                                integration::Union{CovariateFit, Nothing} = nothing,
                                policies::Union{AbstractDict, Nothing} = nothing)

    # ATT(s,t)'s for the assignment `labels` (state => treatment period; states not in it are
    # never treated), in closed form: ATT(s,t) = diff_s - sum_c v_c diff_c / sum_c v_c over the
    # controls of the cell (update spec §2, §3.3-§3.4), which equals the ATT(s,t) regression
    # coefficient. Controls are never-treated states and, with use_pre_controls, states not yet
    # treated at t; v = n under weighting = "diff"/"both", else 1.
    # With `integration` the long differences are the integrated ones (DID-INT eqs. 46-47). The
    # control mean for a treated state s is then a + b' xbar_{s,t}, with a the weighted mean of
    # the controls' intercept differences and b that of their slope changes (update spec §7.4).
    # `policies` (treated state => policy) only labels the ATT(s,t)'s: controls do not depend
    # on policy (§6.2).
    out = DataFrame(state = String[], gvar = Date[], t = Date[], r1 = Date[], periods_post_treat = Int[],
                    policy = String[], att = Float64[], n = Float64[], n_t = Float64[])
    n_covariates = isnothing(integration) ? 0 : length(integration.covariates)
    for cell in eachrow(cells)
        g, t, r1 = cell.gvar, cell.t, cell.r1
        weighted_sum = 0.0
        weighted_slope_change = zeros(n_covariates)
        weight_total = 0.0
        for c in states
            label = get(labels, c, nothing)
            (isnothing(label) || (use_pre_controls && t < label)) || continue
            d = get(lookup, (c, t, r1), nothing)
            isnothing(d) && continue
            v = weighting in ["diff", "both"] ? d.n : 1.0
            weighted_sum += v * d.diff
            if !isnothing(integration)
                weighted_slope_change .+= v .* (integration.slopes[(c, t)] .- integration.slopes[(c, r1)])
            end
            weight_total += v
        end
        weight_total == 0 && continue
        control_mean = weighted_sum / weight_total
        control_slope_change = weighted_slope_change ./ weight_total
        for s in states
            get(labels, s, nothing) == g || continue
            d = get(lookup, (s, t, r1), nothing)
            isnothing(d) && continue
            att = d.diff - control_mean
            if !isnothing(integration)
                att += slope_change_term(integration, s, t, r1, s) - dot(control_slope_change, integration.xbar[(s, t)])
            end
            push!(out, (s, g, t, r1, cell.periods_post_treat, isnothing(policies) ? DEFAULT_POLICY : policies[s], att, d.n, d.n_t))
        end
    end
    return out
end

function point_sub_aggregates(est::DataFrame, agg::AbstractString, weighting::AbstractString,
                              treated_weight::AbstractString)

    # u-weighted mean ATT and Level 2 weight (sum of n_t) of every group with at least one
    # ATT(s,t), keyed by the group's key tuple (same weights as sub_aggregate_results)
    group_cols = group_columns(agg)
    u = treated_unit_weights(est, weighting, treated_weight)
    sums = Dict{Tuple, Vector{Float64}}()
    for i in 1:nrow(est)
        key = Tuple(est[i, col] for col in group_cols)
        acc = get!(sums, key, zeros(3))
        acc[1] += u[i] * est.att[i]
        acc[2] += u[i]
        acc[3] += est.n_t[i]
    end
    return Dict(key => (att = acc[1] / acc[2], weight = acc[3]) for (key, acc) in sums)
end

function point_aggregate(sub::Dict, weighting::AbstractString)

    # Aggregate ATT: Level 2 weighted mean of the defined sub-aggregates (update spec §4.2);
    # NaN if there are none
    isempty(sub) && return NaN
    entries = [sub[key] for key in sort(collect(keys(sub)))]
    atts = [e.att for e in entries]
    if weighting in ["att", "both"]
        weights = [e.weight for e in entries]
        return sum(weights .* atts) / sum(weights)
    else
        return mean(atts)
    end
end

function point_policy_aggregates(sub::Dict, est::DataFrame, policies::AbstractVector, weighting::AbstractString,
                                 policy_weighting::AbstractString)

    # Each policy's ATT_p and the pooled ATT for the policy aggregations, from the point
    # sub-aggregates of one deletion or draw (the policy is the first element of every key).
    # ATT_p is the Level 2 weighted mean of the policy's sub-aggregates; the pooled ATT weights
    # the defined ATT_p's by omega_p, recomputed from this deletion's / draw's ATT(s,t)'s (update
    # spec §6.3-§6.4). A policy without any ATT(s,t) gives NaN and is left out of the pool.
    att_p = [point_aggregate(filter(entry -> entry.first[1] == p, sub), weighting) for p in policies]
    defined = isfinite.(att_p)
    pooled = any(defined) ? sum(policy_pool_weights(est, policies[defined], policy_weighting) .* att_p[defined]) : NaN
    return att_p, pooled
end

function state_level_ri_draw(est::DataFrame, result_keys, observed_labels::Dict, agg::AbstractString,
                             weighting::AbstractString, treated_weight::AbstractString;
                             picks = Dict{String, Union{String, Nothing}}())

    # State-level RI statistics for one draw (agg = "state" or "sgt"; update spec D9, §7.5). For
    # each observed treated state s, take its treatment time g_s and pick at random one of the
    # states the draw assigns g_s; that state stands in for s in every period. Its estimate is
    # the u-weighted mean of its ATT(s,t)'s ("state") or its ATT(s,t) in each of s's periods
    # ("sgt"). The aggregate is the Level 2 weighted mean of these statistics. `picks` (observed
    # treated state => stand-in) fixes the stand-ins instead of drawing them.
    stats = fill(NaN, length(result_keys))
    level2 = zeros(length(result_keys))
    u = treated_unit_weights(est, weighting, treated_weight)
    for (i, key) in enumerate(result_keys)
        s = key[1]
        g = observed_labels[s]
        if !haskey(picks, s)
            candidates = unique(est.state[est.gvar .== g])
            picks[s] = isempty(candidates) ? nothing : shuffle(candidates)[1]
        end
        pick = picks[s]
        isnothing(pick) && continue
        rows = agg == "state" ? (est.state .== pick) .& (est.gvar .== g) : (est.state .== pick) .& (est.gvar .== g) .& (est.t .== key[3])
        any(rows) || continue
        stats[i] = sum(u[rows] .* est.att[rows]) / sum(u[rows])
        level2[i] = sum(est.n_t[rows])
    end
    defined = isfinite.(stats)
    if !any(defined)
        return stats, NaN
    elseif weighting in ["att", "both"]
        return stats, sum(level2[defined] .* stats[defined]) / sum(level2[defined])
    else
        return stats, mean(stats[defined])
    end
end
