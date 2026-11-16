### Batched engine for randomization inference and the fast jackknife (update spec §7, §9.2)
#
# Every estimate didint() reports is a weighted mean of ATT(s,t)'s, and for any assignment of
# treatment times to states
#   ATT(s,t) = Delta_{s,t} - sum_c v_c Delta_{c->s,t} / sum_c v_c        (DID-INT eq. 48)
# is built from first-stage quantities that do not depend on the assignment (the lambdas, the
# slopes and the covariate means; §7.1). The engine stores those once, as (state x cell)
# matrices, and evaluates all draws (or all jackknife deletions) at once: one row of the
# assignment matrix `L` per draw, a few matrix products per (g,t) cell.

# Label of a never-treated state in an assignment matrix. Treated states carry the code of
# their treatment period, which is positive and ordered like the dates.
const NEVER_TREATED = 0
period_code(date::Date) = Dates.value(date)

# First-stage quantities per state (rows, in the fixed order `states`) and cell (columns, the
# rows of `cells`), named as in update spec §7.3. Entries without a long difference are zero
# and flagged in `valid`.
struct EngineData
    states::Vector{String}
    cells::DataFrame                  # gvar, t, r1, periods_post_treat of each (g,t) cell
    cohort::Vector{Int}               # period code of each cell's cohort g
    period::Vector{Int}               # period code of each cell's period t
    Dlambda::Matrix{Float64}          # intercept differences lambda_{s,t} - lambda_{s,r1}
    Dint::Matrix{Float64}             # a state's own integrated long difference (eq. 46)
    DGamma::Vector{Matrix{Float64}}   # per cell, S x p: slope differences gamma_{s,t} - gamma_{s,r1}
    Xb::Vector{Matrix{Float64}}       # per cell, S x p: covariate means xbar_{s,t}
    V::Matrix{Float64}                # weight of a state's long difference as a control (§3.3)
    U::Matrix{Float64}                # Level 1 weight u of a state's ATT(s,t) (§4.1)
    NT::Matrix{Float64}               # n_{s,t}, the Level 2 weights (§4.2)
    valid::BitMatrix                  # the state has a long difference in the cell
end

function build_engine(long_diffs::AbstractDataFrame, cells::DataFrame, states::AbstractVector,
                      weighting::AbstractString, treated_weight::AbstractString;
                      integration::Union{CovariateFit, Nothing} = nothing)

    # Precompute the engine's matrices (update spec §7.3) from the long differences of every
    # state in every cell. `integration` is only given under ccc = "int" with covariates: for
    # every other ccc the covariate term cancels in ATT(s,t) (D5), so DGamma and Xb stay empty
    # and the engine skips them, which changes nothing in the results (§7.4).
    lookup = long_difference_lookup(long_diffs)
    S, K = length(states), nrow(cells)
    p = isnothing(integration) ? 0 : length(integration.covariates)
    Dlambda = zeros(S, K); Dint = zeros(S, K); N = zeros(S, K); NT = zeros(S, K)
    valid = falses(S, K)
    DGamma = [zeros(S, p) for _ in 1:(p > 0 ? K : 0)]
    Xb = [zeros(S, p) for _ in 1:(p > 0 ? K : 0)]
    for k in 1:K, (i, s) in enumerate(states)
        t, r1 = cells.t[k], cells.r1[k]
        d = get(lookup, (s, t, r1), nothing)
        isnothing(d) && continue
        valid[i, k] = true
        Dlambda[i, k] = d.diff
        Dint[i, k] = d.diff
        N[i, k] = d.n
        NT[i, k] = d.n_t
        if p > 0
            DGamma[k][i, :] = integration.slopes[(s, t)] .- integration.slopes[(s, r1)]
            Xb[k][i, :] = integration.xbar[(s, t)]
            Dint[i, k] += dot(DGamma[k][i, :], Xb[k][i, :])
        end
    end

    # Control weights v = n under weighting = "diff"/"both", else 1 (§3.3), and the Level 1
    # weights u of §4.1
    V = weighting in ["diff", "both"] ? copy(N) : ones(S, K)
    U = reshape(treated_unit_weights(DataFrame(n = vec(N), n_t = vec(NT)), weighting, treated_weight), S, K)
    return EngineData(String.(states), cells, period_code.(cells.gvar), period_code.(cells.t),
                      Dlambda, Dint, DGamma, Xb, V, U, NT, valid)
end

function assignment_row(states::AbstractVector, treatment_time::AbstractDict, policy_label::AbstractDict,
                        policy_code::AbstractDict)

    # One assignment in the engine's coding: for each state its treatment period (NEVER_TREATED
    # if it has none) and the code of its policy (0 if it has none)
    L = [haskey(treatment_time, s) ? period_code(treatment_time[s]) : NEVER_TREATED for s in states]
    P = [haskey(treatment_time, s) ? policy_code[policy_label[s]] : 0 for s in states]
    return L, P
end

function policy_coding(results::DataFrame, agg::AbstractString, policy_labels::AbstractDict)

    # Integer codes of the policies. The policies that have estimates come first, in the order
    # of the results, so that code j is the j-th reported policy; a policy without any estimate
    # still gets its own code, so that the draws keep it apart from the others. A non-policy
    # aggregation has the single default policy (every treated state carries it; §6.1).
    names = is_policy_agg(agg) ? unique(results.policy) : String[DEFAULT_POLICY]
    others = sort(setdiff(unique(values(policy_labels)), names))
    return names, Dict(name => code for (code, name) in enumerate(vcat(names, others)))
end

### Draws (update spec §7.2, Appendix A)

# Number of distinct arrangements of a multiset with these label counts (BigInt-safe)
n_arrangements(counts) = factorial(big(sum(counts))) ÷ prod(factorial(big(c)) for c in counts)

function unrank_arrangement(k::Integer, counts::AbstractVector{<:Integer})

    # The k-th arrangement (0-based, lexicographic order) of a multiset with counts[j] copies
    # of label j: a vector with the label of each position (update spec Appendix A). Of the
    # `total` arrangements of what is left, total * counts[j] / remaining start with label j;
    # skip those blocks until k falls inside one. The arithmetic is done in Int128 when
    # total * counts[j] cannot overflow it, which is much faster than BigInt, and in BigInt
    # otherwise.
    counts = collect(counts)
    remaining = sum(counts)
    total = n_arrangements(counts)
    T = total <= typemax(Int128) ÷ max(remaining, 1) ? Int128 : BigInt
    k, total = T(k), T(total)
    arrangement = Vector{Int}(undef, remaining)
    for position in eachindex(arrangement)
        for j in eachindex(counts)
            counts[j] == 0 && continue
            block = total * counts[j] ÷ remaining
            if k < block
                arrangement[position] = j
                counts[j] -= 1
                total = block
                break
            end
            k -= block
        end
        remaining -= 1
    end
    return arrangement
end

function rank_arrangement(arrangement::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer})

    # Position (0-based) of an arrangement in the order unrank_arrangement uses; needed once,
    # to exclude the observed assignment from the draws
    counts = collect(counts)
    remaining = sum(counts)
    total = n_arrangements(counts)
    position = big(0)
    for label in arrangement
        for j in 1:(label - 1)
            position += total * counts[j] ÷ remaining
        end
        total = total * counts[label] ÷ remaining
        counts[label] -= 1
        remaining -= 1
    end
    return position
end

function floyd_sample(M::Integer, B::Integer, rng::AbstractRNG)

    # B distinct integers drawn uniformly from 0:M-1 in exactly B steps (Bentley & Floyd 1987,
    # "Programming pearls: a sample of brilliance", CACM 30(9))
    chosen = Set{BigInt}()
    sample = BigInt[]
    for j in (M - B):(M - 1)
        candidate = rand(rng, big(0):big(j))
        value = candidate in chosen ? big(j) : candidate
        push!(chosen, value)
        push!(sample, value)
    end
    return sample
end

function ri_assignment_draws(L_observed::AbstractVector{Int}, P_observed::AbstractVector{Int}, nperm::Integer,
                             rng::AbstractRNG)

    # Alternative assignments for randomization inference (update spec §7.2): rearrangements of
    # the states' joint (treatment period, policy) labels, so a state's policy and treatment
    # time move together (D14), never equal to the observed assignment. If there are at most
    # `nperm` alternatives they are all enumerated (exact RI); otherwise `nperm` distinct ones
    # are drawn uniformly, as integers by Floyd's algorithm that are then unranked into
    # arrangements (Appendix A). The numbering of arrangements is one-to-one, so uniform
    # integers give uniform assignments.
    # Returns the draws as matrices L and P (one row per draw, one column per state) and
    # whether every alternative was enumerated.
    observed = collect(zip(L_observed, P_observed))
    labels = sort(unique(observed))
    arrangement = [findfirst(==(label), labels) for label in observed]
    counts = [count(==(j), arrangement) for j in eachindex(labels)]
    n_alternatives = n_arrangements(counts) - 1
    observed_rank = rank_arrangement(arrangement, counts)

    exact = n_alternatives <= nperm
    ranks = exact ? collect(big(0):(n_alternatives - 1)) : floyd_sample(n_alternatives, nperm, rng)
    L = Matrix{Int}(undef, length(ranks), length(observed))
    P = Matrix{Int}(undef, length(ranks), length(observed))
    for (b, r) in enumerate(ranks)
        # ranks at or above the observed one are shifted past it
        drawn = unrank_arrangement(r + (r >= observed_rank), counts)
        L[b, :] = [labels[j][1] for j in drawn]
        P[b, :] = [labels[j][2] for j in drawn]
    end
    return L, P, exact
end

### Evaluating assignments (update spec §7.4)

function engine_cell(engine::EngineData, k::Int, L::AbstractMatrix{Int}, active::Union{AbstractMatrix{Bool}, Nothing},
                     use_pre_controls::Bool)

    # ATT(s,t) in cell k = (g, t), for every assignment (row of L) and every state (column)
    # that the assignment treats in period g (update spec §7.4). Controls are the never-treated
    # states and, with use_pre_controls, the states not yet treated at t (§3.4). The control
    # mean for a treated state s is a + bb' xbar_{s,t}: a is the weighted mean of the controls'
    # intercept differences, bb that of their slope differences, so integrating over s's
    # covariates (DID-INT eq. 47) costs one more matrix product. `active` switches states off
    # as treated states and as controls (the jackknife deletions).
    # Returns Tm, which marks the (assignment, state) pairs that have an ATT(s,t), and the
    # ATT(s,t)'s, set to 0 where Tm is false.
    g, t = engine.cohort[k], engine.period[k]
    valid = engine.valid[:, k]'
    controls = (use_pre_controls ? ((L .== NEVER_TREATED) .| (L .> t)) : (L .== NEVER_TREATED)) .& valid
    treated = (L .== g) .& valid
    if !isnothing(active)
        controls = controls .& active
        treated = treated .& active
    end

    Cm = Float64.(controls)
    v = engine.V[:, k]
    WC = Cm * v
    a = (Cm * (v .* engine.Dlambda[:, k])) ./ WC
    ATTst = engine.Dint[:, k]' .- a
    if !isempty(engine.DGamma)
        bb = (Cm * (v .* engine.DGamma[k])) ./ WC
        ATTst .-= bb * engine.Xb[k]'
    end

    # Without a control there is no ATT(s,t)
    Tm = treated .& (WC .> 0)
    return Tm, ifelse.(Tm, ATTst, 0.0)
end

function weighted_mean_of_defined(values::AbstractMatrix{Float64}, weights::AbstractMatrix{Float64})

    # Row by row, the weighted mean of the entries of `values` that are defined (not NaN);
    # NaN for a row without any
    defined = .!isnan.(values)
    w = weights .* defined
    return vec(sum(w .* ifelse.(defined, values, 0.0), dims = 2) ./ sum(w, dims = 2))
end

# Level 2 weights of sub-aggregates (update spec §4.2): their sums of n_t under weighting =
# "att"/"both", equal otherwise
level2_weights(nt::AbstractMatrix{Float64}, weighting::AbstractString) = weighting in ["att", "both"] ? nt : ones(size(nt))

function cell_groups(cells::DataFrame, agg::AbstractString)

    # Group of each (g,t) cell under a cell-based aggregation: its cohort, event time, calendar
    # period or the cell itself (update spec §5.1). For the policy aggregations these are the
    # groups within a policy (§6.3; "policy" has a single one). Returns each cell's group index
    # and the group keys.
    columns = filter(!=(:policy), group_columns(agg))
    cell_keys = [Tuple(cells[k, col] for col in columns) for k in 1:nrow(cells)]
    group_keys = unique(cell_keys)
    return [findfirst(==(key), group_keys) for key in cell_keys], group_keys
end

function grouped_estimates(engine::EngineData, results::DataFrame, agg::AbstractString, weighting::AbstractString,
                           use_pre_controls::Bool, L::AbstractMatrix{Int}, P::AbstractMatrix{Int},
                           active::Union{AbstractMatrix{Bool}, Nothing}, policy_names::AbstractVector,
                           policy_weighting::AbstractString)

    # Estimates of a cell-based aggregation ("simple", "cohort", "time", "calendar" and the
    # policy aggregations) for every assignment (update spec §7.4). An ATT(s,t) belongs to the
    # group of its cell within the policy of its state; the non-policy aggregations are the
    # one-policy case. Per group: num = sum of u ATT(s,t), den = sum of u, nt = sum of n_t.
    # Returns, with one row per assignment:
    #   sub    the sub-aggregates num / den, in the order of the rows of `results` (NaN if undefined)
    #   att_p  each policy's ATT, the Level 2 weighted mean of its defined sub-aggregates
    #   top    the aggregate ATT, or for a policy aggregation the pooled ATT across policies
    B, S = size(L)
    cell_group, group_keys = cell_groups(engine.cells, agg)
    n_groups, n_policies = length(group_keys), length(policy_names)
    num = zeros(B, n_groups, n_policies)
    den = zeros(B, n_groups, n_policies)
    nt = zeros(B, n_groups, n_policies)
    n_att = zeros(B, n_groups, n_policies)
    has_att = falses(B, S)
    in_policy = [P .== p for p in 1:n_policies]
    for k in 1:nrow(engine.cells)
        Tm, ATTst = engine_cell(engine, k, L, active, use_pre_controls)
        has_att .|= Tm
        c = cell_group[k]
        for p in 1:n_policies
            Tp = Tm .& in_policy[p]
            u = Tp .* engine.U[:, k]'
            num[:, c, p] .+= vec(sum(u .* ATTst, dims = 2))
            den[:, c, p] .+= vec(sum(u, dims = 2))
            nt[:, c, p] .+= vec(sum(Tp .* engine.NT[:, k]', dims = 2))
            n_att[:, c, p] .+= vec(sum(Tp, dims = 2))
        end
    end
    sub_groups = ifelse.(den .> 0, num ./ den, NaN)

    # Each policy's ATT_p (§6.3); with one policy, the aggregate ATT (§4.2)
    att_p = Matrix{Float64}(undef, B, n_policies)
    for p in 1:n_policies
        att_p[:, p] = weighted_mean_of_defined(sub_groups[:, :, p], level2_weights(nt[:, :, p], weighting))
    end

    # Pooled ATT across the policies that have an estimate (§6.4), with the weights omega_p
    # recomputed from the assignment's own ATT(s,t)'s: equal, the policy's treated observations,
    # its number of ATT(s,t)'s, or its number of treated states
    if is_policy_agg(agg)
        if policy_weighting == "equal"
            omega = ones(B, n_policies)
        elseif policy_weighting == "treated_obs"
            omega = dropdims(sum(nt, dims = 2), dims = 2)
        elseif policy_weighting == "treated_cells"
            omega = dropdims(sum(n_att, dims = 2), dims = 2)
        else
            omega = Float64.(hcat([vec(sum(has_att .& in_policy[p], dims = 2)) for p in 1:n_policies]...))
        end
        top = weighted_mean_of_defined(att_p, omega)
    else
        top = att_p[:, 1]
    end

    # Sub-aggregates in the order of the results rows
    sub = fill(NaN, B, nrow(results))
    for i in 1:nrow(results)
        key = result_key(results, i, agg)
        p = is_policy_agg(agg) ? findfirst(==(key[1]), policy_names) : 1
        c = findfirst(==(is_policy_agg(agg) ? key[2:end] : key), group_keys)
        (isnothing(p) || isnothing(c)) && continue
        sub[:, i] = sub_groups[:, c, p]
    end
    return (sub = sub, att_p = att_p, top = top)
end

function state_sums(engine::EngineData, L::AbstractMatrix{Int}, active::Union{AbstractMatrix{Bool}, Nothing},
                    use_pre_controls::Bool)

    # For every assignment (row) and state (column): the sums over the state's ATT(s,t)'s of
    # u ATT(s,t), of u and of n_t. A state has ATT(s,t)'s only in the cells of the cohort the
    # assignment gives it, so num / den is its u-weighted mean ATT (the "state" sub-aggregate).
    num = zeros(size(L)); den = zeros(size(L)); nt = zeros(size(L))
    for k in 1:nrow(engine.cells)
        Tm, ATTst = engine_cell(engine, k, L, active, use_pre_controls)
        u = Tm .* engine.U[:, k]'
        num .+= u .* ATTst
        den .+= u
        nt .+= Tm .* engine.NT[:, k]'
    end
    return num, den, nt
end

# Position of the state of each results row in the engine's state order (nothing if absent)
result_states(engine::EngineData, results::DataFrame) = [findfirst(==(s), engine.states) for s in results.state]

# Cell of each results row of agg = "sgt"
result_cells(engine::EngineData, results::DataFrame) =
    [findfirst((engine.cells.gvar .== results.gvar[i]) .& (engine.cells.t .== results.t[i])) for i in 1:nrow(results)]

function state_estimates(engine::EngineData, results::DataFrame, agg::AbstractString, weighting::AbstractString,
                         use_pre_controls::Bool, L::AbstractMatrix{Int}, active::Union{AbstractMatrix{Bool}, Nothing})

    # Estimates of the state-level aggregations for every assignment, each state keeping its
    # own estimates (used by the jackknife; randomization inference uses stand-in states
    # instead, see state_level_ri_estimates). "state": a state's u-weighted mean ATT(s,t);
    # "sgt": the ATT(s,t) itself. Returns the sub-aggregates in the order of the results rows
    # (NaN if undefined) and the aggregate ATT, the Level 2 weighted mean of all defined
    # sub-aggregates.
    B = size(L, 1)
    row_state = result_states(engine, results)
    sub = fill(NaN, B, nrow(results))
    if agg == "state"
        num, den, nt = state_sums(engine, L, active, use_pre_controls)
        sub_states = ifelse.(den .> 0, num ./ den, NaN)
        for i in 1:nrow(results)
            isnothing(row_state[i]) || (sub[:, i] = sub_states[:, row_state[i]])
        end
        return (sub = sub, top = weighted_mean_of_defined(sub_states, level2_weights(nt, weighting)))
    end

    row_cell = result_cells(engine, results)
    top_num = zeros(B); top_den = zeros(B)
    for k in 1:nrow(engine.cells)
        Tm, ATTst = engine_cell(engine, k, L, active, use_pre_controls)
        w = weighting in ["att", "both"] ? Tm .* engine.NT[:, k]' : Float64.(Tm)
        top_num .+= vec(sum(w .* ATTst, dims = 2))
        top_den .+= vec(sum(w, dims = 2))
        for i in 1:nrow(results)
            (row_cell[i] == k && !isnothing(row_state[i])) || continue
            sub[:, i] = ifelse.(Tm[:, row_state[i]], ATTst[:, row_state[i]], NaN)
        end
    end
    return (sub = sub, top = top_num ./ top_den)
end

function engine_estimates(engine::EngineData, results::DataFrame, agg::AbstractString, weighting::AbstractString,
                          use_pre_controls::Bool, L::AbstractMatrix{Int}, P::AbstractMatrix{Int};
                          active::Union{AbstractMatrix{Bool}, Nothing} = nothing,
                          policy_names::AbstractVector = String[DEFAULT_POLICY], policy_weighting::AbstractString = "equal")

    # Every estimate of the chosen aggregation for every assignment: row b of each output
    # belongs to row b of L (treatment periods) and P (policies). At the observed assignment
    # it reproduces the estimates of sub_aggregate_results, aggregate_att and policy_att_results.
    if agg in ["state", "sgt"]
        est = state_estimates(engine, results, agg, weighting, use_pre_controls, L, active)
        return (sub = est.sub, att_p = reshape(est.top, :, 1), top = est.top)
    end
    return grouped_estimates(engine, results, agg, weighting, use_pre_controls, L, P, active, policy_names, policy_weighting)
end

### State-level randomization inference (update spec D9, §7.5)

function draw_stand_ins(L::AbstractMatrix{Int}, has_att::AbstractMatrix{Bool}, cohorts::AbstractVector{Int},
                        rng::AbstractRNG)

    # For each draw (row) and each observed treated state (column, with observed cohort
    # cohorts[j]): one state picked at random among those the draw assigns that cohort and
    # that have an ATT(s,t) in the draw; 0 if there is none. That state stands in for the
    # treated state in every period.
    picks = zeros(Int, size(L, 1), length(cohorts))
    for b in axes(L, 1), (j, g) in enumerate(cohorts)
        candidates = findall(s -> L[b, s] == g && has_att[b, s], axes(L, 2))
        isempty(candidates) || (picks[b, j] = rand(rng, candidates))
    end
    return picks
end

function state_level_ri_estimates(engine::EngineData, results::DataFrame, agg::AbstractString, weighting::AbstractString,
                                  use_pre_controls::Bool, L::AbstractMatrix{Int}, observed_cohort::AbstractDict,
                                  rng::AbstractRNG; picks::Union{AbstractMatrix{Int}, Nothing} = nothing)

    # Placebo statistics of agg = "state" and "sgt" (update spec D9, §7.5). For each observed
    # treated state i with cohort g_i, every draw picks one state at random among those it
    # assigns g_i; that state stands in for i in every period. The statistic is the stand-in's
    # u-weighted mean ATT(s,t) ("state") or its ATT(s,t) in each of i's periods ("sgt"), and the
    # aggregate is the Level 2 weighted mean of the draw's statistics. `picks` fixes the
    # stand-ins instead of drawing them (for tests).
    B = size(L, 1)
    num, den, nt = state_sums(engine, L, nothing, use_pre_controls)
    treated = unique(results.state)
    if isnothing(picks)
        picks = draw_stand_ins(L, den .> 0, [period_code(observed_cohort[s]) for s in treated], rng)
    end
    row_treated = [findfirst(==(s), treated) for s in results.state]

    stats = fill(NaN, B, nrow(results))
    level2 = zeros(B, nrow(results))
    if agg == "state"
        for i in 1:nrow(results), b in 1:B
            s = picks[b, row_treated[i]]
            s == 0 && continue
            stats[b, i] = num[b, s] / den[b, s]
            level2[b, i] = nt[b, s]
        end
    else
        row_cell = result_cells(engine, results)
        for k in 1:nrow(engine.cells)
            Tm, ATTst = engine_cell(engine, k, L, nothing, use_pre_controls)
            for i in 1:nrow(results)
                row_cell[i] == k || continue
                for b in 1:B
                    s = picks[b, row_treated[i]]
                    (s == 0 || !Tm[b, s]) && continue
                    stats[b, i] = ATTst[b, s]
                    level2[b, i] = engine.NT[s, k]
                end
            end
        end
    end
    return (sub = stats, top = weighted_mean_of_defined(stats, level2_weights(level2, weighting)))
end
