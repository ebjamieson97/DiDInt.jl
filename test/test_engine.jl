# Tests of the batched engine behind randomization inference and the jackknife, of the draws
# and of the RI p-values (update spec §7-§9; §12.2 items 7-12).
# The panel is test/data/didintpy_panel.csv: 12 states observed 2001-2008, cohorts 2004 (s01,
# s02), 2006 (s03, s04) and 2007 (s05), seven never-treated states, covariates x1 and x2.

using Dates
using LinearAlgebra
using Random
using Statistics

# The per-draw loops the engine replaced, kept as a slow reference
include("loop_reference.jl")

const ENGINE_PANEL = CSV.read(joinpath(@__DIR__, "data", "didintpy_panel.csv"), DataFrame)
const ENGINE_TREATED = ["s01", "s02", "s03", "s04", "s05"]
const ENGINE_TIMES = Date.([2004, 2004, 2006, 2006, 2007])
# Policy A: s01, s02 (2004) and s05 (2007); policy B: s03, s04 (2006)
const ENGINE_POLICIES = Dict("s01" => "A", "s02" => "A", "s03" => "B", "s04" => "B", "s05" => "A")
const ENGINE_AGGS = ["cohort", "simple", "state", "sgt", "time", "calendar"]
const ENGINE_POLICY_AGGS = ["policy", "policy_gt", "policy_cohort", "policy_calendar", "policy_time"]
const ENGINE_POLICY_WEIGHTINGS = ["equal", "treated_obs", "treated_cells", "treated_states"]

quietly(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

# The first stage of one specification, as didint() computes it (it does not depend on the
# treatment assignment). Integration only applies under ccc = "int" with covariates (D5).
function engine_first_stage(ccc::AbstractString, covariates::Bool)
    data, dates, time_to_index = first_stage_inputs(ENGINE_PANEL)
    lambda_df, vcov_lambda, fit = quietly() do
        DiDInt.iterative_demean(data, ccc, covariates ? ["x1", "x2"] : String[], "hc1", false)
    end
    return (data = data, dates = dates, time_to_index = time_to_index, lambda_df = lambda_df, vcov_lambda = vcov_lambda,
            integration = (ccc == "int" && covariates) ? fit : nothing)
end

# Policy labels and codes as didint() passes them: the policy aggregations use the policies,
# every other aggregation the single default policy
engine_policy_names(agg) = DiDInt.is_policy_agg(agg) ? ["A", "B"] : [DiDInt.DEFAULT_POLICY]
engine_policy_labels(agg, policies) = DiDInt.is_policy_agg(agg) ? policies : Dict(s => DiDInt.DEFAULT_POLICY for s in keys(policies))

# The regression-based second stage for an assignment (treated state => treatment time and
# policy): rebuilds the long differences and the ATT(s,t) table, then aggregates. This is the
# path didint() uses for its point estimates.
function second_stage(stage, times::Dict, policies::Dict, agg, weighting, treated_weight, upc; policy_weighting = "equal")
    treated = sort(collect(keys(times)))
    labels = engine_policy_labels(agg, policies)
    return quietly() do
        diff_df, ri_diff_df = DiDInt.construct_diff_df(stage.lambda_df, treated, [times[s] for s in treated], stage.time_to_index,
                                                        stage.dates, upc, stage.data)
        att_st = DiDInt.build_att_st(diff_df, weighting, "hc1", stage.vcov_lambda; integration = stage.integration, policies = labels)
        results = DiDInt.sub_aggregate_results(att_st, agg, weighting, treated_weight, "hc1"; pooled = DiDInt.is_policy_agg(agg))
        if DiDInt.is_policy_agg(agg)
            results = DiDInt.policy_att_results(results, att_st, agg, weighting, "hc1",
                                                (labels = labels, pooled = true, weighting = policy_weighting))
        else
            results = DiDInt.aggregate_att(results, agg, weighting, "hc1")
        end
        (results = results, att_st = att_st, diff_df = diff_df, long_diffs = vcat(diff_df, ri_diff_df))
    end
end

# Top-level estimates of a results table: the aggregate ATT or the pooled ATT, and each policy's ATT
top_estimate(results, agg) = DiDInt.is_policy_agg(agg) ? results.att_pool[1] : results.agg_att[1]
policy_estimates(results) = Float64[results.att_p[findfirst(results.policy .== p)] for p in ["A", "B"]]

# The engine at one assignment, with its sub-aggregates in the order of the rows of `results`
function engine_at(engine, results, times::Dict, policies::Dict, agg, weighting, upc; policy_weighting = "equal")
    names = engine_policy_names(agg)
    code = Dict(name => i for (i, name) in enumerate(names))
    L, P = DiDInt.assignment_row(engine.states, times, engine_policy_labels(agg, policies), code)
    return DiDInt.engine_estimates(engine, results, agg, weighting, upc, reshape(L, 1, :), reshape(P, 1, :);
                                   policy_names = names, policy_weighting = policy_weighting)
end

close_enough(a, b) = isapprox(a, b; atol = 1e-10)

const ENGINE_OBSERVED_TIMES = Dict(zip(ENGINE_TREATED, ENGINE_TIMES))

@testset "The engine reproduces every observed estimate (§7.4)" begin
    # Spec §12.2 item 7, first half: at the observed assignment the engine gives every
    # sub-aggregate, the aggregate ATT, each policy's ATT and the pooled ATT, for every agg,
    # weighting, treated_weight, use_pre_controls and ccc
    for (ccc, covariates) in [("int", false), ("int", true), ("hom", true), ("time", true)]
        stage = engine_first_stage(ccc, covariates)
        for upc in [false, true], weighting in ["both", "att", "diff", "none"], treated_weight in ["treated", "population"]
            cases = vcat([(agg, "equal") for agg in ENGINE_AGGS],
                         [(agg, pw) for agg in ENGINE_POLICY_AGGS for pw in ENGINE_POLICY_WEIGHTINGS])
            for (agg, policy_weighting) in cases
                observed = second_stage(stage, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, agg, weighting, treated_weight, upc;
                                        policy_weighting = policy_weighting)
                engine = DiDInt.build_engine(observed.long_diffs, DiDInt.att_st_cells(observed.att_st), unique(observed.long_diffs.state),
                                             weighting, treated_weight; integration = stage.integration)
                est = engine_at(engine, observed.results, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, agg, weighting, upc;
                                policy_weighting = policy_weighting)
                suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
                @test all(close_enough.(vec(est.sub), Float64.(observed.results[!, suffix])))
                @test close_enough(est.top[1], top_estimate(observed.results, agg))
                if DiDInt.is_policy_agg(agg)
                    @test all(close_enough.(vec(est.att_p), policy_estimates(observed.results)))
                end
            end
        end
    end
end

@testset "The engine matches a rebuilt ATT(s,t) table and the loop reference on random draws" begin
    # Spec §12.2 item 7, second half: for 20 random assignments, the engine built from the
    # observed assignment gives (a) what the regression-based second stage gives when the
    # ATT(s,t) table is rebuilt for the drawn assignment, and (b) what the per-draw loops it
    # replaced give
    names = ["A", "B"]
    code = Dict("A" => 1, "B" => 2)
    for (ccc, covariates) in [("int", false), ("int", true)]
        stage = engine_first_stage(ccc, covariates)
        for upc in [false, true], (weighting, treated_weight) in [("both", "treated"), ("both", "population"), ("att", "treated"), ("none", "treated")]
            observed = second_stage(stage, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, "sgt", weighting, treated_weight, upc)
            cells = DiDInt.att_st_cells(observed.att_st)
            states = unique(observed.long_diffs.state)
            engine = DiDInt.build_engine(observed.long_diffs, cells, states, weighting, treated_weight; integration = stage.integration)
            lookup = DiDInt.long_difference_lookup(observed.long_diffs)
            L_observed, P_observed = DiDInt.assignment_row(states, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, code)
            L, P, exact = DiDInt.ri_assignment_draws(L_observed, P_observed, 20, Xoshiro(20260501))
            @test !exact && size(L) == (20, length(states))

            for b in 1:20
                drawn = findall(!=(DiDInt.NEVER_TREATED), L[b, :])
                times = Dict(states[s] => Date(Dates.UTD(L[b, s])) for s in drawn)
                policies = Dict(states[s] => names[P[b, s]] for s in drawn)
                for agg in vcat(ENGINE_AGGS, ENGINE_POLICY_AGGS)
                    policy_weighting = agg == "policy_cohort" ? "treated_states" : (agg == "policy_time" ? "treated_cells" : "treated_obs")
                    suffix = DiDInt.SUB_AGG_SPECS[agg].suffix

                    # (a) the regression path on the ATT(s,t) table rebuilt for this assignment
                    rebuilt = second_stage(stage, times, policies, agg, weighting, treated_weight, upc; policy_weighting = policy_weighting)
                    est = engine_at(engine, rebuilt.results, times, policies, agg, weighting, upc; policy_weighting = policy_weighting)
                    @test all(close_enough.(vec(est.sub), Float64.(rebuilt.results[!, suffix])))
                    @test close_enough(est.top[1], top_estimate(rebuilt.results, agg))
                    if DiDInt.is_policy_agg(agg)
                        @test all(close_enough.(vec(est.att_p), policy_estimates(rebuilt.results)))
                    end

                    # (b) the per-draw loop
                    loop = att_st_point_estimates(cells, lookup, times, states, weighting, upc; integration = stage.integration,
                                                  policies = engine_policy_labels(agg, policies))
                    sub = point_sub_aggregates(loop, agg, weighting, treated_weight)
                    @test all(close_enough(est.sub[1, i], sub[DiDInt.result_key(rebuilt.results, i, agg)].att) for i in 1:nrow(rebuilt.results))
                    if DiDInt.is_policy_agg(agg)
                        att_p, pooled = point_policy_aggregates(sub, loop, names, weighting, policy_weighting)
                        @test all(close_enough.(vec(est.att_p), att_p))
                        @test close_enough(est.top[1], pooled)
                    else
                        @test close_enough(est.top[1], point_aggregate(sub, weighting))
                    end
                end
            end
        end
    end
end

@testset "Ranking and unranking arrangements; Floyd sampling (Appendix A)" begin
    # Counts (3, 1, 1): 5!/3! = 20 arrangements, numbered in lexicographic order
    counts = [3, 1, 1]
    @test DiDInt.n_arrangements(counts) == 20
    arrangements = [DiDInt.unrank_arrangement(k, counts) for k in 0:19]
    @test arrangements[1] == [1, 1, 1, 2, 3] && arrangements[end] == [3, 2, 1, 1, 1]
    @test issorted(arrangements) && allunique(arrangements)
    @test all(sort(a) == [1, 1, 1, 2, 3] for a in arrangements)
    @test all(DiDInt.rank_arrangement(arrangements[k + 1], counts) == k for k in 0:19)

    # Counts too large for Int64 factorials still round-trip
    big_counts = [30, 10, 6, 4]
    total = DiDInt.n_arrangements(big_counts)
    @test total == factorial(big(50)) ÷ (factorial(big(30)) * factorial(big(10)) * factorial(big(6)) * factorial(big(4)))
    for k in [big(0), total ÷ 3, total - 1]
        arrangement = DiDInt.unrank_arrangement(k, big_counts)
        @test [count(==(j), arrangement) for j in 1:4] == big_counts
        @test DiDInt.rank_arrangement(arrangement, big_counts) == k
    end

    # Floyd: distinct integers in range; drawing all of them returns every integer once
    rng = Xoshiro(7)
    sample = DiDInt.floyd_sample(big(1000), 50, rng)
    @test length(sample) == 50 && allunique(sample) && all(0 .<= sample .< 1000)
    @test sort(DiDInt.floyd_sample(big(12), 12, rng)) == collect(0:11)
    # Each of the 6 pairs from 0:3 is drawn about equally often
    pair_counts = Dict{Vector{BigInt}, Int}()
    for _ in 1:6000
        pair = sort(DiDInt.floyd_sample(big(4), 2, rng))
        pair_counts[pair] = get(pair_counts, pair, 0) + 1
    end
    @test length(pair_counts) == 6 && all(850 .<= values(pair_counts) .<= 1150)
end

@testset "RI draws: exact enumeration and distinct samples of the joint labels (§7.2)" begin
    # Spec §12.2 item 8. Five states; state 1 treated in period 10 with policy 1, state 2 in
    # period 20 with policy 2: 5!/3! = 20 arrangements, 19 alternatives
    L_observed = [10, 20, 0, 0, 0]
    P_observed = [1, 2, 0, 0, 0]
    joint_counts(L_row, P_row) = sort([(label, count(==(label), collect(zip(L_row, P_row)))) for label in unique(collect(zip(L_row, P_row)))])
    L, P, exact = DiDInt.ri_assignment_draws(L_observed, P_observed, 100, Xoshiro(1))
    @test exact && size(L, 1) == 19
    rows = [(L[b, :], P[b, :]) for b in 1:19]
    @test allunique(rows) && !((L_observed, P_observed) in rows)
    @test all(joint_counts(row...) == joint_counts(L_observed, P_observed) for row in rows)
    # Exactly as many alternatives as nperm are still enumerated; one fewer is sampled
    @test DiDInt.ri_assignment_draws(L_observed, P_observed, 19, Xoshiro(1))[3]
    L, P, exact = DiDInt.ri_assignment_draws(L_observed, P_observed, 18, Xoshiro(1))
    @test !exact && size(L, 1) == 18 && allunique([(L[b, :], P[b, :]) for b in 1:18])

    # A policy never separates from its treatment time: two states share period 10 but differ in
    # policy, so there are 4!/2! = 12 arrangements rather than 4!/(2! 2!) = 6
    L, P, exact = DiDInt.ri_assignment_draws([10, 10, 0, 0], [1, 2, 0, 0], 100, Xoshiro(1))
    @test exact && size(L, 1) == 11
    @test all(sort(collect(zip(L[b, :], P[b, :]))) == [(0, 0), (0, 0), (10, 1), (10, 2)] for b in 1:11)

    # Sampled draws from the panel's assignment: 12!/(7! 2! 2! 1!) = 23760 arrangements with
    # these policies (s05 has its own cohort), distinct and never the observed one
    states = sort(unique(ENGINE_PANEL.state))
    L_observed, P_observed = DiDInt.assignment_row(states, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, Dict("A" => 1, "B" => 2))
    L, P, exact = DiDInt.ri_assignment_draws(L_observed, P_observed, 200, Xoshiro(20260116))
    rows = [(L[b, :], P[b, :]) for b in 1:200]
    @test !exact && allunique(rows) && !((L_observed, P_observed) in rows)
    @test all(joint_counts(row...) == joint_counts(L_observed, P_observed) for row in rows)
    # The same seed gives the same draws
    @test DiDInt.ri_assignment_draws(L_observed, P_observed, 200, Xoshiro(20260116))[1] == L
end

@testset "Kernel-smoothed RI p-value (§8.2)" begin
    # Spec §12.2 item 9: reference numbers from didintpy's kernel_smoothed_ri_pvalue
    placebos = [-0.21, 0.05, 0.33, -0.12, 0.47, -0.38, 0.09, 0.26, -0.02, 0.15, -0.29, 0.41, -0.07]
    @test std(placebos) ≈ 0.2671 atol = 5e-5
    @test 1.575 * std(placebos) * 13^(-4 / 9) ≈ 0.1346 atol = 5e-5
    @test DiDInt.ri_pvalue_kernel(placebos, 0.44) ≈ 0.1412 atol = 5e-5
    @test DiDInt.ri_pvalue(placebos, 0.44) ≈ 2 / 14
    # The smoothed p-value moves smoothly with T, where the ordinary one jumps
    smoothed = [DiDInt.ri_pvalue_kernel(placebos, T) for T in 0.40:0.02:0.50]
    @test all(isapprox.(smoothed, [0.191, 0.165, 0.141, 0.120, 0.100, 0.083]; atol = 5e-4))
    @test DiDInt.ri_pvalue_kernel(placebos, -0.44) == DiDInt.ri_pvalue_kernel(placebos, 0.44)
    @test DiDInt.ri_pvalue_kernel([placebos; NaN], 0.44) == DiDInt.ri_pvalue_kernel(placebos, 0.44)

    # Fallbacks to the ordinary p-value, ties counted as hits (spec §12.2 item 10): a single
    # placebo, or placebos without any spread
    @test DiDInt.ri_pvalue_kernel([0.3], 0.3) == (1 + 1) / (1 + 1)
    @test DiDInt.ri_pvalue_kernel([0.3, 0.3, 0.3], 0.1 + 0.2) == (1 + 3) / (3 + 1)
    @test DiDInt.ri_pvalue_kernel([0.3, 0.3, 0.3], 0.5) == (1 + 0) / (3 + 1)
    @test ismissing(DiDInt.ri_pvalue_kernel([NaN], 0.3))

    # With every alternative enumerated: kernel below RI_KERNEL_CUTOFF placebos, ordinary from it on
    @test DiDInt.RI_KERNEL_CUTOFF == 99
    @test DiDInt.ri_pvalue_enumerated(placebos, 0.44) == DiDInt.ri_pvalue_kernel(placebos, 0.44)
    many = collect(range(-1, 1, length = 99))
    @test DiDInt.ri_pvalue_enumerated(many, 0.44) == DiDInt.ri_pvalue(many, 0.44)
    @test DiDInt.ri_pvalue_enumerated(many[1:98], 0.44) == DiDInt.ri_pvalue_kernel(many[1:98], 0.44)
end

@testset "State-level RI: the engine's statistics match the loop with the same stand-ins (§7.5)" begin
    # Spec §12.2 item 11, on the design without covariates
    stage = engine_first_stage("int", false)
    for upc in [false, true], weighting in ["both", "none"], agg in ["state", "sgt"]
        observed = second_stage(stage, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, agg, weighting, "treated", upc)
        results = observed.results
        cells = DiDInt.att_st_cells(observed.att_st)
        states = unique(observed.long_diffs.state)
        engine = DiDInt.build_engine(observed.long_diffs, cells, states, weighting, "treated")
        lookup = DiDInt.long_difference_lookup(observed.long_diffs)
        L_observed, P_observed = DiDInt.assignment_row(states, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, Dict("A" => 1, "B" => 1))
        L, _, _ = DiDInt.ri_assignment_draws(L_observed, P_observed, 30, Xoshiro(99))

        # Stand-ins: one state per draw and treated state, among the states the draw gives its cohort
        treated = unique(results.state)
        cohorts = [DiDInt.period_code(ENGINE_OBSERVED_TIMES[s]) for s in treated]
        picks = DiDInt.draw_stand_ins(L, trues(size(L)), cohorts, Xoshiro(5))
        @test all(L[b, picks[b, j]] == cohorts[j] for b in 1:30, j in eachindex(treated))
        @test length(unique(picks[:, 1])) > 1

        est = DiDInt.state_level_ri_estimates(engine, results, agg, weighting, upc, L, ENGINE_OBSERVED_TIMES, Xoshiro(5); picks = picks)
        result_keys = [DiDInt.result_key(results, i, agg) for i in 1:nrow(results)]
        for b in 1:30
            drawn = findall(!=(DiDInt.NEVER_TREATED), L[b, :])
            times = Dict(states[s] => Date(Dates.UTD(L[b, s])) for s in drawn)
            loop = att_st_point_estimates(cells, lookup, times, states, weighting, upc)
            stand_ins = Dict{String, Union{String, Nothing}}(s => states[picks[b, j]] for (j, s) in enumerate(treated))
            stats, aggregate = state_level_ri_draw(loop, result_keys, ENGINE_OBSERVED_TIMES, agg, weighting, "treated"; picks = stand_ins)
            @test all(close_enough.(est.sub[b, :], stats))
            @test close_enough(est.top[b], aggregate)
        end

        # Drawing the stand-ins inside gives the same result as passing the same draws in
        same_draws = DiDInt.state_level_ri_estimates(engine, results, agg, weighting, upc, L, ENGINE_OBSERVED_TIMES, Xoshiro(5))
        @test isequal(same_draws.sub, est.sub) && isequal(same_draws.top, est.top)
    end
end

@testset "State-level RI (sgt) uses one stand-in state per treated state in every period" begin
    # s01's sgt statistics in a draw are the ATT(s,t)'s of a single state that the draw treats
    # in 2004, in all of s01's periods (never a mix of the two states treated in 2004)
    stage = engine_first_stage("int", false)
    observed = second_stage(stage, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, "sgt", "both", "treated", false)
    results = observed.results
    states = unique(observed.long_diffs.state)
    engine = DiDInt.build_engine(observed.long_diffs, DiDInt.att_st_cells(observed.att_st), states, "both", "treated")
    L_observed, P_observed = DiDInt.assignment_row(states, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, Dict("A" => 1, "B" => 1))
    L, P, _ = DiDInt.ri_assignment_draws(L_observed, P_observed, 50, Xoshiro(3))
    est = DiDInt.state_level_ri_estimates(engine, results, "sgt", "both", false, L, ENGINE_OBSERVED_TIMES, Xoshiro(11))
    own = DiDInt.engine_estimates(engine, results, "sgt", "both", false, L, P)
    s01_rows = findall(results.state .== "s01")
    g = DiDInt.period_code(Date(2004))
    for b in 1:50
        # each state the draw treats in 2004, with its own ATT(s,t)'s in s01's periods
        candidates = [states[s] for s in findall(==(g), L[b, :])]
        candidate_stats = map(candidates) do c
            relabelled = copy(results)
            relabelled.state[s01_rows] .= c
            DiDInt.engine_estimates(engine, relabelled, "sgt", "both", false, L[b:b, :], P[b:b, :]).sub[1, s01_rows]
        end
        @test any(stats -> all(close_enough.(est.sub[b, s01_rows], stats)), candidate_stats)
    end
    @test size(own.sub) == size(est.sub)
end

@testset "Batched jackknife equals the delete-one loop (§9.2)" begin
    # Spec §12.2 item 12: for ccc = "int" and "state" the engine's row h is the estimate without
    # state h, as the loop over deletions computed it
    names = ["A", "B"]
    for (ccc, covariates) in [("int", true), ("state", true), ("int", false)], upc in [false, true], weighting in ["both", "none"]
        stage = engine_first_stage(ccc, covariates)
        for agg in vcat(ENGINE_AGGS, ENGINE_POLICY_AGGS)
            observed = second_stage(stage, ENGINE_OBSERVED_TIMES, ENGINE_POLICIES, agg, weighting, "treated", upc; policy_weighting = "treated_obs")
            results = observed.results
            cells = DiDInt.att_st_cells(observed.att_st)
            states = unique(observed.diff_df.state)
            n_states = length(states)
            engine = DiDInt.build_engine(observed.diff_df, cells, states, weighting, "treated"; integration = stage.integration)
            policy_names = engine_policy_names(agg)
            labels = engine_policy_labels(agg, ENGINE_POLICIES)
            L, P = DiDInt.assignment_row(states, ENGINE_OBSERVED_TIMES, labels, Dict(name => i for (i, name) in enumerate(policy_names)))
            est = DiDInt.engine_estimates(engine, results, agg, weighting, upc, repeat(L', n_states), repeat(P', n_states);
                                          active = .!Matrix{Bool}(I, n_states, n_states), policy_names = policy_names,
                                          policy_weighting = "treated_obs")

            lookup = DiDInt.long_difference_lookup(observed.diff_df)
            for (h, without) in enumerate(states)
                times = filter(label -> label.first != without, ENGINE_OBSERVED_TIMES)
                loop = att_st_point_estimates(cells, lookup, times, filter(!=(without), states), weighting, upc;
                                              integration = stage.integration, policies = labels)
                sub = point_sub_aggregates(loop, agg, weighting, "treated")
                for i in 1:nrow(results)
                    key = DiDInt.result_key(results, i, agg)
                    @test haskey(sub, key) ? close_enough(est.sub[h, i], sub[key].att) : isnan(est.sub[h, i])
                end
                if DiDInt.is_policy_agg(agg)
                    att_p, pooled = point_policy_aggregates(sub, loop, names, weighting, "treated_obs")
                    @test all(close_enough.(est.att_p[h, :], att_p))
                    @test close_enough(est.top[h], pooled)
                else
                    @test close_enough(est.top[h], point_aggregate(sub, weighting))
                end
            end
        end
    end
end

@testset "didint(): exact and sampled RI, ri_procedure and nperm" begin
    run(data, cohorts; kwargs...) = quietly() do
        treated = sort(collect(keys(cohorts)))
        DiDInt.didint("y", "state", "year", data; treated_states = treated, treatment_times = [cohorts[s] for s in treated],
                      seed = 1234, verbose = false, kwargs...)
    end

    # Six states, one treated: 5 alternatives, all enumerated, kernel-smoothed p-values. The
    # p-value equals the kernel formula applied to the estimates obtained by treating each of the
    # other states instead (through didint() itself).
    data = small_panel(cohorts = Dict("s1" => 2004), n_states = 6, effect = 1.0, noise = 0.6)
    for upc in [false, true]
        exact = run(data, Dict("s1" => 2004); agg = "cohort", nperm = 999, use_pre_controls = upc)
        @test all(exact.ri_procedure .== "kernel_exact")
        @test exact.nperm[1] == 5
        placebos = [run(data, Dict("s$i" => 2004); agg = "cohort", nperm = 2, use_pre_controls = upc).agg_att[1] for i in 2:6]
        @test exact.ri_pval_agg_att[1] ≈ DiDInt.ri_pvalue_kernel(placebos, exact.agg_att[1]) atol = 1e-10
        @test exact.ri_pval_att_cohort[1] ≈ exact.ri_pval_agg_att[1] atol = 1e-10
        # Exact RI does not depend on the seed or on asking for more permutations
        other = quietly() do
            DiDInt.didint("y", "state", "year", data; treated_states = ["s1"], treatment_times = [2004], agg = "cohort",
                          nperm = 50, seed = 99, verbose = false, use_pre_controls = upc)
        end
        @test other.ri_pval_agg_att[1] ≈ exact.ri_pval_agg_att[1] atol = 1e-12
    end

    # Fewer permutations than alternatives: sampled, ordinary p-value on the grid 1/(nperm + 1)
    sampled = run(data, Dict("s1" => 2004); agg = "cohort", nperm = 3)
    @test all(sampled.ri_procedure .== "sampled") && sampled.nperm[1] == 3
    @test sampled.ri_pval_agg_att[1] in [k / 4 for k in 1:4]

    # Ten states, four treated in one period: 210 arrangements, 209 alternatives. They are all
    # enumerated under nperm = 999 with the ordinary p-value, as there are at least 99.
    common = Dict("s1" => 2004, "s2" => 2004, "s3" => 2004, "s4" => 2004)
    data = small_panel(cohorts = common, n_states = 10, effect = 1.0, noise = 0.6)
    for agg in ["cohort", "simple", "state", "sgt", "time", "calendar", "policy_cohort"]
        results = run(data, common; agg = agg, nperm = 999)
        @test all(results.ri_procedure .== "ordinary_exact") && results.nperm[1] == 209
        suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
        # every level has an RI p-value, on the grid 1/210
        top = agg == "policy_cohort" ? results.ri_pval_att_p[1] : results.ri_pval_agg_att[1]
        for p in vcat(top, results[!, "ri_pval_" * suffix])
            @test isapprox(p * 210, round(p * 210); atol = 1e-8) && 1 / 210 <= p <= 1
        end
    end
    # With every alternative enumerated the aggregate's p-value is the rank of |observed| among
    # all 210 assignments, here computed by running didint() on each of them
    results = run(data, common; agg = "cohort", nperm = 999)
    all_states = ["s$i" for i in 1:10]
    estimates = Float64[]
    for a in 1:7, b in (a + 1):8, c in (b + 1):9, d in (c + 1):10
        cohorts = Dict(all_states[i] => 2004 for i in (a, b, c, d))
        push!(estimates, run(data, cohorts; agg = "cohort", nperm = 2).agg_att[1])
    end
    @test length(estimates) == 210
    @test results.ri_pval_agg_att[1] ≈ count(abs.(estimates) .>= abs(results.agg_att[1]) * (1 - 1e-8)) / 210 atol = 1e-12

    # The same seed gives the same results; didint() leaves the global random number stream alone
    Random.seed!(42); before = rand()
    Random.seed!(42)
    first_run = run(data, common; agg = "sgt", nperm = 49)
    @test rand() == before
    @test isequal(first_run, run(data, common; agg = "sgt", nperm = 49))
    @test all(first_run.ri_procedure .== "sampled")
end

@testset "truejack still refits the first stage" begin
    # Spec §12.2 item 12, second half: under ccc = "hom" with covariates a deletion changes the
    # other states' lambdas, so the refit jackknife differs from the fast one; under "int" the
    # argument has no effect
    run(; kwargs...) = quietly() do
        DiDInt.didint("y", "state", "year", ENGINE_PANEL; treated_states = ENGINE_TREATED, treatment_times = [2004, 2004, 2006, 2006, 2007],
                      covariates = ["x1", "x2"], nperm = 9, seed = 1234, verbose = false, kwargs...)
    end
    for agg in ["cohort", "state", "policy_cohort"]
        top = agg == "policy_cohort" ? "jknifese_att_p" : "jknifese_agg_att"
        fast = run(ccc = "hom", agg = agg, truejack = false)
        refit = run(ccc = "hom", agg = agg, truejack = true)
        @test isfinite(refit[1, top]) && isfinite(fast[1, top])
        @test abs(refit[1, top] - fast[1, top]) > 1e-6
        @test isequal(run(ccc = "int", agg = agg, truejack = true)[!, top], run(ccc = "int", agg = agg, truejack = false)[!, top])
    end
end
