# Tests of policies and the policy aggregations.
# The panel is test/data/didintpy_panel.csv: cohorts 2004 (s01, s02), 2006 (s03, s04), 2007 (s05)
# and seven never-treated states.

using Dates
using Random
using Statistics

const POLICY_PANEL = CSV.read(joinpath(@__DIR__, "data", "didintpy_panel.csv"), DataFrame)
const POLICY_TREATED = ["s01", "s02", "s03", "s04", "s05"]
const POLICY_TIMES = [2004, 2004, 2006, 2006, 2007]
# Policy A: s01, s02 (2004) and s05 (2007); policy B: s03, s04 (2006)
const POLICY_LABELS = ["A", "A", "B", "B", "A"]
const POLICY_AGGS = ["policy", "policy_gt", "policy_cohort", "policy_calendar", "policy_time"]

function policy_didint(; data = POLICY_PANEL, treated_states = POLICY_TREATED, treatment_times = POLICY_TIMES, kwargs...)
    return Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        DiDInt.didint("y", "state", "year", data; treated_states = treated_states, treatment_times = treatment_times,
                      nperm = 29, seed = 1234, kwargs...)
    end
end

# Columns of `results` as plain vectors with missing replaced by NaN, for comparisons
plain(results, col) = Float64[ismissing(x) ? NaN : x for x in results[!, col]]
same(a, b) = length(a) == length(b) && all((isnan(x) && isnan(y)) || isapprox(x, y; atol = 1e-10) for (x, y) in zip(a, b))

@testset "Without policy arguments the policy aggregations equal their non-policy counterparts" begin
    # One policy, so policy_cohort = cohort and so on, and the policy's
    # ATT equals the aggregate ATT, with the same inference
    counterparts = [("policy_cohort", "cohort"), ("policy_gt", "simple"), ("policy_calendar", "calendar"), ("policy_time", "time")]
    for (policy_agg, agg) in counterparts, weighting in ["both", "none"], upc in [false, true]
        with_policy = policy_didint(agg = policy_agg, weighting = weighting, use_pre_controls = upc)
        without = policy_didint(agg = agg, weighting = weighting, use_pre_controls = upc)
        @test all(with_policy.policy .== "1")
        suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
        for col in vcat(suffix, DiDInt.inference_columns(suffix))
            @test same(plain(with_policy, col), plain(without, col))
        end
        for (policy_col, agg_col) in zip(vcat("att_p", DiDInt.inference_columns("att_p")), vcat("agg_att", DiDInt.inference_columns("agg_att")))
            @test same(plain(with_policy, policy_col), plain(without, agg_col))
        end
        @test !("agg_att" in names(with_policy))
    end
end

@testset "pvar and treated_policies give the same results; policies do not change other aggregations" begin
    panel = copy(POLICY_PANEL)
    label_of = Dict(zip(POLICY_TREATED, POLICY_LABELS))
    panel.policy_string = [get(label_of, s, "") for s in panel.state]              # "" for never-treated states
    panel.policy_number = [get(Dict("A" => 1, "B" => 2), get(label_of, s, ""), 0) for s in panel.state]   # 0 for never-treated
    by_vector = policy_didint(agg = "policy_cohort", treated_policies = POLICY_LABELS)
    by_column = policy_didint(data = panel, agg = "policy_cohort", pvar = "policy_string")
    @test isequal(by_vector, by_column)
    by_number = policy_didint(data = panel, agg = "policy_cohort", pvar = :policy_number)
    @test by_number.policy == replace.(by_vector.policy, "A" => "1", "B" => "2")
    @test same(plain(by_number, :att_p), plain(by_vector, :att_p))

    # Controls do not depend on policy, so a non-policy aggregation ignores the policies
    for upc in [false, true]
        @test isequal(policy_didint(agg = "sgt", treated_policies = POLICY_LABELS, use_pre_controls = upc),
                      policy_didint(agg = "sgt", use_pre_controls = upc))
    end
end

@testset "A not-yet-treated state of another policy is a control with use_pre_controls only" begin
    # Hand-built design of test_construct_diff_df.jl: A treated 2005 (policy p1), B treated 2007
    # (policy p2), C and D never treated. For A in 2005 and 2006, B is not yet treated.
    policies = Dict("A" => "p1", "B" => "p2")
    for (use_pre_controls, expected_controls) in [(false, ["C", "D"]), (true, ["B", "C", "D"])]
        diff_df, _ = run_construct_diff_df(use_pre_controls)
        att_st = DiDInt.build_att_st(diff_df, "none", "hc1", nothing; policies = policies)
        row = only(eachrow(att_st[(att_st.state .== "A") .& (att_st.t .== Date(2005)), :]))
        @test row.policy == "p1"
        @test row.n_controls == length(expected_controls)
        long_diff(s) = diff_test_lambda(s, Date(2005)) - diff_test_lambda(s, Date(2004))
        @test row.att ≈ long_diff("A") - mean(long_diff.(expected_controls))
        # Once B is treated (2007 on) it is no longer a control for A
        late = only(eachrow(att_st[(att_st.state .== "A") .& (att_st.t .== Date(2007)), :]))
        @test late.n_controls == 2
    end
end

@testset "Each policy's ATT and the pooled ATT under every policy_weighting" begin
    # With weighting = "both" and treated_weight = "treated", ATT_p is the
    # n_t-weighted mean of the policy's ATT(s,t)'s, whatever the policy aggregation
    sgt = policy_didint(agg = "sgt")
    n_t = [count((POLICY_PANEL.state .== sgt.state[i]) .& (POLICY_PANEL.year .== year(sgt.t[i]))) for i in 1:nrow(sgt)]
    label_of = Dict(zip(POLICY_TREATED, POLICY_LABELS))
    row_policy = [label_of[s] for s in sgt.state]
    expected_att_p = Dict(p => sum(n_t[row_policy .== p] .* sgt.att_sgt[row_policy .== p]) / sum(n_t[row_policy .== p]) for p in ["A", "B"])
    expected_omega = Dict(
        "equal" => [0.5, 0.5],
        "treated_obs" => [sum(n_t[row_policy .== p]) for p in ["A", "B"]] ./ sum(n_t),
        "treated_cells" => [count(row_policy .== p) for p in ["A", "B"]] ./ nrow(sgt),
        "treated_states" => [3, 2] ./ 5,
    )
    for agg in POLICY_AGGS, policy_weighting in ["equal", "treated_obs", "treated_cells", "treated_states"]
        results = policy_didint(agg = agg, treated_policies = POLICY_LABELS, policy_aggregation = true, policy_weighting = policy_weighting)
        first_rows = [findfirst(results.policy .== p) for p in ["A", "B"]]
        att_p = Float64.(results.att_p[first_rows])
        @test att_p ≈ [expected_att_p["A"], expected_att_p["B"]] atol = 1e-10
        @test Float64.(results.policy_weights[first_rows]) ≈ expected_omega[policy_weighting] atol = 1e-10
        @test results.att_pool[1] ≈ sum(expected_omega[policy_weighting] .* att_p) atol = 1e-10

        # HC SE of the pool: weighted regression of the two ATT_p's on a constant (one degree of freedom)
        regression = DiDInt.final_regression_results(ones(2, 1), att_p, W = copy(expected_omega[policy_weighting]), hc = "hc1")
        @test results.se_att_pool[1] ≈ regression["beta_hat_se"]
        @test !ismissing(results.jknifese_att_pool[1]) && !ismissing(results.ri_pval_att_pool[1])
    end
end

@testset "Policy argument checks" begin
    # No pooled columns without policy_aggregation; a non-default policy_weighting errors
    results = policy_didint(agg = "policy_cohort", treated_policies = POLICY_LABELS)
    @test !any(occursin("pool", col) || col == "policy_weights" for col in names(results))
    @test_throws r"'policy_weighting' only applies" policy_didint(agg = "policy_cohort", treated_policies = POLICY_LABELS, policy_weighting = "treated_obs")
    # With a non-policy aggregation, policy_aggregation has nothing to pool and is ignored
    @test isequal(policy_didint(agg = "cohort", policy_aggregation = true), policy_didint(agg = "cohort"))
    @test isequal(policy_didint(agg = "cohort", treated_policies = POLICY_LABELS, policy_aggregation = true,
                                policy_weighting = "treated_obs"), policy_didint(agg = "cohort"))
    # With a policy aggregation but no policy inputs there is one policy, so the pool is that policy's ATT
    single = policy_didint(agg = "policy_cohort", policy_aggregation = true)
    @test single.att_pool[1] ≈ single.att_p[1]
    @test_throws r"'policy_weighting' must be one of" policy_didint(agg = "policy", policy_aggregation = true, policy_weighting = "inverse_variance")

    # pvar / treated_policies validation, as for gvar / treatment_times
    panel = copy(POLICY_PANEL)
    label_of = Dict(zip(POLICY_TREATED, POLICY_LABELS))
    panel.policy = [get(label_of, s, "") for s in panel.state]
    panel.gvar = [get(Dict(zip(POLICY_TREATED, POLICY_TIMES)), s, missing) for s in panel.state]
    @test_throws r"either 'pvar' or 'treated_policies'" policy_didint(data = panel, agg = "policy", pvar = "policy", treated_policies = POLICY_LABELS)
    @test_throws r"cannot be used with 'gvar'" policy_didint(data = panel, agg = "policy", gvar = "gvar", treated_states = nothing,
                                                             treatment_times = nothing, treated_policies = POLICY_LABELS)
    @test_throws r"same length" policy_didint(agg = "policy", treated_policies = ["A", "B"])
    @test_throws r"Every treated state needs a policy" policy_didint(agg = "policy", treated_policies = ["A", "A", "", "B", "A"])
    inconsistent = copy(panel)
    inconsistent.policy[findfirst(inconsistent.state .== "s01")] = "B"
    @test_throws r"not constant within state" policy_didint(data = inconsistent, agg = "policy", pvar = "policy")
    no_policy = copy(panel)
    no_policy.policy[no_policy.state .== "s03"] .= ""
    @test_throws r"No policy found" policy_didint(data = no_policy, agg = "policy", pvar = "policy")

    # gvar with pvar works and matches the vector inputs
    @test same(plain(policy_didint(data = panel, agg = "policy", gvar = "gvar", pvar = "policy", treated_states = nothing, treatment_times = nothing), :att_p),
               plain(policy_didint(agg = "policy", treated_policies = POLICY_LABELS), :att_p))
end

@testset "RI draws permute the joint (treatment time, policy) labels" begin
    # Draws are distinct, exclude the observed assignment, and keep the
    # count of every joint label, so a state's policy moves with its treatment time
    states = ["s" * lpad(i, 2, '0') for i in 1:12]
    times = Dict(zip(POLICY_TREATED, Date.(POLICY_TIMES)))
    policies = Dict(zip(POLICY_TREATED, POLICY_LABELS))
    L_observed, P_observed = DiDInt.assignment_row(states, times, policies, Dict("A" => 1, "B" => 2))
    L, P, exact = DiDInt.ri_assignment_draws(L_observed, P_observed, 200, Xoshiro(20260116))
    observed = collect(zip(L_observed, P_observed))
    joint(b) = collect(zip(L[b, :], P[b, :]))
    label_counts(labels) = sort([(label, count(==(label), labels)) for label in unique(labels)])
    @test !exact && size(L, 1) == 200
    @test all(joint(b) != observed for b in 1:200)
    @test all(label_counts(joint(b)) == label_counts(observed) for b in 1:200)
    @test length(unique(joint(b) for b in 1:200)) == 200
end

@testset "The batched engine reproduces the policy estimates" begin
    # The jackknife and RI recompute ATT_p and the pooled ATT with the batched engine; at the
    # observed assignment it must reproduce the regression-based results (hand-built design).
    # More engine tests are in test_engine.jl.
    labels = Dict("A" => Date(2005), "B" => Date(2007))
    policies = Dict("A" => "p1", "B" => "p2")
    for use_pre_controls in [false, true], weighting in ["both", "none"], policy_weighting in ["equal", "treated_obs", "treated_cells", "treated_states"]
        diff_df, ri_diff_df = run_construct_diff_df(use_pre_controls)
        long_diffs = vcat(diff_df, ri_diff_df)
        states = unique(long_diffs.state)
        att_st = DiDInt.build_att_st(diff_df, weighting, "hc1", nothing; policies = policies)
        engine = DiDInt.build_engine(long_diffs, DiDInt.att_st_cells(att_st), states, weighting, "treated")
        L, P = DiDInt.assignment_row(states, labels, policies, Dict("p1" => 1, "p2" => 2))
        policy = (labels = policies, pooled = true, weighting = policy_weighting)
        for agg in POLICY_AGGS
            results = DiDInt.sub_aggregate_results(att_st, agg, weighting, "treated", "hc1"; pooled = true)
            results = DiDInt.policy_att_results(results, att_st, agg, weighting, "hc1", policy)
            est = DiDInt.engine_estimates(engine, results, agg, weighting, use_pre_controls, reshape(L, 1, :), reshape(P, 1, :);
                                          policy_names = ["p1", "p2"], policy_weighting = policy_weighting)
            first_rows = [findfirst(results.policy .== p) for p in ["p1", "p2"]]
            @test vec(est.sub) ≈ Float64.(results[!, DiDInt.SUB_AGG_SPECS[agg].suffix]) atol = 1e-10
            @test vec(est.att_p) ≈ Float64.(results.att_p[first_rows]) atol = 1e-10
            @test est.top[1] ≈ results.att_pool[1] atol = 1e-10
        end
    end
end

@testset "Jackknife for policies needs at least 2 treated states per policy" begin
    # Policy B has a single treated state (s05): no jackknife SE for its ATT, but the other
    # policy's ATT and the pooled ATT still get one
    results = policy_didint(agg = "policy_cohort", treated_policies = ["A", "A", "A", "A", "B"], policy_aggregation = true)
    @test !ismissing(results.jknifese_att_p[findfirst(results.policy .== "A")])
    @test ismissing(results.jknifese_att_p[findfirst(results.policy .== "B")])
    @test !ismissing(results.jknifese_att_pool[1])
end
