# Tests of the ATT(s,t)-first second stage (update spec §3.3-§5, §12.2).
# Uses a small deterministic panel: states s1..s8 observed 2001-2006 with cell sizes that vary by
# state and period; y = state effect + period effect + treatment effect + within-cell noise.

using Dates

const SECOND_STAGE_AGGS = ["cohort", "simple", "state", "sgt", "time", "calendar"]

function small_panel(; cohorts::Dict{String, Int}, n_states::Int = 8, years = 2001:2006,
                     effect::Float64 = 0.0, noise::Float64 = 0.0, shock = nothing)
    state_col = String[]; year_col = Int[]; y_col = Float64[]
    for si in 1:n_states, (k, yr) in enumerate(years)
        s = "s$si"
        g = get(cohorts, s, nothing)
        n_st = 2 + mod(si + 2 * k, 4) + (si == 1 && k >= 4 ? 5 : 0)
        for j in 1:n_st
            y = si + 0.5 * k^2 + noise * sin(37 * si + 11 * k + 3 * j)
            if !isnothing(g) && yr >= g
                y += effect * (1 + 0.5 * (yr - g))
            end
            if !isnothing(shock) && s == shock.state && yr == shock.year
                y += shock.size
            end
            push!(state_col, s); push!(year_col, yr); push!(y_col, y)
        end
    end
    return DataFrame(state = state_col, year = year_col, y = y_col)
end

function quiet_didint(data, cohorts::Dict{String, Int}; kwargs...)
    treated = sort(collect(keys(cohorts)))
    return Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        DiDInt.didint("y", "state", "year", data; treated_states = treated,
                      treatment_times = [cohorts[s] for s in treated], nperm = 19, seed = 1234,
                      verbose = false, kwargs...)
    end
end

const STAG_COHORTS = Dict("s1" => 2003, "s2" => 2003, "s3" => 2005, "s4" => 2005)
const COMMON_COHORTS = Dict("s1" => 2004, "s2" => 2004, "s3" => 2004)

@testset "The batched engine reproduces the regression-based estimates" begin
    # The jackknife and RI recompute estimates with the batched engine; at the observed
    # assignment it must reproduce every sub-aggregate and aggregate (update spec §7.4).
    # Uses the hand-built design of test_construct_diff_df.jl (A treated 2005, B 2007).
    # More engine tests are in test_engine.jl.
    labels = Dict("A" => Date(2005), "B" => Date(2007))
    policies = Dict("A" => DiDInt.DEFAULT_POLICY, "B" => DiDInt.DEFAULT_POLICY)
    for use_pre_controls in [false, true]
        diff_df, ri_diff_df = run_construct_diff_df(use_pre_controls)
        long_diffs = vcat(diff_df, ri_diff_df)
        states = unique(long_diffs.state)
        for weighting in ["both", "att", "diff", "none"], treated_weight in ["treated", "population"]
            att_st = DiDInt.build_att_st(diff_df, weighting, "hc1", nothing)
            engine = DiDInt.build_engine(long_diffs, DiDInt.att_st_cells(att_st), states, weighting, treated_weight)
            L, P = DiDInt.assignment_row(states, labels, policies, Dict(DiDInt.DEFAULT_POLICY => 1))
            for agg in SECOND_STAGE_AGGS
                results = DiDInt.aggregate_att(DiDInt.sub_aggregate_results(att_st, agg, weighting, treated_weight, "hc1"),
                                               agg, weighting, "hc1")
                est = DiDInt.engine_estimates(engine, results, agg, weighting, use_pre_controls, reshape(L, 1, :), reshape(P, 1, :))
                suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
                @test vec(est.sub) ≈ Float64.(results[!, suffix]) atol = 1e-10
                @test est.top[1] ≈ results.agg_att[1] atol = 1e-10
            end
        end
    end
end

@testset "Zero effect: parallel trends, no noise, unequal cell sizes" begin
    # Spec §12.2 item 3: every agg returns 0 for every weighting and use_pre_controls
    data = small_panel(cohorts = STAG_COHORTS)
    for agg in SECOND_STAGE_AGGS, weighting in ["both", "att", "diff", "none"], upc in [false, true]
        results = quiet_didint(data, STAG_COHORTS; agg = agg, weighting = weighting, use_pre_controls = upc)
        suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
        @test all(abs.(results[!, suffix]) .< 1e-10)
        @test abs(results.agg_att[1]) < 1e-10
    end
end

@testset "treated_weight: treated = sum n_t ATT(s,t) / sum n_t, population = sum n ATT(s,t) / sum n" begin
    # Spec §12.2 item 4: one cohort, two treated states with different n_{s,r1}/n_{s,t}
    cohorts = Dict("s1" => 2004, "s2" => 2004)
    data = small_panel(cohorts = cohorts, n_states = 5, effect = 1.0, noise = 0.3)
    cell_n(s, yr) = count((data.state .== s) .& (data.year .== yr))
    sgt = quiet_didint(data, cohorts; agg = "sgt")
    n_t = [cell_n(sgt.state[i], year(sgt.t[i])) for i in 1:nrow(sgt)]
    n = n_t .+ [cell_n(sgt.state[i], 2003) for i in 1:nrow(sgt)]
    @test n_t[sgt.state .== "s1"] ./ n[sgt.state .== "s1"] != n_t[sgt.state .== "s2"] ./ n[sgt.state .== "s2"]
    for tw in ["treated", "population"]
        u = tw == "treated" ? n_t : n
        cohort = quiet_didint(data, cohorts; agg = "cohort", treated_weight = tw)
        @test cohort.att_cohort[1] ≈ sum(u .* sgt.att_sgt) / sum(u)
        simple = quiet_didint(data, cohorts; agg = "simple", treated_weight = tw)
        for i in 1:nrow(simple)
            in_cell = sgt.t .== simple.time[i]
            @test simple.att_gt[i] ≈ sum(u[in_cell] .* sgt.att_sgt[in_cell]) / sum(u[in_cell])
        end
    end
    # With one treated state per cell the two agree for ATT(g,t). (Across periods they still
    # differ, e.g. for agg = "cohort", because n_{s,r1}/n_{s,t} changes with t.)
    one = Dict("s1" => 2004)
    data_one = filter(r -> r.state != "s2", data)
    @test quiet_didint(data_one, one; agg = "simple", treated_weight = "treated").att_gt ≈
          quiet_didint(data_one, one; agg = "simple", treated_weight = "population").att_gt
end

@testset "Consistency: under att/both with treated_weight = treated, agg_att is the same for every agg" begin
    # Spec §5.3 / §12.2 item 5: DID-INT eq. 49 with the eq. 84 weights
    data = small_panel(cohorts = STAG_COHORTS, effect = 1.0, noise = 0.3)
    for weighting in ["att", "both"], upc in [false, true]
        agg_atts = [quiet_didint(data, STAG_COHORTS; agg = agg, weighting = weighting, use_pre_controls = upc).agg_att[1]
                    for agg in SECOND_STAGE_AGGS]
        @test all(isapprox.(agg_atts, agg_atts[1]; atol = 1e-10))
    end
end

@testset "SE fallbacks: single ATT(s,t) groups (A1) and HC2-HC4 (D17)" begin
    # One treated state per cohort, so every simple cell holds a single ATT(s,t)
    cohorts = Dict("s1" => 2003, "s3" => 2005)
    data = small_panel(cohorts = cohorts, effect = 1.0, noise = 0.3)
    simple = quiet_didint(data, cohorts; agg = "simple")
    sgt = quiet_didint(data, cohorts; agg = "sgt")
    for i in 1:nrow(simple)
        j = findfirst((sgt.gvar .== simple.gvar[i]) .& (sgt.t .== simple.time[i]))
        @test simple.se_att_gt[i] == sgt.se_att_sgt[j]
        @test simple.pval_att_gt[i] == sgt.pval_att_sgt[j]
    end
    # A single treated row has leverage 1, so HC2, HC3 and HC4 SEs of ATT(s,t) are NaN
    for hc in ["hc2", "hc3", "hc4"]
        @test all(isnan, quiet_didint(data, cohorts; agg = "sgt", hc = hc).se_att_sgt)
    end
end

@testset "Common adoption is the one-cohort case of the same estimator (D18)" begin
    data = small_panel(cohorts = COMMON_COHORTS, effect = 1.0, noise = 0.3)

    # Every agg runs and reports an aggregate ATT
    for agg in SECOND_STAGE_AGGS
        @test isfinite(quiet_didint(data, COMMON_COHORTS; agg = agg).agg_att[1])
    end

    # A single treatment time is broadcast to every treated state
    treated = sort(collect(keys(COMMON_COHORTS)))
    broadcast = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        DiDInt.didint("y", "state", "year", data; treated_states = treated, treatment_times = 2004,
                      agg = "simple", nperm = 19, seed = 1234, verbose = false)
    end
    @test isequal(broadcast, quiet_didint(data, COMMON_COHORTS; agg = "simple"))

    # The base period is g - 1: a shock in an earlier pre period leaves every estimate unchanged
    shocked = small_panel(cohorts = COMMON_COHORTS, effect = 1.0, noise = 0.3, shock = (state = "s1", year = 2002, size = 5.0))
    for agg in ["sgt", "cohort"]
        suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
        @test quiet_didint(shocked, COMMON_COHORTS; agg = agg)[!, suffix] ≈ quiet_didint(data, COMMON_COHORTS; agg = agg)[!, suffix]
    end

    # A one-cohort panel gives the same ATT(s,t)'s and ATT(g,t)'s as that cohort in a staggered panel
    stag_data = small_panel(cohorts = STAG_COHORTS, effect = 1.0, noise = 0.3)
    one_cohort = Dict("s1" => 2003, "s2" => 2003)
    one_data = filter(r -> r.state ∉ ["s3", "s4"], stag_data)
    for agg in ["sgt", "simple"]
        suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
        full = quiet_didint(stag_data, STAG_COHORTS; agg = agg)
        full = full[full.gvar .== Date(2003), :]
        single = quiet_didint(one_data, one_cohort; agg = agg)
        @test single[!, suffix] ≈ full[!, suffix]
        @test single[!, "se_" * suffix] ≈ full[!, "se_" * suffix]
    end

    # The edge case works with two states (one treated, one control)
    two = filter(r -> r.state in ["s1", "s5"], data)
    @test all(isfinite, quiet_didint(two, Dict("s1" => 2004); agg = "sgt", edgecase = true).se_att_sgt)
end

@testset "Jackknife needs at least 2 treated and 2 control states, overall and per sub-aggregate" begin
    # Cohort 2003 has two treated states, cohort 2005 one: only cohort 2003 gets a jackknife SE
    cohorts = Dict("s1" => 2003, "s2" => 2003, "s3" => 2005)
    data = small_panel(cohorts = cohorts, n_states = 6, effect = 1.0, noise = 0.3)
    results = quiet_didint(data, cohorts; agg = "cohort")
    @test !ismissing(results.jknifese_att_cohort[results.treatment_time .== Date(2003)][1])
    @test ismissing(results.jknifese_att_cohort[results.treatment_time .== Date(2005)][1])
    @test !ismissing(results.jknifese_agg_att[1])

    # With a single control state no jackknife SE or p-value is computed at all
    one_control = filter(r -> r.state in ["s1", "s2", "s3", "s4"], data)
    results = quiet_didint(one_control, cohorts; agg = "cohort")
    @test all(ismissing, results.jknifese_att_cohort)
    @test ismissing(results.jknifese_agg_att[1])
    @test ismissing(results.jknifepval_agg_att[1])
end
