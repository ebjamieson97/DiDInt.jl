# Tests of DiDIntResult, the object that didint() returns: its tables hold the same numbers as
# the wide table of results, and its display, accessors and table conversion work.
# The panel is test/data/didintpy_panel.csv (see test_policies.jl).

using Dates

result_didint(; kwargs...) = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
    DiDInt.didint("y", "state", "year", POLICY_PANEL; treated_states = POLICY_TREATED, treatment_times = POLICY_TIMES,
                  nperm = 29, seed = 1234, kwargs...)
end

# Two columns hold the same values, missing in the same places
same_column(a, b) = length(a) == length(b) && all((ismissing(x) && ismissing(y)) || isequal(x, y) for (x, y) in zip(a, b))

const RESULT_COLUMNS = ["att" => "", "se" => "se_", "pval" => "pval_", "jknife_se" => "jknifese_",
                        "jknife_pval" => "jknifepval_", "ri_pval" => "ri_pval_"]

@testset "didint() returns a DiDIntResult whose tables match the wide table" begin
    for agg in ["cohort", "simple", "state", "sgt", "time", "calendar"], weighting in ["both", "none"]
        result = result_didint(agg = agg, weighting = weighting, covariates = ["x1", "x2"])
        @test result isa DiDIntResult
        table = result.table
        suffix = DiDInt.SUB_AGG_SPECS[agg].suffix

        # Aggregate ATT
        @test result.agg.estimate == ["Aggregate ATT"]
        for (column, prefix) in RESULT_COLUMNS
            @test same_column(result.agg[!, column], table[1:1, prefix * "agg_att"])
            @test same_column(result.sub[!, column], table[!, prefix * suffix])
        end

        # Sub-aggregates: the key columns of the aggregation, then the estimates and the weight
        keys = [string(col) for (col, _, _) in DiDInt.SUB_AGG_SPECS[agg].keys]
        @test names(result.sub) == vcat(keys, first.(RESULT_COLUMNS), "weight")
        @test all(isequal(result.sub[!, key], table[!, key]) for key in keys)
        @test weighting == "both" ? result.sub.weight ≈ table.weights : all(ismissing, result.sub.weight)

        # The ATT of each treated state and period, and the specification
        @test names(result.att_st) == ["state", "gvar", "t", "r1", "periods_post_treat", "att", "se", "pval", "n_t", "n", "n_controls"]
        @test Set(result.att_st.state) == Set(POLICY_TREATED)
        if agg == "sgt"
            @test result.att_st.att ≈ result.sub.att
        end
        @test result.specs.agg == agg && result.specs.ccc == "int" && result.specs.weighting == weighting
        @test result.specs.covariates == ["x1", "x2"] && result.specs.nperm == 29 && result.specs.seed == 1234
        @test result.specs.period == "1 year" && result.specs.start_date == "2001" && result.specs.end_date == "2008"
    end
end

@testset "DiDIntResult for the policy aggregations" begin
    for agg in POLICY_AGGS
        # Without pooling: one row per policy
        result = result_didint(agg = agg, treated_policies = POLICY_LABELS)
        table = result.table
        first_rows = [findfirst(table.policy .== p) for p in ["A", "B"]]
        @test result.agg.estimate == ["Policy A", "Policy B"]
        @test !("policy_weight" in names(result.agg))
        for (column, prefix) in RESULT_COLUMNS
            @test same_column(result.agg[!, column], table[first_rows, prefix * "att_p"])
        end
        @test "policy" in names(result.att_st) && first(names(result.sub)) == "policy"

        # With pooling: the pooled ATT first, and each policy's weight in it
        pooled = result_didint(agg = agg, treated_policies = POLICY_LABELS, policy_aggregation = true, policy_weighting = "treated_obs")
        @test pooled.agg.estimate == ["Pooled ATT", "Policy A", "Policy B"]
        @test pooled.agg.att[1] ≈ pooled.table.att_pool[1]
        @test pooled.agg.ri_pval[1] ≈ pooled.table.ri_pval_att_pool[1]
        @test ismissing(pooled.agg.policy_weight[1]) && sum(pooled.agg.policy_weight[2:3]) ≈ 1
        @test pooled.agg.att[1] ≈ sum(pooled.agg.policy_weight[2:3] .* pooled.agg.att[2:3])
    end
end

@testset "The Stata and R wrappers receive the wide table" begin
    result = result_didint(agg = "cohort")
    for wrapper in ["stata", "r"]
        table = result_didint(agg = "cohort", wrapper = wrapper)
        @test table isa DataFrame
        @test names(table) == names(result.table)
        @test same_column(table.att_cohort, result.table.att_cohort) && same_column(table.ri_pval_agg_att, result.table.ri_pval_agg_att)
    end
    # Dates are passed to R as strings
    @test eltype(result_didint(agg = "cohort", wrapper = "r").treatment_time) <: AbstractString
    @test eltype(result.table.treatment_time) == Date
    @test result_didint(agg = "cohort", wrapper = "julia") isa DiDIntResult
end

@testset "coef, stderror, coefnames and conversion to a table" begin
    result = result_didint(agg = "cohort")
    @test coef(result) == result.agg.att && coef(result; level = :sub) == result.sub.att
    @test stderror(result) == result.agg.se
    @test same_column(stderror(result; level = :sub, type = :jackknife), result.sub.jknife_se)
    @test coefnames(result) == ["Aggregate ATT"]
    @test coefnames(result; level = :sub) == ["2004-01-01", "2006-01-01", "2007-01-01"]
    @test coefnames(result_didint(agg = "time"); level = :sub)[1:2] == ["0", "1"]
    @test coefnames(result_didint(agg = "sgt"); level = :sub)[1] == "s01; 2004-01-01; 2004-01-01"
    @test_throws r"level" coef(result; level = :all)
    @test_throws r"type" stderror(result; type = :bootstrap)

    # DataFrame(result) is a copy of the wide table, and the result is a Tables.jl table
    table = DataFrame(result)
    @test table isa DataFrame && isequal(table, result.table) && table !== result.table
    @test DiDInt.Tables.istable(result)
    path = tempname() * ".csv"
    CSV.write(path, result)
    @test names(CSV.read(path, DataFrame)) == names(result.table)

    # Without weights for the sub-aggregates (weighting = "none" or "diff") the weights are
    # missing, and the result can still be written to a file; the wrappers' table is unchanged
    for weighting in ["none", "diff"], agg in ["cohort", "policy_cohort"]
        unweighted = result_didint(agg = agg, weighting = weighting)
        @test all(ismissing, unweighted.table.weights) && all(ismissing, unweighted.sub.weight)
        CSV.write(path, unweighted)
        @test nrow(CSV.read(path, DataFrame)) == nrow(unweighted.table)
        @test all(isnothing, result_didint(agg = agg, weighting = weighting, wrapper = "stata").weights)
    end
end

@testset "DiDIntResult display" begin
    result = result_didint(agg = "cohort", covariates = ["x1"], ccc = "hom", notyet = true)
    text = sprint(show, MIME("text/plain"), result)
    for expected in ["DID-INT results", "Outcome:          y", "Homogeneous DID-INT (ccc = \"hom\")", "Covariates:       x1",
                     "\"cohort\" (by treatment time)", "never-treated and not-yet-treated states", "1 year, from 2001 to 2008",
                     "HC1", "29 (seed 1234)", "Aggregate ATT", "Sub-aggregate ATTs by treatment time (3)", "jknife_se", "ri_pval"]
        @test occursin(expected, text)
    end
    @test occursin("No covariates", sprint(show, MIME("text/plain"), result_didint(agg = "cohort")))
    @test startswith(sprint(show, result), "DiDIntResult(agg = \"cohort\", ccc = \"hom\", 3 sub-aggregate ATTs, aggregate ATT = ")

    # Long tables of sub-aggregates are cut off, with a pointer to the full table
    merit_result = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        DiDInt.didint("coll", "state", "year", TEST_DATA; treated_states = TREATED_STATES, treatment_times = TREATED_TIMES,
                      agg = "sgt", nperm = 9, seed = 1234)
    end
    n_sub = nrow(merit_result.sub)
    @test n_sub > 20
    @test occursin("... and $(n_sub - 20) more; see the `sub` field", sprint(show, MIME("text/plain"), merit_result))
    @test !occursin("more; see the `sub` field", sprint(show, MIME("text/plain"), result_didint(agg = "sgt")))

    # Policy aggregations
    pooled = sprint(show, MIME("text/plain"), result_didint(agg = "policy", treated_policies = POLICY_LABELS, policy_aggregation = true))
    @test occursin("Pooled ATT and ATT by policy", pooled) && !occursin("Sub-aggregate ATTs", pooled)
    by_policy = sprint(show, MIME("text/plain"), result_didint(agg = "policy_time", treated_policies = POLICY_LABELS))
    @test occursin("ATT by policy", by_policy) && occursin("Sub-aggregate ATTs by policy and periods since treatment", by_policy)
end
