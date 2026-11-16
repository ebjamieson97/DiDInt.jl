# The edgecase procedure: an ATT(s,t) whose regression has exactly two rows
# (the treated state and a single control) gets its SE from the lambdas' covariance matrix
# instead of leaving it missing. Estimates must not depend on edgecase.
@testset "edgecase returns finite SEs for two-row ATT(s,t) regressions, with identical estimates" begin
    for ccc in ["int", "state", "time", "hom"]
        # Staggered: two treated states in different cohorts and one control, so every ATT(s,t)
        # regression has two rows. Sub-aggregates with a single ATT(s,t) report its SE;
        # those with several get their SE from the spread of their ATT(s,t)'s either way.
        for agg in ["cohort", "state", "simple", "sgt", "time", "calendar"]
            suffix = DiDInt.SUB_AGG_SPECS[agg].suffix
            sub_att = Symbol(suffix)
            sub_se = Symbol("se_" * suffix)
            edge = DiDInt.didint("coll", "state", "year", TEST_DATA_EC_STAG;
                                 treated_states = EC_STATES_STAG, treatment_times = EC_TIMES_STAG, edgecase = true,
                                 covariates = [:male, :asian, :black], ccc = ccc, agg = agg, nperm = 1)
            no_edge = DiDInt.didint("coll", "state", "year", TEST_DATA_EC_STAG;
                                    treated_states = EC_STATES_STAG, treatment_times = EC_TIMES_STAG, edgecase = false,
                                    covariates = [:male, :asian, :black], ccc = ccc, agg = agg, nperm = 1)

            # The edgecase run never leaves a sub-aggregate SE missing; the normal run does
            @test all(x -> !ismissing(x) && isfinite(x), edge[:, sub_se])
            @test any(ismissing, no_edge[:, sub_se])

            # Sub-aggregate and aggregate ATTs are the same
            @test all(isapprox.(edge[:, sub_att], no_edge[:, sub_att]; atol = 1e-10))
            @test isapprox(edge.agg_att[1], no_edge.agg_att[1]; atol = 1e-10)

            # Where both report an SE, they agree
            both = .!ismissing.(no_edge[:, sub_se])
            @test all(isapprox.(edge[both, sub_se], no_edge[both, sub_se]; atol = 1e-10))
        end

        # Common adoption (one cohort) with one control: every ATT(s,t) is an edge case
        for treated_states in [["71"], ["71", "34"]]
            data = filter(r -> r.state ∈ [treated_states..., EC_CONTROL], TEST_DATA_FULL)
            edge = DiDInt.didint("coll", "state", "year", data;
                                 treated_states = treated_states, treatment_times = 1991, edgecase = true,
                                 covariates = [:male, :asian, :black], ccc = ccc, agg = "sgt", nperm = 1)
            no_edge = DiDInt.didint("coll", "state", "year", data;
                                    treated_states = treated_states, treatment_times = 1991, edgecase = false,
                                    covariates = [:male, :asian, :black], ccc = ccc, agg = "sgt", nperm = 1)
            @test all(x -> !ismissing(x) && isfinite(x), edge.se_att_sgt)
            @test all(ismissing, no_edge.se_att_sgt)
            @test all(isapprox.(edge.att_sgt, no_edge.att_sgt; atol = 1e-10))
            @test isapprox(edge.agg_att[1], no_edge.agg_att[1]; atol = 1e-10)
        end
    end
end
