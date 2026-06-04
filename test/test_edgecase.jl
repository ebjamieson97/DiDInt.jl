# 1) The edgecase procedure returns a finite SE where the normal procedure leaves missing
@testset "edgecase returns finite SE when diff regression is saturated but produces equivalent lambda (and consequently equal ATTs)" begin
    for ccc in ["int", "state", "time", "hom"]
        # staggered
        for agg in ["cohort", "state", "simple", "sgt", "time"]
            sub_att = if agg == "cohort"
                :att_cohort
            elseif agg == "state"
                :att_s
            elseif agg == "simple"
                :att_gt
            elseif agg == "sgt"
                :att_sgt
            elseif agg == "time"
                :att_t
            end
            sub_se = if agg == "cohort"
                :se_att_cohort
            elseif agg == "state"
                :se_att_s
            elseif agg == "simple"
                :se_att_gt
            elseif agg == "sgt"
                :se_att_sgt
            elseif agg == "time"
                :se_att_t
            end
            edge = DiDInt.didint("coll", "state", "year", TEST_DATA_EC_STAG;
                              treated_states = EC_STATES_STAG, treatment_times = EC_TIMES_STAG, edgecase = true, 
                              covariates = [:male, :asian, :black], ccc = ccc, agg = agg, nperm = 1)
            no_edge = DiDInt.didint("coll", "state", "year", TEST_DATA_EC_STAG;
                              treated_states = EC_STATES_STAG, treatment_times = EC_TIMES_STAG, edgecase = false, 
                              covariates = [:male, :asian, :black], ccc = ccc, agg = agg, nperm = 1)

            # First test that edge is never missing standard errors for the subaggregate ATTs
            @test all(x -> !ismissing(x) && isfinite(x), edge[:,sub_se])

            # Then test that some for the non edge case are missing
            @test any(x -> ismissing(x), no_edge[:,sub_se])

            # Test that the subaggregate ATTs are the same (and aggregate ATT while we are here)
            @test all(isapprox.(edge[:,sub_att], no_edge[:,sub_att]; atol = 1e-10))
            @test isapprox(edge.agg_att[1], no_edge.agg_att[1]; atol = 1e-10)

            # If there are any instances when they both report se, they should be equal
            edge = edge[.!ismissing.(no_edge[:,sub_se]),:]
            if (nrow(edge) > 0)
                no_edge = no_edge[.!ismissing.(no_edge[:,sub_se]),:]
                @test all(isapprox.(edge[:,sub_se], no_edge[:,sub_se]; atol = 1e-10))
            end
        end

        # common adoption
        for agg in ["state", "none"]
            if agg == "state"
                treated_states = ["71", "34"]
                treated_times = [1991, 1991]
                sub_se = :se_att_s
                sub_att = :att_s
            else
                treated_states = ["71"]
                treated_times = [1991]
            end
            edge = DiDInt.didint("coll", "state", "year", filter(r -> r.state ∈ [treated_states..., EC_CONTROL], TEST_DATA_FULL);
                              treated_states = treated_states, treatment_times = treated_times, edgecase = true,
                              covariates = [:male, :asian, :black], ccc = ccc, agg = agg, nperm = 1)
            no_edge = DiDInt.didint("coll", "state", "year", filter(r -> r.state ∈ [treated_states..., EC_CONTROL], TEST_DATA_FULL);
                              treated_states = treated_states, treatment_times = treated_times, edgecase = false,
                              covariates = [:male, :asian, :black], ccc = ccc, agg = agg, nperm = 1)
            
            @test isapprox(edge.agg_att[1], no_edge.agg_att[1]; atol = 1e-10)
            if agg == "state"
                # First test that edge is never missing standard errors for the subaggregate ATTs
                @test all(x -> !ismissing(x) && isfinite(x), edge[:,sub_se])

                # Then test that some for the non edge case are missing
                @test any(x -> ismissing(x), no_edge[:,sub_se])

                # Test that the subaggregate ATTs are the same (and aggregate ATT while we are here)
                @test all(isapprox.(edge[:,sub_att], no_edge[:,sub_att]; atol = 1e-10))
            end
        end
    end
end