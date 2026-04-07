# Tests for didint()

# Try different aggregation methods for common adoption and different ccc
@testset "Agg check common - none (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check common - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234,
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399)
    @test !isnothing(result)
end

# Try different aggregation methods for staggered adoption
@testset "Agg check staggered - cohort (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "state",
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - simple" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "hom",
                           covariates = [:male, :asian, :black],
                           agg = "simple", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - sgt" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234,  ccc = "add",
                           covariates = [:male, :asian, :black],
                           agg = "sgt", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "time",
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - time" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "time", nperm = 399)
    @test !isnothing(result)
end

@testset "Agg check staggered - none" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399)
    @test !isnothing(result)
end

# Check that gvar works for staggered adoption
@testset "Gvar check - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_GVAR,
                           gvar = :gvar, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black], nperm = 399)
    @test !isnothing(result)
end

# Test that string time column works, also check start_date, end_date, and freq arguments
@testset "String time column, start_date, end_date, freq - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_STR,
                           treatment_times = string.(TREATED_TIMES),
                           date_format = "yyyy", start_date = "1989", end_date = "2000", freq = "yearly",
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black], nperm = 399)
    @test !isnothing(result)
end

# Test usage of notyet and different weighting option
@testset "notyet and weighting = none - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, notyet = true, weighting = "none",
                           covariates = [:male, :asian, :black], nperm = 399)
    @test !isnothing(result)
end

# Test that common or staggered work with missing data
@testset "Missing data - staggered" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399)
    @test !isnothing(result)
end

@testset "Missing data - common" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399)
    @test !isnothing(result)
end

# Test that the truejack option works 
@testset "truejack - common none" begin 
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399, truejack = true)
    @test !isnothing(result)
end 

@testset "truejack - common state" begin 
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399, truejack = true)
    @test !isnothing(result)
end 

@testset "truejack check staggered - cohort (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "state",
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399, truejack = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - simple" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "hom",
                           covariates = [:male, :asian, :black],
                           agg = "simple", nperm = 399, truejack = true,
                           notyet = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - sgt" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234,  ccc = "add",
                           covariates = [:male, :asian, :black],
                           agg = "sgt", nperm = 399, truejack = true,
                           notyet = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - state" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, ccc = "time",
                           covariates = [:male, :asian, :black],
                           agg = "state", nperm = 399, truejack = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - time" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "time", nperm = 399, truejack = true)
    @test !isnothing(result)
end


@testset "truejack and notyet check staggered - time" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "time", nperm = 399, truejack = true, notyet = true)
    @test !isnothing(result)
end

@testset "truejack check staggered - none" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES,
                           treated_states = TREATED_STATES, 
                           seed = 1234, 
                           covariates = [:male, :asian, :black],
                           agg = "none", nperm = 399, truejack = true)
    @test !isnothing(result)
end

@testset "CCC collinearity recovery" begin

    @testset "ccc=state recover=true" begin
        result = DiDInt.didint("coll", "state", "year", TEST_DATA_STATE_COLLINEAR,
                               treatment_times = TREATED_TIMES,
                               treated_states = TREATED_STATES,
                               seed = 1234, ccc = "state", recover = true,
                               covariates = [:male],
                               agg = "cohort", nperm = 399)
        @test !isnothing(result)
        @test !ismissing(result.agg_att[1])
    end

    @testset "ccc=state recover=false" begin
        result = DiDInt.didint("coll", "state", "year", TEST_DATA_STATE_COLLINEAR,
                               treatment_times = TREATED_TIMES,
                               treated_states = TREATED_STATES,
                               seed = 1234, ccc = "state", recover = false,
                               covariates = [:male],
                               agg = "cohort", nperm = 399)
        @test !isnothing(result)
        @test !ismissing(result.agg_att[1])
    end

    @testset "ccc=state recover true vs false differ" begin
        result_true = DiDInt.didint("coll", "state", "year", TEST_DATA_STATE_COLLINEAR,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "state", recover = true,
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399)
        result_false = DiDInt.didint("coll", "state", "year", TEST_DATA_STATE_COLLINEAR,
                                     treatment_times = TREATED_TIMES,
                                     treated_states = TREATED_STATES,
                                     seed = 1234, ccc = "state", recover = false,
                                     covariates = [:male, :asian, :black],
                                     agg = "cohort", nperm = 399)

        @test isapprox(result_true.agg_att[1], result_false.agg_att[1], atol=1e-10)
    end

    @testset "ccc=time recover=true" begin
        result = DiDInt.didint("coll", "state", "year", TEST_DATA_TIME_COLLINEAR,
                               treatment_times = TREATED_TIMES,
                               treated_states = TREATED_STATES,
                               seed = 1234, ccc = "time", recover = true,
                               covariates = [:male],
                               agg = "cohort", nperm = 399)
        @test !isnothing(result)
        @test !ismissing(result.agg_att[1])
    end

    @testset "ccc=time recover=false" begin
        result = DiDInt.didint("coll", "state", "year", TEST_DATA_TIME_COLLINEAR,
                               treatment_times = TREATED_TIMES,
                               treated_states = TREATED_STATES,
                               seed = 1234, ccc = "time", recover = false,
                               covariates = [:male],
                               agg = "cohort", nperm = 399)
        @test !isnothing(result)
    end

    @testset "ccc=time recover true vs false differ" begin
        result_true = DiDInt.didint("coll", "state", "year", TEST_DATA_TIME_COLLINEAR,
                                    treatment_times = TREATED_TIMES_T_COLLINEAR,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "time", recover = true,
                                    covariates = [:male],
                                    agg = "cohort", nperm = 399)
        result_false = DiDInt.didint("coll", "state", "year", TEST_DATA_TIME_COLLINEAR,
                                     treatment_times = TREATED_TIMES_T_COLLINEAR,
                                     treated_states = TREATED_STATES,
                                     seed = 1234, ccc = "time", recover = false,
                                     covariates = [:male],
                                     agg = "cohort", nperm = 399)
        @test isapprox(result_true.agg_att[1], result_false.agg_att[1], atol=1e-10)
    end

    @testset "ccc=hom recover=true" begin
        result = DiDInt.didint("coll", "state", "year", TEST_DATA_HOM_COLLINEAR,
                               treatment_times = TREATED_TIMES,
                               treated_states = TREATED_STATES,
                               seed = 1234, ccc = "hom", recover = true,
                               covariates = [:male],
                               agg = "cohort", nperm = 399)
        @test !isnothing(result)
    end

    @testset "ccc=hom recover=false" begin
        result = DiDInt.didint("coll", "state", "year", TEST_DATA_HOM_COLLINEAR,
                               treatment_times = TREATED_TIMES,
                               treated_states = TREATED_STATES,
                               seed = 1234, ccc = "hom", recover = false,
                               covariates = [:male],
                               agg = "cohort", nperm = 399)
        @test !isnothing(result)
    end

    @testset "ccc=int recover true vs false differ" begin
        result_true = DiDInt.didint("coll", "state", "year", TEST_DATA_FULL,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "int", recover = true,
                                    fem = true, iterative = false,
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399)
        result_false = DiDInt.didint("coll", "state", "year", TEST_DATA_FULL,
                                     treatment_times = TREATED_TIMES,
                                     treated_states = TREATED_STATES,
                                     seed = 1234, ccc = "int", recover = false,
                                     fem = true, iterative = false,
                                     covariates = [:male, :asian, :black],
                                     agg = "cohort", nperm = 399)
        @test (result_true.agg_att[1] != result_false.agg_att[1])
    end

end

@testset "make sure that edge case se works" begin
    
    @testset "common adoption" begin
        result = DiDInt.didint("coll", "state", "year", TEST_DATA_FULL[(TEST_DATA_FULL.state .== "71") .| (TEST_DATA_FULL.state .== "73"), :],
                            treatment_times = 1991,
                            treated_states = "71",
                            seed = 1234, ccc = "state", recover = true,
                            covariates = nothing,
                            fem = true, iterative = false,
                            agg = "cohort", nperm = 399)
        @test (!ismissing(result.se_agg_att))
    end

    @testset "staggered adoption" begin
        result_sgt = DiDInt.didint("coll", "state", "year", filter(row -> row.state ∈ vcat(CONTROL_STATES[1], TREATED_STATES), TEST_DATA_FULL),
                                     treatment_times = TREATED_TIMES,
                                     treated_states = TREATED_STATES,
                                     seed = 1234, ccc = "int", recover = false,
                                     fem = true, iterative = false,
                                     covariates = [:male, :asian, :black],
                                     agg = "sgt", nperm = 399)

        result_simple = DiDInt.didint("coll", "state", "year", filter(row -> row.state ∈ vcat(CONTROL_STATES[1], TREATED_STATES), TEST_DATA_FULL),
                                     treatment_times = TREATED_TIMES,
                                     treated_states = TREATED_STATES,
                                     seed = 1234, ccc = "add", recover = false,
                                     fem = true, iterative = false,
                                     covariates = [:male, :asian, :black],
                                     agg = "simple", nperm = 399)
        @test (all(x -> !ismissing(x), result_sgt.se_att_sgt) && all(x -> !ismissing(x), result_simple.se_att_gt))
    end

end

@testset "Iterative vs FEM comparison" begin
    
    @testset "ccc=int: iterative matches FEM" begin
        result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "int",
                                   covariates = [:male, :asian, :black],
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "int",
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        # Test that ATT estimates are approximately equal
        @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-8)
        
    end
    
    @testset "ccc=time: iterative matches FEM" begin
        result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "time",
                                   covariates = [:male, :asian, :black],
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "time",
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-8)
    end
    
    @testset "ccc=state: iterative matches FEM" begin
        result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "state",
                                   covariates = [:male, :asian, :black],
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "state",
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-8)
    end
    
    @testset "ccc=add: iterative matches FEM" begin
        result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "add",
                                   covariates = [:male, :asian, :black],
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "add",
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-8)
    end
    
    @testset "ccc=hom: iterative matches FEM" begin
        result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "hom",
                                   covariates = [:male],
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "hom",
                                    covariates = [:male],
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-8)
    end
    
    @testset "No covariates: iterative matches FEM" begin
        result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "int",
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "int",
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-8)
    end
    
end

@testset "Missing data - iterative vs FEM" begin
    @testset "Staggered with missing: iterative matches FEM" begin
        result_fem_rec = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "int",
                                   covariates = [:male, :asian, :black],
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "int",
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        @test !isnothing(result_fem_rec)
        @test !isnothing(result_iter)
        @test isapprox(result_fem_rec.agg_att[1], result_iter.agg_att[1], atol=1e-8)
    end
    
    @testset "Different CCC with missing data" begin
        for ccc_type in ["int", "time", "state", "add", "hom"]
            result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                                       treatment_times = TREATED_TIMES,
                                       treated_states = TREATED_STATES,
                                       seed = 1234, ccc = ccc_type,
                                       covariates = [:male, :asian, :black],
                                       agg = "cohort", nperm = 399,
                                       iterative = false, fem = true)
            
            result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                                        treatment_times = TREATED_TIMES,
                                        treated_states = TREATED_STATES,
                                        seed = 1234, ccc = ccc_type,
                                        covariates = [:male, :asian, :black],
                                        agg = "cohort", nperm = 399,
                                        iterative = true, fem = false)
            
            @test !isnothing(result_fem)
            @test !isnothing(result_iter)
            @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-7)
        end
    end
end

@testset "Collinearity handling - iterative vs FEM" begin
    @testset "ccc=state collinear: iterative matches FEM" begin
        result_fem = DiDInt.didint("coll", "state", "year", TEST_DATA_STATE_COLLINEAR,
                                   treatment_times = TREATED_TIMES,
                                   treated_states = TREATED_STATES,
                                   seed = 1234, ccc = "state", recover = false,
                                   covariates = [:male, :asian, :black],
                                   agg = "cohort", nperm = 399,
                                   iterative = false, fem = true)
        
        result_iter = DiDInt.didint("coll", "state", "year", TEST_DATA_STATE_COLLINEAR,
                                    treatment_times = TREATED_TIMES,
                                    treated_states = TREATED_STATES,
                                    seed = 1234, ccc = "state", recover = true,
                                    covariates = [:male, :asian, :black],
                                    agg = "cohort", nperm = 399,
                                    iterative = true, fem = false)
        
        @test isapprox(result_fem.agg_att[1], result_iter.agg_att[1], atol=1e-8)
    end
end