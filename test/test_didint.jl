# Tests for didint()

# Try different aggregation methods for common adoption and different ccc
@testset "Agg check common - cohort (default)" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES,
                           seed = 1234,
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399)
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

# agg = "none" has been removed
@testset "agg = none raises an informative error" begin
    @test_throws r"agg = \"none\" has been removed" DiDInt.didint("coll", "state", "year", TEST_DATA,
                                                                  treatment_times = TREATED_TIMES,
                                                                  treated_states = TREATED_STATES,
                                                                  seed = 1234, agg = "none", nperm = 399)
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
@testset "truejack - common cohort" begin
    result = DiDInt.didint("coll", "state", "year", TEST_DATA_MISSING,
                           treatment_times = TREATED_TIMES[1],
                           treated_states = TREATED_STATES,
                           seed = 1234,
                           covariates = [:male, :asian, :black],
                           agg = "cohort", nperm = 399, truejack = true)
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



# RI p-values: p = (1 + #{|T*| >= |T|}) / (S + 1); placebo estimates that
# tie the observed estimate, exactly or up to rounding, count as hits
@testset "RI p-value counts ties, including rounding-level ties" begin
    observed = 0.3
    # 0.3 - 5.55e-17 and 0.1 + 0.2 are rounding-level ties; -0.3 is an exact tie in |T*|
    att_ri = [0.3 - 5.551115123125783e-17, 0.1 + 0.2, -0.3, 0.31, -0.5, 0.29, 0.0]
    @test DiDInt.ri_pvalue(att_ri, observed) == (1 + 5) / (7 + 1)
    @test DiDInt.ri_pvalue(att_ri, -observed) == (1 + 5) / (7 + 1)
    # Non-finite placebo estimates are dropped from S
    @test DiDInt.ri_pvalue([att_ri; NaN; Inf], observed) == (1 + 5) / (7 + 1)
    @test ismissing(DiDInt.ri_pvalue([NaN], observed))
    # Ordinary p-value for these placebos and T = 0.44 is (1 + 1) / 14
    placebos = [-0.21, 0.05, 0.33, -0.12, 0.47, -0.38, 0.09, 0.26, -0.02, 0.15, -0.29, 0.41, -0.07]
    @test DiDInt.ri_pvalue(placebos, 0.44) ≈ 2 / 14
end

# A two-category string covariate becomes one 0/1 dummy, the same as the numeric version
@testset "Two-category string covariate equals its numeric dummy" begin
    panel = CSV.read(joinpath(@__DIR__, "data", "didintpy_panel.csv"), DataFrame)
    panel.sex = ifelse.(panel.x2 .== 1.0, "female", "male")
    run(covariates) = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        DiDInt.didint("y", "state", "year", panel; treated_states = ["s01", "s02", "s03", "s04", "s05"],
                      treatment_times = [2004, 2004, 2006, 2006, 2007], covariates = covariates,
                      ccc = "int", agg = "sgt", nperm = 9, seed = 1234)
    end
    # "male" is the reference category (first seen), so the dummy is 1 for "female", i.e. x2
    @test run(["x1", "sex"]).att_sgt ≈ run(["x1", "x2"]).att_sgt atol = 1e-10
end
