# Tests of how collinear covariates are handled in the first stage.
# Under ccc = "hom", "state", "time" and "int", a covariate that is a linear combination of
# covariates listed before it is dropped wherever the slopes are estimated, with a warning,
# and the regular and the edge-case procedure use the same covariates. Under ccc = "add" no
# covariate is dropped, and a warning is given if the design has more rank deficiency than
# the split of the slopes between states and periods.
# The panel is INTEGRATION_PANEL (see test_integration.jl): 12 states observed 2001-2008.

using Dates
using LinearAlgebra
using Statistics

# The warnings that a call raises, as strings
function warnings_of(f)
    logger = Test.TestLogger(min_level = Base.CoreLogging.Warn)
    Base.CoreLogging.with_logger(f, logger)
    return String[string(record.message) for record in logger.logs]
end
dropped_warnings(f) = filter(w -> occursin("were dropped", w), warnings_of(f))
rank_warnings(f) = filter(w -> occursin("has rank", w), warnings_of(f))

# 2 states x 2 periods with 20 observations per cell, state A treated in 2002, and
# x2 = 2 x1 + c_{s,t} with c = (0, 1, 3, -2): after demeaning within cells x2 is exactly
# 2 x1, but the two covariates differ by constants that are not additive in state and period.
function two_by_two_panel(; scale = 1)
    shift = Dict(("A", 2001) => 0.0, ("A", 2002) => 1.0, ("B", 2001) => 3.0, ("B", 2002) => -2.0)
    panel = DataFrame(state = String[], year = Int[], x1 = Float64[], x2 = Float64[], y = Float64[])
    i = 0
    for state in ["A", "B"], year in [2001, 2002], _ in 1:20
        i += 1
        x1 = 2 * sin(12.9898 * i) + (state == "A") + 0.5 * (year == 2002)
        y = 1 + (state == "A") + 0.4 * (year == 2002) + 0.7 * x1 + 0.5 * (state == "A" && year == 2002) + 0.3 * sin(78.233 * i)
        push!(panel, (state, year, x1, scale * (2 * x1 + shift[(state, year)]), y))
    end
    return panel
end

two_by_two_didint(panel; kwargs...) = quiet_call() do
    didint_table("y", "state", "year", panel; treated_states = "A", treatment_times = 2002, agg = "sgt",
                 nperm = 9, seed = 1234, kwargs...)
end

# The lambdas and covariate slopes of the first stage, called as didint() calls it
function first_stage_of(panel, ccc, covariates; edgecase = false, treated_states = ["A"], treatment_times = [2002])
    data = DataFrame(state_71X9yTx = String.(panel.state), time_71X9yTx = string.(Date.(panel.year)),
                     time_dmG5fpM = Date.(panel.year), outcome_71X9yTx = Float64.(panel.y))
    for cov in covariates
        data[!, cov] = Float64.(panel[!, cov])
    end
    dates = Date.(sort(unique(panel.year)))
    lambda_df, _, fit = DiDInt.first_stage(data, ccc, covariates, "hc1", edgecase;
                                           time_to_index = Dict(d => i for (i, d) in enumerate(dates)),
                                           treatment_times = Date.(treatment_times), match_to_these_dates = dates,
                                           treated_states = treated_states, use_pre_controls = false)
    return sort(lambda_df, [:state, :time]), fit
end

@testset "2 x 2 example: a covariate collinear with another up to cell-specific constants" begin
    for ccc in ["hom", "state", "time", "int"]
        base = two_by_two_didint(two_by_two_panel(); ccc = ccc, covariates = ["x1", "x2"])
        @test isfinite(base.att_sgt[1])

        # x2 is dropped: the estimates are those with x1 alone, whatever the units of x2
        @test base.att_sgt ≈ two_by_two_didint(two_by_two_panel(); ccc = ccc, covariates = ["x1"]).att_sgt rtol = 1e-9 atol = 1e-10
        for scale in [10, 100]
            rescaled = two_by_two_didint(two_by_two_panel(scale = scale); ccc = ccc, covariates = ["x1", "x2"])
            @test rescaled.att_sgt ≈ base.att_sgt rtol = 1e-9 atol = 1e-10
        end

        # The order of the covariates decides which one is dropped, and the estimate differs
        reversed = two_by_two_didint(two_by_two_panel(); ccc = ccc, covariates = ["x2", "x1"])
        @test reversed.att_sgt ≈ two_by_two_didint(two_by_two_panel(); ccc = ccc, covariates = ["x2"]).att_sgt rtol = 1e-9 atol = 1e-10
        @test abs(reversed.att_sgt[1] - base.att_sgt[1]) > 1e-3

        # One warning, which names the dropped covariate
        warnings = dropped_warnings() do
            didint("y", "state", "year", two_by_two_panel(); treated_states = "A", treatment_times = 2002, ccc = ccc,
                   covariates = ["x1", "x2"], nperm = 9, seed = 1234)
        end
        @test length(warnings) == 1
        @test occursin("x2", warnings[1]) && !occursin("x1", warnings[1])
        @test occursin("order of 'covariates'", warnings[1])

        # The edge-case procedure uses the same covariates: same lambdas and ATT, and its
        # standard error is computed, equal to the one with x1 alone
        edge = two_by_two_didint(two_by_two_panel(); ccc = ccc, covariates = ["x1", "x2"], edgecase = true)
        @test edge.att_sgt ≈ base.att_sgt rtol = 1e-9 atol = 1e-10
        @test !ismissing(edge.se_att_sgt[1]) && isfinite(edge.se_att_sgt[1]) && edge.se_att_sgt[1] > 0
        @test edge.se_att_sgt ≈ two_by_two_didint(two_by_two_panel(); ccc = ccc, covariates = ["x1"], edgecase = true).se_att_sgt rtol = 1e-8
        for scale in [1, 10, 100]
            regular, regular_fit = quiet_call(() -> first_stage_of(two_by_two_panel(scale = scale), ccc, ["x1", "x2"]))
            edgecase, edgecase_fit = quiet_call(() -> first_stage_of(two_by_two_panel(scale = scale), ccc, ["x1", "x2"]; edgecase = true))
            @test edgecase.lambda ≈ regular.lambda rtol = 1e-9 atol = 1e-10
            for cell in keys(regular_fit.slopes)
                @test edgecase_fit.slopes[cell] ≈ regular_fit.slopes[cell] rtol = 1e-8 atol = 1e-10
                @test regular_fit.slopes[cell][2] == 0
            end
        end
    end
end

@testset "Without collinear covariates the estimates are those of version 0.10.0" begin
    # Aggregate ATT, cohort ATTs and the standard error of the aggregate ATT from version 0.10.0
    reference = Dict(
        "hom" => (0.32255093073673835, [0.4064516769203454, 0.3608603795575715, -0.1861596066093344], 0.08964016312937668),
        "state" => (0.32434005687427997, [0.4118183131818078, 0.35860275841851486, -0.18725819074089992], 0.09095480699246158),
        "time" => (0.34692742156698286, [0.40656155496421037, 0.3534909673855189, 0.05376169718720666], 0.055297047869979285),
        "add" => (0.3455123412203705, [0.40744108833446363, 0.34248610921914313, 0.07364763960166432], 0.05414024433071098),
        "int" => (0.30545835353863815, [0.3437628706680567, 0.3432756093706496, 0.005918601091867316], 0.05092743240024002))
    for (ccc, (agg_att, att_cohort, se_agg_att)) in reference
        call() = didint_table("y", "state", "year", INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                              ccc = ccc, agg = "cohort", covariates = ["x1", "x2"], nperm = 9, seed = 1234)
        result = quiet_call(call)
        @test result.agg_att[1] ≈ agg_att rtol = 1e-9
        @test result.att_cohort ≈ att_cohort rtol = 1e-9
        @test result.se_agg_att[1] ≈ se_agg_att rtol = 1e-9
        # No covariate is dropped, and "add" finds no rank deficiency beyond the expected one
        warnings = warnings_of(call)
        @test !any(w -> occursin("were dropped", w) || occursin("has rank", w), warnings)

        # A covariate may have any name, also one that looks like an internal column
        renamed = rename(INTEGRATION_PANEL, :x2 => :x1_demeaned)
        result = quiet_call() do
            didint_table("y", "state", "year", renamed; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                         ccc = ccc, agg = "cohort", covariates = ["x1", "x1_demeaned"], nperm = 9, seed = 1234)
        end
        @test result.att_cohort ≈ att_cohort rtol = 1e-9
    end
end

@testset "ccc = int: the covariate listed last is the one dropped in a cell" begin
    in_cell = (INTEGRATION_PANEL.state .== "s06") .& (INTEGRATION_PANEL.year .== 2005)
    cell = ("s06", Date(2005))
    run(panel, covariates; kwargs...) = input_didint(panel; ccc = "int", agg = "cohort", covariates = covariates, kwargs...)

    # x2 is an exact linear function of x1 in one cell, with a much smaller scale than x1.
    # Version 0.10.0 chose the covariate to drop by its norm, which also dropped x2 here, so
    # the estimates are unchanged (aggregate and cohort ATTs from version 0.10.0).
    small = copy(INTEGRATION_PANEL)
    small.x2 = Float64.(small.x2)
    small.x2[in_cell] = 2 .+ 0.01 .* small.x1[in_cell]
    result = run(small, ["x1", "x2"])
    @test result.agg_att[1] ≈ 0.3055077428073193 rtol = 1e-9
    @test result.att_cohort ≈ [0.343887425165474, 0.3432365661247329, 0.005918601091867316] rtol = 1e-9

    # With a larger scale than x1 it is still x2 that is dropped in that cell (version 0.10.0
    # dropped x1), and with the covariates listed the other way round it is x1. Elsewhere
    # both covariates are kept.
    large = copy(INTEGRATION_PANEL)
    large.x2 = Float64.(large.x2)
    large.x2[in_cell] = 2 .+ 3 .* large.x1[in_cell]
    _, fit = quiet_call(() -> first_stage_of(large, "int", ["x1", "x2"]))
    _, unchanged_fit = quiet_call(() -> first_stage_of(INTEGRATION_PANEL, "int", ["x1", "x2"]))
    @test fit.slopes[cell][1] != 0 && fit.slopes[cell][2] == 0
    @test unchanged_fit.slopes[cell][2] != 0
    @test all(slopes == unchanged_fit.slopes[key] for (key, slopes) in fit.slopes if key != cell)
    _, fit_reversed = quiet_call(() -> first_stage_of(large, "int", ["x2", "x1"]))
    @test fit_reversed.slopes[cell][1] != 0 && fit_reversed.slopes[cell][2] == 0
    warnings = dropped_warnings(() -> first_stage_of(large, "int", ["x1", "x2"]))
    @test length(warnings) == 1
    @test occursin("x2 (state s06 in period 2005-01-01)", warnings[1])
    # Rescaling x2 changes nothing
    rescaled = copy(large)
    rescaled.x2 = 1000 .* rescaled.x2
    @test run(rescaled, ["x1", "x2"]).att_cohort ≈ run(large, ["x1", "x2"]).att_cohort rtol = 1e-9 atol = 1e-10

    # A factor whose reference category "u" does not occur in one cell: the dummies of "v" and
    # "w" add up to one there, and the dummy listed last ("w") is dropped in that cell only
    factor = copy(INTEGRATION_PANEL)
    factor.cat = [["u", "v", "w"][1 + mod(i, 3)] for i in 1:nrow(factor)]
    factor.cat[in_cell] = [["v", "w", "w"][1 + mod(i, 3)] for i in 1:count(in_cell)]
    factor.cat_v = Float64.(factor.cat .== "v")
    factor.cat_w = Float64.(factor.cat .== "w")
    _, fit = quiet_call(() -> first_stage_of(factor, "int", ["cat_v", "cat_w"]))
    @test fit.slopes[cell][1] != 0 && fit.slopes[cell][2] == 0
    @test all(all(!iszero, slopes) for (key, slopes) in fit.slopes if key != cell)
    warnings = dropped_warnings() do
        didint("y", "state", "year", factor; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES, ccc = "int",
               covariates = ["cat"], ref = Dict("cat" => "u"), nperm = 9, seed = 1234)
    end
    @test length(warnings) == 1
    @test occursin("cat_w (state s06 in period 2005-01-01)", warnings[1]) && !occursin("cat_v", warnings[1])

    # A cell with fewer observations than covariates plus one keeps at most (observations - 1)
    # covariates, the ones listed first
    tiny = INTEGRATION_PANEL[.!in_cell .| (cumsum(in_cell) .<= 2), :]
    tiny.x2 = Float64.(tiny.x2)
    tiny.x3 = tiny.x1 .^ 2
    _, fit = quiet_call(() -> first_stage_of(tiny, "int", ["x1", "x2", "x3"]))
    @test fit.slopes[cell][1] != 0 && fit.slopes[cell][2:3] == [0, 0]
end

@testset "Covariates without within-cell variation are dropped without a warning" begin
    # A covariate measured at the state-period level, and one that only varies within the
    # cells of some states: neither is collinear with another covariate
    panel = copy(INTEGRATION_PANEL)
    cells = Dict(cell => i for (i, cell) in enumerate(unique(zip(panel.state, panel.year))))
    panel.unemp = [round(3 + 6 * abs(sin(12.9898 * cells[cell])), digits = 1) for cell in zip(panel.state, panel.year)]
    for ccc in ["hom", "state", "time", "int"]
        call() = didint("y", "state", "year", panel; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES, ccc = ccc,
                        covariates = ["x1", "unemp", "x2"], nperm = 9, seed = 1234)
        @test isempty(dropped_warnings(call))
        # unemp has slope 0 everywhere, and the other slopes are those without it
        _, fit = quiet_call(() -> first_stage_of(panel, ccc, ["x1", "unemp", "x2"]))
        _, without_fit = quiet_call(() -> first_stage_of(panel, ccc, ["x1", "x2"]))
        @test all(slopes[2] == 0 && slopes[[1, 3]] ≈ without_fit.slopes[key] for (key, slopes) in fit.slopes)
    end

    # The pruning rule itself: a column of zeros is neither kept nor reported as collinear
    X = hcat([1.0, -1.0, 2.0, -2.0], zeros(4), [2.0, -2.0, 4.0, -4.0], [1.0, 1.0, -1.0, -1.0])
    @test DiDInt.prune_covariate_columns(X) == ([1, 4], [3])
    # ... and it does not depend on the units of the columns
    @test DiDInt.prune_covariate_columns(X .* [1e-6 1.0 1e6 1e3]) == ([1, 4], [3])
    # A column that is nearly, but not within the tolerance, a combination of others is kept
    @test DiDInt.prune_covariate_columns(hcat(X[:, 1], X[:, 1] .+ 1e-5 .* X[:, 4])) == ([1, 2], Int[])
    @test DiDInt.prune_covariate_columns(hcat(X[:, 1], X[:, 1] .+ 1e-9 .* X[:, 4])) == ([1], [2])
end

@testset "ccc = add: rank deficiency beyond the split between states and periods" begin
    add_call(panel, covariates) = () -> didint("y", "state", "year", panel; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                                               ccc = "add", covariates = covariates, nperm = 9, seed = 1234)

    # Well-behaved covariates: no warning (the estimates are checked against version 0.10.0 above)
    @test isempty(rank_warnings(add_call(INTEGRATION_PANEL, ["x1", "x2"])))

    # A covariate that only varies within the cells of one treated state: the period slopes
    # and the slopes of the other states are not all determined
    one_state = copy(INTEGRATION_PANEL)
    one_state.z = [state == "s01" ? x1 : 0.0 for (state, x1) in zip(one_state.state, one_state.x1)]
    warnings = rank_warnings(add_call(one_state, ["x1", "z"]))
    @test length(warnings) == 1
    @test occursin("ccc = \"add\"", warnings[1])
    @test isempty(dropped_warnings(add_call(one_state, ["x1", "z"])))

    # Covariates that are collinear within cells
    collinear = copy(INTEGRATION_PANEL)
    collinear.z = 2 .* collinear.x1 .+ [sin(year + length(state)) for (state, year) in zip(collinear.state, collinear.year)]
    @test length(rank_warnings(add_call(collinear, ["x1", "z"]))) == 1

    # A covariate without any within-cell variation costs no rank: no warning
    cell_level = copy(INTEGRATION_PANEL)
    cell_level.z = [Float64(year) for year in cell_level.year]
    @test isempty(rank_warnings(add_call(cell_level, ["x1", "z"])))

    # The expected rank allows for a panel that falls apart into groups of states observed in
    # disjoint sets of periods: one rank less per covariate and group
    @test DiDInt.count_panel_components(["A", "A", "B", "B"], ["1", "2", "1", "2"]) == 1
    @test DiDInt.count_panel_components(["A", "A", "B", "C", "C"], ["1", "2", "2", "3", "4"]) == 2
    @test DiDInt.count_panel_components(["A", "B", "C"], ["1", "2", "3"]) == 3
    split_panel = filter(r -> (r.state <= "s06") == (r.year <= 2004), INTEGRATION_PANEL)
    @test isempty(rank_warnings(() -> first_stage_of(split_panel, "add", ["x1", "x2"])))
end

@testset "The true jackknife re-estimates with its own pruning, without repeating warnings" begin
    # x3 is collinear with x1 within cells; each jackknife replicate drops it as the full
    # sample does, so the jackknife standard error is the one with x1 alone
    panel = copy(INTEGRATION_PANEL)
    panel.x3 = 2 .* panel.x1 .+ [sin(year + length(state)) for (state, year) in zip(panel.state, panel.year)]
    for ccc in ["hom", "time"]
        call(covariates) = () -> didint_table("y", "state", "year", panel; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                                              ccc = ccc, agg = "cohort", covariates = covariates, truejack = true, nperm = 9, seed = 1234)
        @test length(dropped_warnings(call(["x1", "x3"]))) == 1
        with_x3 = quiet_call(call(["x1", "x3"]))
        only_x1 = quiet_call(call(["x1"]))
        @test with_x3.agg_att[1] ≈ only_x1.agg_att[1] rtol = 1e-9
        @test with_x3.jknifese_agg_att[1] ≈ only_x1.jknifese_agg_att[1] rtol = 1e-8
    end
end
