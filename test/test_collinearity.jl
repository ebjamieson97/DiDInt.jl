# Tests of how collinear covariates are handled in the first stage.
# Under ccc = "hom", "state", "time" and "int", a covariate that is a linear combination of
# covariates listed before it is dropped wherever the slopes are estimated, and the regular
# and the edge-case procedure use the same covariates. Under "hom", "state" and "time" a
# warning names the dropped covariates. Under "int" a warning is given only for the ATT(s,t)'s
# that depend on a dropped covariate. Under ccc = "add" no covariate is dropped, and a warning
# is given if the lambdas are not uniquely determined, saying whether the ATTs are affected.
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
atts_affected(warning) = occursin("ATT estimates are affected", warning)
atts_unaffected(warning) = occursin("ATT estimates are not affected", warning)

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

# Under ccc = "int": the ATT(s,t)'s that depend on a covariate dropped from a cell's regression,
# as tuples (treated state, year of t, state of the cell, year of the cell, covariates)
function dropped_effects(panel, covariates; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES)
    data = DataFrame(state_71X9yTx = String.(panel.state), time_71X9yTx = string.(Date.(panel.year)),
                     time_dmG5fpM = Date.(panel.year), outcome_71X9yTx = Float64.(panel.y))
    for cov in covariates
        data[!, cov] = Float64.(panel[!, cov])
    end
    dates = Date.(sort(unique(panel.year)))
    time_to_index = Dict(d => i for (i, d) in enumerate(dates))
    lambda_df, _, fit = DiDInt.first_stage(data, "int", covariates, "hc1", false)
    diff_df, _ = DiDInt.construct_diff_df(lambda_df, treated_states, Date.(treatment_times), time_to_index, dates, false, data)
    return [(e.state, year(e.t), e.cell[1], year(e.cell[2]), e.covariates) for e in DiDInt.dropped_covariate_effects(fit, diff_df)]
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
        # No covariate is dropped as collinear, and "add" finds no rank deficiency beyond the
        # expected one. Under "int" one cell cannot use x2, which is binary and takes a single
        # value there, and the ATT(s,t)'s that use that cell are reported.
        warnings = filter(w -> occursin("were dropped", w) || occursin("has rank", w), warnings_of(call))
        if ccc == "int"
            @test length(warnings) == 1
            @test occursin("x2 in state s12, period 2006-01-01", warnings[1]) && !occursin("x1", warnings[1])
        else
            @test isempty(warnings)
        end

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
    # The first stage itself does not warn. The ATT(s,t)'s that depend on the dropped covariate
    # are those that use the cell: s06 is a control, in period 2005 for the 2004 cohort and in
    # the base period 2005 of the 2006 cohort
    @test isempty(warnings_of(() -> first_stage_of(large, "int", ["x1", "x2"])))
    in_cell_effects(effects) = filter(e -> (e[3], e[4]) == ("s06", 2005), effects)
    effects = in_cell_effects(dropped_effects(large, ["x1", "x2"]))
    @test Set((e[1], e[2]) for e in effects) ==
          Set(vcat([(s, 2005) for s in ["s01", "s02"]], [(s, yr) for s in ["s03", "s04"] for yr in 2006:2008]))
    @test all(e -> e[5] == ["x2"], effects)
    @test all(e -> e[5] == ["x1"], in_cell_effects(dropped_effects(large, ["x2", "x1"])))
    @test isempty(in_cell_effects(dropped_effects(INTEGRATION_PANEL, ["x1", "x2"])))
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
    @test occursin("cat_w in state s06, period 2005-01-01", warnings[1]) && !occursin("cat_v", warnings[1])

    # A cell with fewer observations than covariates plus one keeps at most (observations - 1)
    # covariates, the ones listed first
    tiny = INTEGRATION_PANEL[.!in_cell .| (cumsum(in_cell) .<= 2), :]
    tiny.x2 = Float64.(tiny.x2)
    tiny.x3 = tiny.x1 .^ 2
    _, fit = quiet_call(() -> first_stage_of(tiny, "int", ["x1", "x2", "x3"]))
    @test fit.slopes[cell][1] != 0 && fit.slopes[cell][2:3] == [0, 0]
end

@testset "Covariates without within-cell variation are dropped" begin
    # A covariate measured at the state-period level, and one that only varies within the
    # cells of some states: neither is collinear with another covariate
    panel = copy(INTEGRATION_PANEL)
    cells = Dict(cell => i for (i, cell) in enumerate(unique(zip(panel.state, panel.year))))
    panel.unemp = [round(3 + 6 * abs(sin(12.9898 * cells[cell])), digits = 1) for cell in zip(panel.state, panel.year)]
    for ccc in ["hom", "state", "time", "int"]
        # unemp has slope 0 everywhere, and the other slopes are those without it
        _, fit = quiet_call(() -> first_stage_of(panel, ccc, ["x1", "unemp", "x2"]))
        _, without_fit = quiet_call(() -> first_stage_of(panel, ccc, ["x1", "x2"]))
        @test all(slopes[2] == 0 && slopes[[1, 3]] ≈ without_fit.slopes[key] for (key, slopes) in fit.slopes)
    end
    # Under "hom", "state" and "time" the cell dummies absorb it: no warning
    for ccc in ["hom", "state", "time"]
        call() = didint("y", "state", "year", panel; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES, ccc = ccc,
                        covariates = ["x1", "unemp", "x2"], nperm = 9, seed = 1234)
        @test isempty(dropped_warnings(call))
    end
    # Under "int" every cell other than the treated one is evaluated at the treated cell's
    # covariate means, where unemp takes another value: the ATT(s,t)'s depend on it
    effects = dropped_effects(panel, ["x1", "unemp", "x2"])
    @test any(e -> "unemp" in e[5], effects)
    @test !any(e -> "x1" in e[5], effects)

    # The pruning rule itself: a column of zeros is neither kept nor reported as collinear
    X = hcat([1.0, -1.0, 2.0, -2.0], zeros(4), [2.0, -2.0, 4.0, -4.0], [1.0, 1.0, -1.0, -1.0])
    @test DiDInt.prune_covariate_columns(X) == ([1, 4], [3])
    # ... and it does not depend on the units of the columns
    @test DiDInt.prune_covariate_columns(X .* [1e-6 1.0 1e6 1e3]) == ([1, 4], [3])
    # A column that is nearly, but not within the tolerance, a combination of others is kept
    @test DiDInt.prune_covariate_columns(hcat(X[:, 1], X[:, 1] .+ 1e-5 .* X[:, 4])) == ([1, 2], Int[])
    @test DiDInt.prune_covariate_columns(hcat(X[:, 1], X[:, 1] .+ 1e-9 .* X[:, 4])) == ([1], [2])
end

@testset "ccc = int: a warning only when a dropped covariate changes an ATT(s,t)" begin
    # A cell's slopes are evaluated at the covariate means of the treated cell. A covariate
    # dropped from the cell's regression matters only if those means lie outside what the
    # cell's own data determine. x1 varies in every cell; z is added to it.
    int_warnings(panel, covariates) = filter(w -> occursin("ccc = \"int\"", w), warnings_of() do
        didint("y", "state", "year", panel; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES, ccc = "int",
               covariates = covariates, nperm = 9, seed = 1234)
    end)
    n_att = nrow(input_didint(; ccc = "int", agg = "sgt", covariates = ["x1"]))
    state_number(s) = parse(Int, s[2:end])

    # No covariate is dropped: no warning
    @test isempty(dropped_effects(INTEGRATION_PANEL, ["x1"]))
    @test isempty(int_warnings(INTEGRATION_PANEL, ["x1"]))

    # A covariate that is constant in every cell and equal to the treated cell's mean everywhere:
    # dropped in every cell, but no ATT(s,t) depends on it
    same = copy(INTEGRATION_PANEL)
    same.z .= 4.0
    @test isempty(dropped_effects(same, ["x1", "z"]))
    @test isempty(int_warnings(same, ["x1", "z"]))
    @test input_didint(same; ccc = "int", agg = "sgt", covariates = ["x1", "z"]).att_sgt ≈
          input_didint(; ccc = "int", agg = "sgt", covariates = ["x1"]).att_sgt rtol = 1e-10 atol = 1e-12

    # ... and one that takes a different value in every cell: every ATT(s,t) depends on it
    differ = copy(INTEGRATION_PANEL)
    differ.z = Float64.(differ.year .- 2000) .+ state_number.(differ.state)
    effects = dropped_effects(differ, ["x1", "z"])
    @test length(unique((e[1], e[2]) for e in effects)) == n_att
    @test all(e -> e[5] == ["z"], effects)
    # the treated cell itself is evaluated at its own means and is never listed
    @test !any(e -> (e[1], e[2]) == (e[3], e[4]), effects)
    warnings = int_warnings(differ, ["x1", "z"])
    @test length(warnings) == 1
    @test occursin("$n_att ATT(s,t) estimates depend", warnings[1])
    @test occursin("z in state", warnings[1]) && !occursin("x1", warnings[1])

    # A covariate that only depends on the period: the cells of period t share the treated
    # cell's value, the cells of the base period do not
    by_period = copy(INTEGRATION_PANEL)
    by_period.z = Float64.(by_period.year)
    effects = dropped_effects(by_period, ["x1", "z"])
    @test !isempty(effects)
    @test all(e -> e[4] != e[2], effects)

    # Collinear covariates whose relation also holds at the treated cell's means: no warning,
    # whatever the units and however large the constant
    for (slope, constant) in [(2.0, 3.0), (2000.0, 3.0), (2.0, 1e6)]
        holds = copy(INTEGRATION_PANEL)
        holds.z = slope .* holds.x1 .+ constant
        @test isempty(dropped_effects(holds, ["x1", "z"]))
        @test isempty(int_warnings(holds, ["x1", "z"]))
    end

    # Collinear within every cell, but with a constant that differs from cell to cell: the
    # relation does not hold at the treated cell's means, and the covariate listed last is named
    breaks = copy(INTEGRATION_PANEL)
    breaks.z = 2 .* breaks.x1 .+ [sin(yr + 3 * state_number(s)) for (s, yr) in zip(breaks.state, breaks.year)]
    effects = dropped_effects(breaks, ["x1", "z"])
    @test length(unique((e[1], e[2]) for e in effects)) == n_att
    @test all(e -> e[5] == ["z"], effects)
    @test all(e -> e[5] == ["x1"], dropped_effects(breaks, ["z", "x1"]))
    @test length(int_warnings(breaks, ["x1", "z"])) == 1

    # The basis of what a cell leaves undetermined: a unit vector for a column of zeros, and
    # the coefficients of the linear combination for a collinear column
    X = hcat([1.0, -1.0, 2.0, -2.0], zeros(4), [2.0, -2.0, 4.0, -4.0], [1.0, 1.0, -1.0, -1.0])
    kept, _ = DiDInt.prune_covariate_columns(X)
    basis = DiDInt.null_space_basis(X, kept)
    @test basis ≈ [0.0 -2.0; 1.0 0.0; 0.0 1.0; 0.0 0.0] atol = 1e-12
    @test maximum(abs, X * basis) < 1e-12
end

@testset "ccc = add: rank deficiency beyond the split between states and periods" begin
    add_call(panel, covariates) = () -> didint("y", "state", "year", panel; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                                               ccc = "add", covariates = covariates, nperm = 9, seed = 1234)
    add_att(panel, covariates) = input_didint(panel; ccc = "add", agg = "sgt", covariates = covariates).att_sgt
    # A covariate that varies within every cell and is not related to x1
    noise = [2 * sin(12.9898 * i) for i in 1:nrow(INTEGRATION_PANEL)]
    with_z(z) = (panel = copy(INTEGRATION_PANEL); panel.z = Float64.(z); panel)

    # The ATT(s,t)'s of another least-squares solution of the first stage: the package's
    # solution with `move[state]` added to the slopes of that state. The move must leave every
    # fitted value unchanged (checked here), so both solutions fit equally well.
    function att_of_moved_solution(panel, covariates, move = Dict{String, Vector{Float64}}())
        data = DataFrame(state_71X9yTx = String.(panel.state), time_71X9yTx = string.(Date.(panel.year)),
                         time_dmG5fpM = Date.(panel.year), outcome_71X9yTx = Float64.(panel.y))
        for cov in covariates
            data[!, cov] = Float64.(panel[!, cov])
        end
        demeaned = DiDInt.with_demeaned_columns(data, covariates)
        additive = DiDInt.estimate_additive_slopes(demeaned, DiDInt.demeaned_name.(covariates); null_space = false)
        moved(state) = get(move, state, zeros(length(covariates)))
        X = Matrix{Float64}(demeaned[:, DiDInt.demeaned_name.(covariates)])
        @test maximum(abs(dot(moved(data.state_71X9yTx[i]), X[i, :])) for i in 1:nrow(data)) < 1e-8
        lambdas = combine(groupby(data, [:state_71X9yTx, :time_71X9yTx, :time_dmG5fpM])) do cell
            state = cell.state_71X9yTx[1]
            slopes = additive.gamma_states[:, additive.state_idx[state]] .+ moved(state) .+
                     additive.gamma_times[:, additive.time_idx[cell.time_71X9yTx[1]]]
            DataFrame(lambda = mean(cell.outcome_71X9yTx) - dot(slopes, [mean(cell[!, cov]) for cov in covariates]))
        end
        lambda_df = DataFrame(state = lambdas.state_71X9yTx, time = lambdas.time_dmG5fpM, lambda = lambdas.lambda,
                              lambda_index = 1:nrow(lambdas))
        dates = sort(unique(data.time_dmG5fpM))
        diff_df, _ = DiDInt.construct_diff_df(lambda_df, INTEGRATION_TREATED, Date.(INPUT_TIMES),
                                              Dict(d => i for (i, d) in enumerate(dates)), dates, false, data)
        return Float64.(DiDInt.build_att_st(diff_df, "both", "hc1", nothing).att)
    end
    # Without a move it is the package's estimate
    @test att_of_moved_solution(INTEGRATION_PANEL, ["x1", "x2"]) ≈ add_att(INTEGRATION_PANEL, ["x1", "x2"]) rtol = 1e-8 atol = 1e-10

    # Well-behaved covariates: no warning (the estimates are checked against version 0.10.0 above)
    @test isempty(rank_warnings(add_call(INTEGRATION_PANEL, ["x1", "x2"])))
    @test isempty(rank_warnings(add_call(with_z(noise), ["x1", "z"])))

    # A covariate that is constant within the cells of one state and over time. That state's
    # slope is not determined, and any value of it shifts all of the state's lambdas by the
    # same amount, which cancels in its long differences: the lambdas are not unique, the ATTs
    # are. This holds for a treated state and for a control state.
    for state in ["s01", "s07"]
        flat = with_z([s == state ? 5.0 : z for (s, z) in zip(INTEGRATION_PANEL.state, noise)])
        warnings = rank_warnings(add_call(flat, ["x1", "z"]))
        @test length(warnings) == 1
        @test occursin("ccc = \"add\"", warnings[1]) && occursin("not uniquely determined", warnings[1])
        @test atts_unaffected(warnings[1]) && !atts_affected(warnings[1])
        @test occursin("Covariates involved: z.", warnings[1])
        @test isempty(dropped_warnings(add_call(flat, ["x1", "z"])))
        # another solution, with any slope of z in that state, gives the same ATTs
        @test att_of_moved_solution(flat, ["x1", "z"], Dict(state => [0.0, 1.0])) ≈ add_att(flat, ["x1", "z"]) rtol = 1e-8 atol = 1e-10
    end

    # A covariate that is constant within each cell of one state but changes over time: the
    # lambdas of that state shift by different amounts, and the ATTs depend on the solution.
    # For a control state, every ATT(s,t) that uses it is affected.
    moving = with_z([s == "s07" ? Float64(yr - 2000)^2 : z for (s, yr, z) in zip(INTEGRATION_PANEL.state, INTEGRATION_PANEL.year, noise)])
    warnings = rank_warnings(add_call(moving, ["x1", "z"]))
    @test length(warnings) == 1
    @test atts_affected(warnings[1]) && !atts_unaffected(warnings[1])
    @test occursin("Covariates involved: z.", warnings[1])
    @test maximum(abs, att_of_moved_solution(moving, ["x1", "z"], Dict("s07" => [0.0, 1.0])) .- add_att(moving, ["x1", "z"])) > 1e-3

    # A covariate without any within-cell variation that differs across periods
    cell_level = with_z(Float64.(INTEGRATION_PANEL.year))
    warnings = rank_warnings(add_call(cell_level, ["x1", "z"]))
    @test length(warnings) == 1 && atts_affected(warnings[1])
    # ... and one that is the same in every cell: it shifts the lambdas by amounts that are
    # additive in state and period, so the ATTs are unaffected
    constant = with_z(fill(4.0, nrow(INTEGRATION_PANEL)))
    warnings = rank_warnings(add_call(constant, ["x1", "z"]))
    @test length(warnings) == 1 && atts_unaffected(warnings[1])
    @test add_att(constant, ["x1", "z"]) ≈ add_att(INTEGRATION_PANEL, ["x1"]) rtol = 1e-8 atol = 1e-10

    # Covariates that are collinear within cells: unaffected if the relation is the same
    # everywhere, affected if its constant differs from cell to cell
    same_relation = with_z(2 .* INTEGRATION_PANEL.x1 .+ 3)
    warnings = rank_warnings(add_call(same_relation, ["x1", "z"]))
    @test length(warnings) == 1 && atts_unaffected(warnings[1])
    @test occursin("Covariates involved: x1, z.", warnings[1])
    # moving a state's slopes along the collinear direction: z = 2 x1 + 3, so (-2, 1) fits equally well
    @test att_of_moved_solution(same_relation, ["x1", "z"], Dict("s07" => [-2.0, 1.0])) ≈ add_att(same_relation, ["x1", "z"]) rtol = 1e-8 atol = 1e-10
    cell_constants = with_z(2 .* INTEGRATION_PANEL.x1 .+ [sin(yr + 3 * parse(Int, s[2:end])) for (s, yr) in zip(INTEGRATION_PANEL.state, INTEGRATION_PANEL.year)])
    warnings = rank_warnings(add_call(cell_constants, ["x1", "z"]))
    @test length(warnings) == 1 && atts_affected(warnings[1])
    @test maximum(abs, att_of_moved_solution(cell_constants, ["x1", "z"], Dict("s07" => [-2.0, 1.0])) .- add_att(cell_constants, ["x1", "z"])) > 1e-3

    # The warning is given once, and not again by the true jackknife's re-estimations
    truejack_call() = didint("y", "state", "year", cell_level; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                             ccc = "add", covariates = ["x1", "z"], truejack = true, nperm = 9, seed = 1234)
    @test length(rank_warnings(truejack_call)) == 1

    # didint_plot() shows the lambdas themselves and says when they are not uniquely determined
    plot_call(panel) = () -> didint_plot("y", "state", "year", panel; treated_states = INTEGRATION_TREATED, treatment_times = INPUT_TIMES,
                                         ccc = "add", covariates = ["x1", "z"])
    @test isempty(rank_warnings(plot_call(with_z(noise))))
    warnings = rank_warnings(plot_call(cell_level))
    @test length(warnings) == 1 && occursin("not uniquely determined", warnings[1])

    # The change in the ATT(s,t)'s when the lambdas move: a shift that is the same in every
    # period of a state cancels, and a shift of a control's lambda in period t moves the ATTs
    # of that period by its share of the control weights
    lambda_df, fit = quiet_call(() -> first_stage_of(INTEGRATION_PANEL, "add", ["x1"]; treated_states = INTEGRATION_TREATED,
                                                     treatment_times = INPUT_TIMES))
    data = DataFrame(state_71X9yTx = String.(INTEGRATION_PANEL.state), time_dmG5fpM = Date.(INTEGRATION_PANEL.year))
    dates = Date.(sort(unique(INTEGRATION_PANEL.year)))
    diff_df, _ = DiDInt.construct_diff_df(lambda_df, INTEGRATION_TREATED, Date.(INPUT_TIMES), Dict(d => i for (i, d) in enumerate(dates)),
                                          dates, false, data)
    cells = collect(zip(lambda_df.state, lambda_df.time))
    by_state = Dict(cell => (cell[1] == "s07" ? 3.0 : cell[1] == "s01" ? -2.0 : 0.0) for cell in cells)
    @test maximum(abs, DiDInt.additive_att_shifts(by_state, diff_df, "both")) < 1e-12
    one_cell = Dict(cell => (cell == ("s07", Date(2006)) ? 1.0 : 0.0) for cell in cells)
    changes = DiDInt.additive_att_shifts(one_cell, diff_df, "none")
    treated_rows = diff_df[diff_df.treat .== 1, :]
    n_controls = length(setdiff(unique(INTEGRATION_PANEL.state), INTEGRATION_TREATED))
    # s07 is a control: period 2006 is t for the ATT(s,2006)'s and the base period of the 2007 cohort
    expected = [row.t == Date(2006) ? -1 / n_controls : row.r1 == Date(2006) ? 1 / n_controls : 0.0 for row in eachrow(treated_rows)]
    @test changes ≈ expected atol = 1e-12

    # The null space read from the sparse QR decomposition: Z N = 0, with as many columns as
    # the design lacks in rank
    Z = DiDInt.sparse([1.0 2.0 0.0 1.0; 2.0 4.0 0.0 0.5; 3.0 6.0 0.0 -1.0; 1.0 2.0 0.0 2.0])
    Zqr = qr(Z)
    N = DiDInt.sparse_qr_null_space(Zqr, 4)
    @test size(N) == (4, 2) && rank(N) == 2
    @test maximum(abs, Z * N) < 1e-12

    # The expected rank allows for a panel that falls apart into groups of states observed in
    # disjoint sets of periods: one rank less per covariate and group
    @test DiDInt.count_panel_components(["A", "A", "B", "B"], ["1", "2", "1", "2"]) == 1
    @test DiDInt.count_panel_components(["A", "A", "B", "C", "C"], ["1", "2", "2", "3", "4"]) == 2
    @test DiDInt.count_panel_components(["A", "B", "C"], ["1", "2", "3"]) == 3
    split_panel = filter(r -> (r.state <= "s06") == (r.year <= 2004), INTEGRATION_PANEL)
    @test isempty(rank_warnings(() -> first_stage_of(split_panel, "add", ["x1", "x2"])))
    _, split_fit = quiet_call(() -> first_stage_of(split_panel, "add", ["x1", "x2"]))
    @test isempty(split_fit.lambda_shifts)
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
