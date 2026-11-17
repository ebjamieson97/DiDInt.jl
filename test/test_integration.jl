# Tests of the integrated estimator.
# The panel is test/data/didintpy_panel.csv: 12 states observed 2001-2008, cohorts 2004 (s01, s02),
# 2006 (s03, s04) and 2007 (s05), covariates x1 (continuous) and x2 (binary) whose distributions
# drift over time at state-specific rates.

using Dates
using LinearAlgebra
using Random
using Statistics

const INTEGRATION_PANEL = CSV.read(joinpath(@__DIR__, "data", "didintpy_panel.csv"), DataFrame)
const INTEGRATION_COHORTS = Dict("s01" => 2004, "s02" => 2004, "s03" => 2006, "s04" => 2006, "s05" => 2007)
const INTEGRATION_TREATED = sort(collect(keys(INTEGRATION_COHORTS)))
const INTEGRATION_COVARIATES = ["x1", "x2"]

function integration_didint(panel = INTEGRATION_PANEL; cohorts = INTEGRATION_COHORTS, kwargs...)
    treated = sort(collect(keys(cohorts)))
    return Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        didint_table("y", "state", "year", panel; treated_states = treated,
                      treatment_times = [cohorts[s] for s in treated], covariates = INTEGRATION_COVARIATES,
                      nperm = 19, seed = 1234, kwargs...)
    end
end

# The first-stage inputs as didint() prepares them, to call the internal stages directly
function first_stage_inputs(panel)
    data = DataFrame(state_71X9yTx = String.(panel.state), time_71X9yTx = string.(Date.(panel.year)),
                     time_dmG5fpM = Date.(panel.year), outcome_71X9yTx = Float64.(panel.y),
                     x1 = Float64.(panel.x1), x2 = Float64.(panel.x2))
    dates = Date.(sort(unique(panel.year)))
    return data, dates, Dict(d => i for (i, d) in enumerate(dates))
end

cell_rows(panel, s, yr) = panel[(panel.state .== s) .& (panel.year .== yr), :]
cell_xbar(panel, s, yr) = [mean(cell_rows(panel, s, yr)[!, cov]) for cov in INTEGRATION_COVARIATES]

# OLS of y on [1, covariates] within one cell, dropping covariates that are constant in the cell.
# Returns the coefficients in full length (0 for a dropped covariate), their HC covariance in
# full size (zero rows and columns for a dropped covariate), and the residuals.
function cell_ols(panel, s, yr; hc_scale = "hc1")
    rows = cell_rows(panel, s, yr)
    X_all = Matrix(Float64.(rows[:, INTEGRATION_COVARIATES]))
    kept = [j for j in 1:size(X_all, 2) if length(unique(X_all[:, j])) > 1]
    X = hcat(ones(nrow(rows)), X_all[:, kept])
    beta = X \ rows.y
    resid = rows.y .- X * beta
    n, k = size(X)
    scale = hc_scale == "hc1" ? n / (n - k) : 1.0
    XXinv = inv(X' * X)
    V = XXinv * (X' * ((scale .* resid .^ 2) .* X)) * XXinv
    positions = vcat(1, 1 .+ kept)
    beta_full = zeros(1 + size(X_all, 2)); beta_full[positions] = beta
    V_full = zeros(1 + size(X_all, 2), 1 + size(X_all, 2)); V_full[positions, positions] = V
    return (beta = beta_full, V = V_full, resid = resid, X = X_all, y = Vector{Float64}(rows.y))
end

@testset "ATT(s,t) and aggregates match didintpy's integrated estimator" begin
    # Reference numbers from didintpy's fitted_didint.run_fitted_core
    # (weighting = "both"), stored in test/data/didintpy_reference.csv
    reference = CSV.read(joinpath(@__DIR__, "data", "didintpy_reference.csv"), DataFrame)
    for ccc in ["hom", "state", "time", "int"], upc in [false, true]
        ref = reference[(reference.ccc .== ccc) .& (reference.not_yet_treated .== upc), :]
        level(name) = ref[ref.level .== name, :]
        estimate(agg) = integration_didint(ccc = ccc, agg = agg, weighting = "both", use_pre_controls = upc)

        sgt = estimate("sgt")
        ref_st = level("att_st")
        @test nrow(sgt) == nrow(ref_st)
        for r in eachrow(ref_st)
            ours = sgt.att_sgt[(sgt.state .== r.state) .& (year.(sgt.t) .== r.t)]
            @test only(ours) ≈ r.estimate atol = 1e-8
        end

        simple = estimate("simple")
        for r in eachrow(level("att_gt"))
            @test only(simple.att_gt[(year.(simple.gvar) .== r.g) .& (year.(simple.time) .== r.t)]) ≈ r.estimate atol = 1e-8
        end

        cohort = estimate("cohort")
        for r in eachrow(level("att_g"))
            @test only(cohort.att_cohort[year.(cohort.treatment_time) .== r.g]) ≈ r.estimate atol = 1e-8
        end
        @test cohort.agg_att[1] ≈ only(level("att").estimate) atol = 1e-8

        state = estimate("state")
        for r in eachrow(level("att_s"))
            @test only(state.att_s[state.state .== r.state]) ≈ r.estimate atol = 1e-8
        end
    end
end

@testset "Integration changes ATT(s,t) only under ccc = int" begin
    data, dates, time_to_index = first_stage_inputs(INTEGRATION_PANEL)
    times = [Date(INTEGRATION_COHORTS[s]) for s in INTEGRATION_TREATED]
    labels = Dict(zip(INTEGRATION_TREATED, times))
    for ccc in ["hom", "state", "time", "add", "int"], upc in [false, true]
        lambda_df, vcov_lambda, fit = DiDInt.first_stage(data, ccc, INTEGRATION_COVARIATES, "hc1", false)
        diff_df, ri_diff_df = DiDInt.construct_diff_df(lambda_df, INTEGRATION_TREATED, times, time_to_index, dates, upc, data)
        plain = DiDInt.build_att_st(diff_df, "both", "hc1", vcov_lambda)
        integrated = DiDInt.build_att_st(diff_df, "both", "hc1", vcov_lambda; integration = fit)

        # Size of the covariate terms that integration adds to the treated states' long differences
        slope_terms = [DiDInt.slope_change_term(fit, r.state, r.t, r.r1, r.state) for r in eachrow(plain)]
        if ccc in ["hom", "state"]
            # Slopes do not change over time: the term is zero in every long difference
            @test maximum(abs.(slope_terms)) < 1e-10
        else
            @test maximum(abs.(slope_terms)) > 1e-3
        end

        if ccc == "int"
            @test maximum(abs.(integrated.att .- plain.att)) > 1e-3
        else
            # For "time" and "add" the terms are nonzero but cancel between treated and controls
            @test integrated.att ≈ plain.att atol = 1e-10
            @test integrated.se ≈ plain.se atol = 1e-10
        end

        # The batched engine used by the jackknife and RI reproduces the integrated ATT(s,t)'s
        long_diffs = vcat(diff_df, ri_diff_df)
        states = unique(long_diffs.state)
        engine = DiDInt.build_engine(long_diffs, DiDInt.att_st_cells(integrated), states, "both"; integration = fit)
        L, P = DiDInt.assignment_row(states, labels, Dict(s => DiDInt.DEFAULT_POLICY for s in INTEGRATION_TREATED),
                                     Dict(DiDInt.DEFAULT_POLICY => 1))
        sgt_rows = DataFrame(state = integrated.state, gvar = integrated.gvar, t = integrated.t)
        est = DiDInt.engine_estimates(engine, sgt_rows, "sgt", "both", upc, reshape(L, 1, :), reshape(P, 1, :))
        @test vec(est.sub) ≈ Float64.(integrated.att) atol = 1e-10
    end

    # didint() itself: every ccc except "int" gives the same estimates with and without the
    # integration switch, so only "int" is integrated; without covariates "add" runs and
    # every ccc gives the cell-mean estimator
    no_covariates(ccc) = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        didint_table("y", "state", "year", INTEGRATION_PANEL; treated_states = INTEGRATION_TREATED,
                      treatment_times = [INTEGRATION_COHORTS[s] for s in INTEGRATION_TREATED],
                      ccc = ccc, agg = "sgt", nperm = 19, seed = 1234).att_sgt
    end
    @test no_covariates("add") ≈ no_covariates("int")
    @test no_covariates("time") ≈ no_covariates("int")
end

@testset "Integrated long differences match a direct computation of DID-INT eqs. 46-47" begin
    data, dates, _ = first_stage_inputs(INTEGRATION_PANEL)
    states = sort(unique(INTEGRATION_PANEL.state))
    years = sort(unique(INTEGRATION_PANEL.year))
    n_cov = length(INTEGRATION_COVARIATES)

    # Direct first stage: one regression of y on cell dummies and the covariates interacted
    # with period dummies ("time") or with state dummies and period dummies ("add")
    function dummy_regression(ccc)
        cells = [(s, yr) for s in states for yr in years]
        cell_index = Dict(c => i for (i, c) in enumerate(cells))
        n = nrow(INTEGRATION_PANEL)
        n_slopes = ccc == "time" ? length(years) * n_cov : (length(states) + length(years)) * n_cov
        Z = zeros(n, length(cells) + n_slopes)
        for (i, row) in enumerate(eachrow(INTEGRATION_PANEL))
            Z[i, cell_index[(row.state, row.year)]] = 1.0
            s_i = findfirst(==(row.state), states); t_i = findfirst(==(row.year), years)
            for (j, cov) in enumerate(INTEGRATION_COVARIATES)
                if ccc == "time"
                    Z[i, length(cells) + (t_i - 1) * n_cov + j] = row[cov]
                else
                    Z[i, length(cells) + (s_i - 1) * n_cov + j] = row[cov]
                    Z[i, length(cells) + (length(states) + t_i - 1) * n_cov + j] = row[cov]
                end
            end
        end
        beta = pinv(Z) * INTEGRATION_PANEL.y
        lambda(s, yr) = beta[cell_index[(s, yr)]]
        function gamma(s, yr)
            s_i = findfirst(==(s), states); t_i = findfirst(==(yr), years)
            time_part = beta[length(cells) .+ ((ccc == "time" ? 0 : length(states)) + t_i - 1) * n_cov .+ (1:n_cov)]
            return ccc == "time" ? time_part : time_part .+ beta[length(cells) .+ (s_i - 1) * n_cov .+ (1:n_cov)]
        end
        return lambda, gamma
    end

    # Direct first stage for "int": a separate regression in every cell
    int_lambda(s, yr) = cell_ols(INTEGRATION_PANEL, s, yr).beta[1]
    int_gamma(s, yr) = cell_ols(INTEGRATION_PANEL, s, yr).beta[2:end]

    for ccc in ["time", "add", "int"]
        lambda_df, _, fit = DiDInt.first_stage(data, ccc, INTEGRATION_COVARIATES, "hc1", false)
        package_lambda(s, yr) = only(lambda_df.lambda[(lambda_df.state .== s) .& (lambda_df.time .== Date(yr))])
        lambda, gamma = ccc == "int" ? (int_lambda, int_gamma) : dummy_regression(ccc)

        # Long differences of the treated state itself (eq. 46) and of controls evaluated at
        # the treated state's covariate means (eq. 47)
        for (target, source, t, r1) in [("s01", "s01", 2006, 2003), ("s01", "s07", 2006, 2003), ("s03", "s11", 2008, 2005),
                                         ("s05", "s05", 2007, 2006), ("s05", "s09", 2008, 2006)]
            direct = (lambda(source, t) - lambda(source, r1)) +
                     dot(gamma(source, t) .- gamma(source, r1), cell_xbar(INTEGRATION_PANEL, target, t))
            package = (package_lambda(source, t) - package_lambda(source, r1)) +
                      DiDInt.slope_change_term(fit, source, Date(t), Date(r1), target)
            @test package ≈ direct atol = 1e-8
        end
    end

    # didint() under "int": ATT(s,t) is the treated state's integrated long difference minus the
    # n-weighted mean of the never-treated states' integrated long differences (eq. 48)
    sgt = integration_didint(ccc = "int", agg = "sgt", weighting = "both")
    controls = setdiff(states, INTEGRATION_TREATED)
    predict(s, yr, xbar) = int_lambda(s, yr) + dot(int_gamma(s, yr), xbar)
    for row in eachrow(sgt[in.(sgt.state, Ref(["s01", "s05"])), :])
        t, r1 = year(row.t), year(row.gvar) - 1
        xbar = cell_xbar(INTEGRATION_PANEL, row.state, t)
        long_diff(k) = predict(k, t, xbar) - predict(k, r1, xbar)
        weights = [nrow(cell_rows(INTEGRATION_PANEL, c, t)) + nrow(cell_rows(INTEGRATION_PANEL, c, r1)) for c in controls]
        @test row.att_sgt ≈ long_diff(row.state) - sum(weights .* long_diff.(controls)) / sum(weights) atol = 1e-8
    end
end

@testset "Edge-case SE under int with covariates: variance conditional on the covariates" begin
    # One treated state (s01, treated 2004) and one control (s06)
    panel = INTEGRATION_PANEL[in.(INTEGRATION_PANEL.state, Ref(["s01", "s06"])), :]
    cohorts = Dict("s01" => 2004)
    for hc in ["hc0", "hc1"]
        edge = integration_didint(panel; cohorts = cohorts, ccc = "int", agg = "sgt", edgecase = true, hc = hc)
        no_edge = integration_didint(panel; cohorts = cohorts, ccc = "int", agg = "sgt", edgecase = false, hc = hc)
        @test all(ismissing, no_edge.se_att_sgt)
        @test edge.att_sgt ≈ no_edge.att_sgt atol = 1e-10
        for row in eachrow(edge)
            t, r1 = year(row.t), 2003
            cells = [cell_ols(panel, "s01", t; hc_scale = hc), cell_ols(panel, "s01", r1; hc_scale = hc),
                     cell_ols(panel, "s06", t; hc_scale = hc), cell_ols(panel, "s06", r1; hc_scale = hc)]
            xbar = cell_xbar(panel, "s01", t)

            # The ATT is the combination of the four cells' fitted values at xbar
            fitted = [cell.beta[1] + dot(cell.beta[2:end], xbar) for cell in cells]
            @test row.att_sgt ≈ fitted[1] - fitted[2] - fitted[3] + fitted[4] atol = 1e-8

            # Var(ATT(s,t) | X) = sum over the four cells of
            #   Var(lambda_k) + 2 xbar' Cov(gamma_k, lambda_k) + xbar' Var(gamma_k) xbar
            variance = sum(cell.V[1, 1] + 2 * dot(xbar, cell.V[2:end, 1]) + xbar' * cell.V[2:end, 2:end] * xbar for cell in cells)
            @test row.se_att_sgt ≈ sqrt(variance) rtol = 1e-8
            # An edge-case SE has no p-value
            @test ismissing(row.pval_att_sgt) || isnan(row.pval_att_sgt)
        end
    end
end

@testset "Edge-case SE under int without covariates: the sum of the four Var(lambda)" begin
    # Without covariates the conditional variance reduces to the sum of the variances of the
    # four cell means, the formula used for every ccc before covariates enter
    panel = INTEGRATION_PANEL[in.(INTEGRATION_PANEL.state, Ref(["s01", "s06"])), :]
    cohorts = Dict("s01" => 2004)
    function mean_variance(s, yr, hc)
        y = cell_rows(panel, s, yr).y
        n = length(y)
        return (hc == "hc1" ? n / (n - 1) : 1.0) * sum((y .- mean(y)) .^ 2) / n^2
    end
    for hc in ["hc0", "hc1"]
        edge = integration_didint(panel; cohorts = cohorts, ccc = "int", agg = "sgt", edgecase = true, hc = hc, covariates = nothing)
        for row in eachrow(edge)
            t, r1 = year(row.t), 2003
            variance = sum(mean_variance(s, yr, hc) for (s, yr) in [("s01", t), ("s01", r1), ("s06", t), ("s06", r1)])
            @test row.se_att_sgt ≈ sqrt(variance) rtol = 1e-8
        end
        # ... and under HC0 it is the same for every ccc, as none of them has covariates to
        # adjust for. (Under HC1 the degrees-of-freedom factor is that of the regression that
        # estimates the cells together: one per cell for "int", fewer for the others.)
        if hc == "hc0"
            for ccc in ["hom", "state", "time"]
                other = integration_didint(panel; cohorts = cohorts, ccc = ccc, agg = "sgt", edgecase = true, hc = hc, covariates = nothing)
                @test other.se_att_sgt ≈ edge.se_att_sgt rtol = 1e-8
            end
        end
    end

    # The function itself: with covariance matrices that only have an intercept entry, the
    # result is the square root of the sum of the four intercept variances
    fit = DiDInt.CovariateFit(["x1"])
    variances = Dict(("T", Date(2002)) => 0.4, ("T", Date(2001)) => 0.3, ("C", Date(2002)) => 0.2, ("C", Date(2001)) => 0.1)
    for (cell, variance) in variances
        fit.cell_vcov[cell] = [variance 0.0; 0.0 0.0]
        fit.xbar[cell] = [5.0]
    end
    @test DiDInt.edge_case_se_integrated(fit, "T", "C", Date(2002), Date(2001)) ≈ sqrt(1.0)
    # A missing cell gives a missing SE
    delete!(fit.cell_vcov, ("C", Date(2001)))
    @test ismissing(DiDInt.edge_case_se_integrated(fit, "T", "C", Date(2002), Date(2001)))
end

@testset "Edge-case SE under int tracks the sampling SD of ATT(s,t) given the covariates" begin
    # A short Monte Carlo. Two states (T treated in 2002, C control), two periods, and
    # cell-specific slopes that do not cancel. The covariates are drawn once and held fixed;
    # only the errors are redrawn, so the sampling SD is the one conditional on the covariates.
    rng = Random.MersenneTwister(20260115)
    cells = [("T", 2001, 0.0, 0.5, 0.0), ("T", 2002, 1.0, 1.5, 0.8), ("C", 2001, 0.3, 1.0, 0.2), ("C", 2002, 0.6, 0.2, 0.5)]
    design = DataFrame(state = String[], year = Int[], x1 = Float64[], x2 = Float64[], signal = Float64[])
    for (s, yr, intercept, slope, x_mean) in cells, _ in 1:60
        x1 = x_mean + randn(rng)
        x2 = randn(rng)
        push!(design, (s, yr, x1, x2, intercept + slope * x1 + 0.3 * x2))
    end
    function draw_panel()
        panel = copy(design)
        panel.y = panel.signal .+ randn(rng, nrow(panel))
        return panel
    end
    function direct_att(panel)
        fit(s, yr) = cell_ols(panel, s, yr).beta
        a = vcat(1.0, cell_xbar(panel, "T", 2002))
        return dot(a, fit("T", 2002)) - dot(a, fit("T", 2001)) - (dot(a, fit("C", 2002)) - dot(a, fit("C", 2001)))
    end

    atts = [direct_att(draw_panel()) for _ in 1:600]
    sampling_sd = std(atts)

    ses = Float64[]
    for _ in 1:15
        panel = draw_panel()
        result = integration_didint(panel; cohorts = Dict("T" => 2002), ccc = "int", agg = "sgt", edgecase = true)
        @test only(result.att_sgt) ≈ direct_att(panel) atol = 1e-8
        push!(ses, only(result.se_att_sgt))
    end
    @test 0.85 < mean(ses) / sampling_sd < 1.15
end
