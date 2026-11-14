# Snapshot test cases (update spec §12.1, Phase 0).
#
# Defines the datasets, the grid of didint() calls whose full outputs are recorded in
# test/snapshots/*.csv, and the helpers that flatten a results DataFrame into
# (case, row, column, value) rows and compare two such sets.
#
# Everything here is deterministic and independent of the global RNG and of the
# globals built in runtests.jl: the merit subset is hard-coded, and the synthetic
# panels use their own SplitMix64 generator so they are identical on every Julia version.
#
# Regenerate snapshots with
#   julia --project -e 'using Pkg; Pkg.test(test_args=["--update-snapshots"])'
# or only the cases whose id contains a substring with
#   julia --project -e 'using Pkg; Pkg.test(test_args=["--update-snapshots", "--only=merit_stag"])'

using DataFrames
using CSV

const SNAPSHOT_DIR = joinpath(@__DIR__, "snapshots")
const SNAPSHOT_NPERM = 49
const SNAPSHOT_SEED = 1234
const SNAPSHOT_RTOL = 1e-8
const SNAPSHOT_ATOL = 1e-10

# Columns that depend on Julia's RNG stream (RI draws via Random.seed!/shuffle). Their
# values are only compared on the Julia minor version that recorded the snapshots.
is_ri_column(col::AbstractString) = startswith(col, "ri_pval")

# ── Deterministic random numbers (SplitMix64) ─────────────────────────────────────────

mutable struct SplitMix64
    state::UInt64
end

function next_u64!(r::SplitMix64)
    r.state += 0x9e3779b97f4a7c15
    z = r.state
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# Uniform on [0, 1) from the top 53 bits
unif!(r::SplitMix64) = (next_u64!(r) >> 11) * 2.0^-53

# Standard normal via Box-Muller (1 - u keeps the log argument in (0, 1])
normal!(r::SplitMix64) = sqrt(-2 * log(1 - unif!(r))) * cos(2π * unif!(r))

# ── Synthetic panel generator ─────────────────────────────────────────────────────────

"""
    make_synthetic_panel(; n_states, years, cohorts, seed)

Individual-level panel with `n_states` states (labelled "s01", "s02", ...) observed in
every year of `years`. `cohorts` maps a state label to its treatment year; all other
states are never treated. Features (spec §12.1):
- cell sizes n_{s,t} that vary by state and period,
- state and year fixed effects plus a state-by-period shock (state-level noise),
- covariates `x1` (continuous) and `x2` (binary) whose distributions drift over time at
  state-specific rates, with a time-varying coefficient on `x1` so `ccc` matters,
- effect path 0.2 + 0.1 e in event time e = t - g ≥ 0.
"""
function make_synthetic_panel(; n_states::Int, years::AbstractVector{Int},
                              cohorts::Dict{String, Int}, seed::Integer)
    rng = SplitMix64(UInt64(seed))
    state_col = String[]; year_col = Int[]; y_col = Float64[]; x1_col = Float64[]; x2_col = Float64[]
    for s in 1:n_states
        label = "s" * lpad(string(s), 2, '0')
        alpha_s = 0.5 * normal!(rng)
        x1_drift = 0.15 * ((s % 4) - 1.5) / 1.5
        x2_drift = 0.04 * (s % 3)
        g = get(cohorts, label, nothing)
        for (k, yr) in enumerate(years)
            n_st = 6 + mod(3 * s + 5 * k, 11) + floor(Int, 8 * unif!(rng))
            shock_st = 0.05 * normal!(rng)
            beta_t = 0.1 * k
            effect = (isnothing(g) || yr < g) ? 0.0 : 0.2 + 0.1 * (yr - g)
            p_x2 = clamp(0.3 + x2_drift * (k - 1), 0.05, 0.95)
            for _ in 1:n_st
                x1 = 0.2 * s / n_states + x1_drift * (k - 1) + normal!(rng)
                x2 = unif!(rng) < p_x2 ? 1.0 : 0.0
                y = alpha_s + beta_t + (0.5 + 0.1 * k) * x1 + 0.3 * x2 + effect + shock_st + 0.5 * normal!(rng)
                push!(state_col, label); push!(year_col, yr); push!(y_col, y)
                push!(x1_col, x1); push!(x2_col, x2)
            end
        end
    end
    return DataFrame(state = state_col, year = year_col, y = y_col, x1 = x1_col, x2 = x2_col)
end

# ── Datasets ──────────────────────────────────────────────────────────────────────────

# The merit subset used by the existing tests (runtests.jl), hard-coded so it does not
# depend on the global RNG: the 10 controls drawn there by shuffle() under seed 1234,
# with the same edit to state 95 in 2000.
const MERIT_SNAPSHOT_CONTROLS = ["31", "12", "92", "91", "95", "73", "83", "33", "62", "74"]
const MERIT_SNAPSHOT_TREATED = ["34", "57", "58", "59", "61", "64", "71", "72", "85", "88"]
# Each treated state's first year with merit == 1 (TREATED_TIMES in runtests.jl)
const MERIT_SNAPSHOT_TIMES = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000]
# Common-adoption tests pass TREATED_TIMES[1] (state 34's year) for every treated state
const MERIT_SNAPSHOT_COMMON_TIME = 2000

function load_merit_snapshot_data()
    data = CSV.read(joinpath(@__DIR__, "data", "merit.csv"), DataFrame)
    data.state = string.(data.state)
    sel = (data.state .== "95") .& (data.year .== 2000)
    data[sel, :asian] .= 1
    data[sel, :male] .= 1
    data[sel, :black] .= 1
    keep = vcat(MERIT_SNAPSHOT_CONTROLS, MERIT_SNAPSHOT_TREATED)
    return filter(row -> row.state ∈ keep, data)
end

const SYNTH_YEARS = collect(2001:2008)
# Three cohorts (2004, 2006, 2007) and seven never-treated states
const SYNTH_STAG_COHORTS = Dict("s01" => 2004, "s02" => 2004, "s03" => 2006,
                                "s04" => 2006, "s05" => 2007)
# One cohort (2005) of four treated states and six never-treated states
const SYNTH_COMMON_COHORTS = Dict("s01" => 2005, "s02" => 2005, "s03" => 2005, "s04" => 2005)

"""
    snapshot_datasets()

Dict from dataset name to (data, outcome, keyword arguments that set the treatment and
covariates). Each entry describes one snapshot file `test/snapshots/<name>.csv`.
"""
function snapshot_datasets()
    merit = load_merit_snapshot_data()
    synth_stag = make_synthetic_panel(n_states = 12, years = SYNTH_YEARS,
                                      cohorts = SYNTH_STAG_COHORTS, seed = 20260101)
    synth_common = make_synthetic_panel(n_states = 10, years = SYNTH_YEARS,
                                        cohorts = SYNTH_COMMON_COHORTS, seed = 20260202)
    stag_states = sort(collect(keys(SYNTH_STAG_COHORTS)))
    common_states = sort(collect(keys(SYNTH_COMMON_COHORTS)))
    return Dict(
        "merit_stag" => (data = merit, outcome = "coll",
                         treat = (treated_states = MERIT_SNAPSHOT_TREATED,
                                  treatment_times = MERIT_SNAPSHOT_TIMES),
                         covariates = [:male, :asian, :black]),
        "merit_common" => (data = merit, outcome = "coll",
                           treat = (treated_states = MERIT_SNAPSHOT_TREATED,
                                    treatment_times = MERIT_SNAPSHOT_COMMON_TIME),
                           covariates = [:male, :asian, :black]),
        "synth_stag" => (data = synth_stag, outcome = "y",
                         treat = (treated_states = stag_states,
                                  treatment_times = [SYNTH_STAG_COHORTS[s] for s in stag_states]),
                         covariates = [:x1, :x2]),
        "synth_common" => (data = synth_common, outcome = "y",
                           treat = (treated_states = common_states,
                                    treatment_times = [SYNTH_COMMON_COHORTS[s] for s in common_states]),
                           covariates = [:x1, :x2]),
    )
end

# ── Case grid ─────────────────────────────────────────────────────────────────────────

const SNAPSHOT_AGGS = ["cohort", "simple", "state", "sgt", "time"]
const SNAPSHOT_WEIGHTINGS = ["both", "att", "diff", "none"]

"""
    snapshot_cases()

Vector of cases, each a NamedTuple (id, file, dataset, kwargs). The main grid (spec §12.1)
crosses every dataset with: no covariates (ccc irrelevant, run as "int") or covariates
with ccc ∈ ("int", "hom"); every `agg`; every `weighting`; `use_pre_controls` false/true.
The "extras" file adds `truejack = true` and `edgecase = true` runs.
"""
function snapshot_cases()
    cases = NamedTuple[]
    for dataset in ["merit_stag", "merit_common", "synth_stag", "synth_common"]
        for (covs, ccc) in [(false, "int"), (true, "int"), (true, "hom")]
            for agg in SNAPSHOT_AGGS, weighting in SNAPSHOT_WEIGHTINGS, upc in [false, true]
                id = join([dataset, covs ? "cov" : "nocov", "ccc=$ccc", "agg=$agg",
                           "w=$weighting", "upc=$upc"], "|")
                kwargs = (ccc = ccc, agg = agg, weighting = weighting, use_pre_controls = upc)
                push!(cases, (id = id, file = dataset, dataset = dataset, covariates = covs, kwargs = kwargs))
            end
        end
    end

    # truejack: the refit jackknife path, only used for ccc ∉ ("int", "state")
    for agg in SNAPSHOT_AGGS, upc in [false, true]
        id = join(["synth_stag", "cov", "ccc=hom", "agg=$agg", "w=both", "upc=$upc", "truejack"], "|")
        kwargs = (ccc = "hom", agg = agg, weighting = "both", use_pre_controls = upc, truejack = true)
        push!(cases, (id = id, file = "extras", dataset = "synth_stag", covariates = true, kwargs = kwargs))
    end

    # edgecase: two-row long-difference regressions (as in test_edgecase.jl)
    for ccc in ["int", "hom"], agg in SNAPSHOT_AGGS
        id = join(["ec_stag", "cov", "ccc=$ccc", "agg=$agg", "w=both", "upc=false", "edgecase"], "|")
        kwargs = (ccc = ccc, agg = agg, edgecase = true)
        push!(cases, (id = id, file = "extras", dataset = "ec_stag", covariates = true, kwargs = kwargs))
    end
    for ccc in ["int", "hom"], dataset in ["ec_common_1t", "ec_common_2t"]
        agg = dataset == "ec_common_1t" ? "cohort" : "state"
        id = join([dataset, "cov", "ccc=$ccc", "agg=$agg", "w=both", "upc=false", "edgecase"], "|")
        kwargs = (ccc = ccc, agg = agg, edgecase = true)
        push!(cases, (id = id, file = "extras", dataset = dataset, covariates = true, kwargs = kwargs))
    end
    return cases
end

# Datasets used only by the "extras" edgecase cases: merit states 88 and 71 (and 34) with
# one control, as in test_edgecase.jl, but taken from the full merit file directly.
function edgecase_datasets()
    full = CSV.read(joinpath(@__DIR__, "data", "merit.csv"), DataFrame)
    full.state = string.(full.state)
    sel = (full.state .== "95") .& (full.year .== 2000)
    full[sel, :asian] .= 1
    full[sel, :male] .= 1
    full[sel, :black] .= 1
    control = "11"
    covs = [:male, :asian, :black]
    return Dict(
        "ec_stag" => (data = filter(r -> r.state ∈ ["88", "71", control], full), outcome = "coll",
                      treat = (treated_states = ["88", "71"], treatment_times = [2000, 1991]),
                      covariates = covs),
        "ec_common_1t" => (data = filter(r -> r.state ∈ ["71", control], full), outcome = "coll",
                           treat = (treated_states = ["71"], treatment_times = [1991]),
                           covariates = covs),
        "ec_common_2t" => (data = filter(r -> r.state ∈ ["71", "34", control], full), outcome = "coll",
                           treat = (treated_states = ["71", "34"], treatment_times = [1991, 1991]),
                           covariates = covs),
    )
end

all_snapshot_datasets() = merge(snapshot_datasets(), edgecase_datasets())

# ── Running a case and flattening its output ──────────────────────────────────────────

"""
    run_snapshot_case(case, datasets)

Run didint() for one case with warnings and progress output suppressed. Returns the
results DataFrame, or the error message (first line) as a String if didint() throws.
"""
function run_snapshot_case(case, datasets)
    ds = datasets[case.dataset]
    covariates = case.covariates ? ds.covariates : nothing
    try
        return Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            DiDInt.didint(ds.outcome, "state", "year", ds.data;
                          ds.treat..., covariates = covariates, seed = SNAPSHOT_SEED,
                          nperm = SNAPSHOT_NPERM, verbose = false, case.kwargs...)
        end
    catch err
        return first(split(sprint(showerror, err), '\n'))
    end
end

# Text form of one cell: shortest round-trip repr for floats, explicit tokens otherwise
snapshot_string(x::Missing) = "missing"
snapshot_string(x::Nothing) = "nothing"
snapshot_string(x::AbstractFloat) = string(Float64(x))
snapshot_string(x) = string(x)

"""
    flatten_result(id, result)

Long-format DataFrame (case, row, column, value) of every cell of a results DataFrame,
or a single (case, 0, "error", message) row if the case threw.
"""
function flatten_result(id::AbstractString, result)
    if result isa AbstractString
        return DataFrame(case = [id], row = [0], column = ["error"], value = [result])
    end
    out = DataFrame(case = String[], row = Int[], column = String[], value = String[])
    for i in 1:nrow(result), col in names(result)
        push!(out, (id, i, col, snapshot_string(result[i, col])))
    end
    return out
end

snapshot_path(file::AbstractString) = joinpath(SNAPSHOT_DIR, file * ".csv")
const SNAPSHOT_VERSION_FILE = joinpath(SNAPSHOT_DIR, "julia_version.txt")

function read_snapshot_file(file::AbstractString)
    df = CSV.read(snapshot_path(file), DataFrame; types = Dict(:case => String, :row => Int,
                  :column => String, :value => String), pool = false, stringtype = String)
    return df
end

function write_snapshot_file(file::AbstractString, df::DataFrame)
    mkpath(SNAPSHOT_DIR)
    CSV.write(snapshot_path(file), df)
end

"""
    snapshot_values_match(expected, actual; compare_ri)

Compare two text cells: numbers with rtol/atol, NaN equal to NaN, anything else exactly.
"""
function snapshot_values_match(expected::AbstractString, actual::AbstractString)
    expected == actual && return true
    a = tryparse(Float64, expected)
    b = tryparse(Float64, actual)
    (isnothing(a) || isnothing(b)) && return false
    (isnan(a) || isnan(b)) && return isnan(a) && isnan(b)
    return isapprox(a, b; rtol = SNAPSHOT_RTOL, atol = SNAPSHOT_ATOL)
end

"""
    compare_snapshot(expected, actual; compare_ri = true)

Compare the long-format rows of one case. Returns a vector of human-readable mismatch
descriptions (empty when the case matches).
"""
function compare_snapshot(expected::AbstractDataFrame, actual::AbstractDataFrame; compare_ri::Bool = true)
    problems = String[]
    exp_cols = Set(expected.column); act_cols = Set(actual.column)
    for col in setdiff(exp_cols, act_cols)
        push!(problems, "column $col missing from output")
    end
    for col in setdiff(act_cols, exp_cols)
        push!(problems, "unexpected column $col in output")
    end
    exp_rows = maximum(expected.row); act_rows = maximum(actual.row)
    exp_rows != act_rows && push!(problems, "expected $exp_rows rows, got $act_rows")
    !isempty(problems) && return problems

    act_lookup = Dict((r.row, r.column) => r.value for r in eachrow(actual))
    for r in eachrow(expected)
        !compare_ri && is_ri_column(r.column) && continue
        got = act_lookup[(r.row, r.column)]
        if !snapshot_values_match(r.value, got)
            push!(problems, "row $(r.row), $(r.column): expected $(r.value), got $got")
        end
    end
    return problems
end
