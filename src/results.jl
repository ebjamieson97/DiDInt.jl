###########################################################################################
# results.jl: the object that didint() returns                                            #
#                                                                                         #
#   Section 1: DiDIntResult                                                               #
#   Section 2: Building it from the table of results                                      #
#   Section 3: Display                                                                    #
#   Section 4: Accessors (coef, stderror, coefnames) and conversion to a table            #
###########################################################################################


###########################################################################################
# Section 1: DiDIntResult                                                                 #
###########################################################################################

"""
    DiDIntResult

The results of [`didint`](@ref). Displaying it prints a summary of the specification and of the
estimates. Its fields hold the estimates as tables:

- `agg::DataFrame`: the top-level estimates, one row each. For most aggregations this is the single
  aggregate ATT. For the policy aggregations it is the ATT of each policy, preceded by the pooled ATT
  across policies when `policy_aggregation = true` (the column `policy_weight` then holds each policy's
  weight in the pooled ATT).
- `sub::DataFrame`: the sub-aggregate ATTs of the chosen `agg`, one row each, identified by the leading
  column(s) (for example `treatment_time` under `agg = "cohort"`). `weight` is the weight of each
  sub-aggregate in the aggregate ATT (or in its policy's ATT).
- `att_st::DataFrame`: the ATT of each treated state (`state`) in each period (`t`) from its treatment time
  (`gvar`) on, from which every other estimate is computed, with the number of control states behind it
  (`n_controls`) and the treated state's numbers of observations (`n_t` in period `t`; `n` in period `t`
  and in the period before treatment, `r1`, together).
- `specs::NamedTuple`: the specification that was estimated (`agg`, `ccc`, `weighting`, `nperm`, `seed`,
  the period length and date range, and so on).
- `table::DataFrame`: every result in one wide table, as passed to the Stata and R wrappers.

In `agg` and `sub`, the estimate is in `att` and its inference in `se` and `pval` (heteroskedasticity-consistent
standard error and p-value), `jknife_se` and `jknife_pval` (jackknife standard error and p-value) and `ri_pval`
(randomization inference p-value). Entries that could not be computed are `missing`.

# Methods
- `coef(result; level = :agg)`, `stderror(result; level = :agg, type = :hc)` and `coefnames(result; level = :agg)`
  return the estimates, their standard errors (`type = :hc` or `:jackknife`) and their labels as vectors, for the
  top-level estimates (`level = :agg`) or the sub-aggregate ATTs (`level = :sub`).
- `DataFrame(result)` returns the wide table (`result.table`). A `DiDIntResult` is also a Tables.jl table, so
  it can be passed directly to functions such as `CSV.write`.

# Example
```julia
result = didint("coll", "state", "year", data; treated_states = states, treatment_times = times)
result                 # prints a summary
result.agg             # the aggregate ATT with its standard errors and p-values
result.sub             # the sub-aggregate ATTs
coef(result)           # the aggregate ATT as a vector
DataFrame(result)      # everything in one wide table
```
"""
struct DiDIntResult
    agg::DataFrame
    sub::DataFrame
    att_st::DataFrame
    specs::NamedTuple
    table::DataFrame
end


###########################################################################################
# Section 2: Building it from the table of results                                        #
###########################################################################################

# Columns of the `agg` and `sub` tables, with the prefix that the same quantity carries in the
# wide table of results (for example "se" is "se_att_cohort" there)
const ESTIMATE_COLUMNS = ["att" => "", "se" => "se_", "pval" => "pval_", "jknife_se" => "jknifese_",
                          "jknife_pval" => "jknifepval_", "ri_pval" => "ri_pval_"]

# The estimate named `name` in the wide table (for example "agg_att") with its inference, taken
# from the given rows, as columns named as in ESTIMATE_COLUMNS
function estimate_columns(table::DataFrame, name::AbstractString, rows)
    return [Symbol(column) => Vector{Union{Missing, Float64}}(table[rows, prefix * name]) for (column, prefix) in ESTIMATE_COLUMNS]
end

function DiDIntResult(table::DataFrame, att_st::DataFrame, specs::NamedTuple)

    # Arrange the wide table of results (one row per sub-aggregate, top-level estimates on the
    # first row or on the first row of each policy) into the tables of a DiDIntResult
    agg_option = specs.agg
    spec = SUB_AGG_SPECS[agg_option]

    # Sub-aggregates: their key columns, then the estimate with its inference and its weight
    sub = select(table, [string(col) for (col, _, _) in spec.keys])
    for (column, values) in estimate_columns(table, spec.suffix, 1:nrow(table))
        sub[!, column] = values
    end
    sub.weight = Union{Missing, Float64}[isnothing(w) ? missing : w for w in table.weights]

    # Top-level estimates
    if !is_policy_agg(agg_option)
        agg = DataFrame(:estimate => ["Aggregate ATT"], estimate_columns(table, "agg_att", 1:1)...)
    else
        policies = unique(table.policy)
        first_rows = [findfirst(table.policy .== p) for p in policies]
        agg = DataFrame(:estimate => ["Policy " * p for p in policies], estimate_columns(table, "att_p", first_rows)...)
        if specs.policy_aggregation
            agg.policy_weight = Vector{Union{Missing, Float64}}(table.policy_weights[first_rows])
            pooled = DataFrame(:estimate => ["Pooled ATT"], estimate_columns(table, "att_pool", 1:1)...)
            pooled.policy_weight = Union{Missing, Float64}[missing]
            agg = vcat(pooled, agg)
        end
    end

    # The ATT of each treated state and period; the policy only matters for the policy aggregations
    att_st = is_policy_agg(agg_option) ? copy(att_st) : select(att_st, Not(:policy))

    return DiDIntResult(agg, sub, att_st, specs, table)
end


###########################################################################################
# Section 3: Display                                                                      #
###########################################################################################

const MODEL_LABELS = Dict("hom" => "Homogeneous DID-INT", "state" => "State-Varying DID-INT",
                          "time" => "Time-Varying DID-INT", "add" => "Two One-Way DID-INT",
                          "int" => "Two-Way Intersection DID-INT")

const GROUP_LABELS = Dict("cohort" => "treatment time", "simple" => "treatment time and period",
                          "state" => "treated state", "sgt" => "treated state and period",
                          "time" => "periods since treatment", "calendar" => "period",
                          "policy" => "policy", "policy_gt" => "policy, treatment time and period",
                          "policy_cohort" => "policy and treatment time", "policy_calendar" => "policy and period",
                          "policy_time" => "policy and periods since treatment")

# Number of sub-aggregate ATTs printed in the summary; the rest are in the `sub` field
const MAX_ROWS_SHOWN = 20

function Base.show(io::IO, result::DiDIntResult)

    # One-line form, used inside collections
    specs = result.specs
    print(io, "DiDIntResult(agg = \"", specs.agg, "\", ccc = \"", specs.ccc, "\", ", nrow(result.sub), " sub-aggregate ATT",
          nrow(result.sub) == 1 ? "" : "s")
    if !is_policy_agg(specs.agg)
        print(io, ", aggregate ATT = ", round(result.agg.att[1], sigdigits = 5))
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", result::DiDIntResult)

    # Summary: the specification, the top-level estimates and the sub-aggregate ATTs
    specs = result.specs
    show_table(table) = show(io, table; summary = false, eltypes = false, show_row_number = false, allrows = true, allcols = true)
    has_covariates = !isempty(specs.covariates)

    println(io, "DID-INT results")
    println(io, "  Outcome:          ", specs.outcome)
    println(io, "  Model:            ", has_covariates ? MODEL_LABELS[specs.ccc] * " (ccc = \"" * specs.ccc * "\")" : "No covariates")
    if has_covariates
        println(io, "  Covariates:       ", join(specs.covariates, ", "))
    end
    println(io, "  Aggregation:      \"", specs.agg, "\" (by ", GROUP_LABELS[specs.agg], ")")
    println(io, "  Weighting:        \"", specs.weighting, "\", treated_weight = \"", specs.treated_weight, "\"")
    println(io, "  Controls:         ", specs.notyet ? "never-treated and not-yet-treated states" : "never-treated states")
    println(io, "  Periods:          ", specs.period, ", from ", specs.start_date, " to ", specs.end_date)
    println(io, "  Standard errors:  ", uppercase(specs.hc), " (se) and jackknife (jknife_se)")
    println(io, "  Permutations:     ", ismissing(specs.nperm) ? "randomization inference was not performed" : string(specs.nperm, " (seed ", specs.seed, ")"))

    println(io)
    println(io, !is_policy_agg(specs.agg) ? "Aggregate ATT" : specs.policy_aggregation ? "Pooled ATT and ATT by policy" : "ATT by policy")
    show_table(result.agg)

    # For agg = "policy" the sub-aggregates are the policies' ATTs shown above
    if specs.agg != "policy"
        n_sub = nrow(result.sub)
        println(io)
        println(io)
        println(io, "Sub-aggregate ATTs by ", GROUP_LABELS[specs.agg], " (", n_sub, ")")
        show_table(first(result.sub, MAX_ROWS_SHOWN))
        if n_sub > MAX_ROWS_SHOWN
            println(io)
            print(io, "  ... and ", n_sub - MAX_ROWS_SHOWN, " more; see the `sub` field for all of them")
        end
    end
end


###########################################################################################
# Section 4: Accessors (coef, stderror, coefnames) and conversion to a table              #
###########################################################################################

function estimates_table(result::DiDIntResult, level::Symbol)

    # The table of top-level estimates (:agg) or of sub-aggregate ATTs (:sub)
    level in (:agg, :sub) || error("'level' must be :agg or :sub.")
    return level == :agg ? result.agg : result.sub
end

"""
    coef(result::DiDIntResult; level = :agg)

The ATT estimates as a vector: the top-level estimates (`level = :agg`, the aggregate ATT or, for the policy
aggregations, the pooled ATT and the ATT of each policy) or the sub-aggregate ATTs (`level = :sub`).
"""
StatsAPI.coef(result::DiDIntResult; level::Symbol = :agg) = copy(estimates_table(result, level).att)

"""
    stderror(result::DiDIntResult; level = :agg, type = :hc)

The standard errors of the estimates returned by `coef(result; level = level)`: heteroskedasticity-consistent
(`type = :hc`) or jackknife (`type = :jackknife`).
"""
function StatsAPI.stderror(result::DiDIntResult; level::Symbol = :agg, type::Symbol = :hc)
    type in (:hc, :jackknife) || error("'type' must be :hc or :jackknife.")
    table = estimates_table(result, level)
    return copy(type == :hc ? table.se : table.jknife_se)
end

"""
    coefnames(result::DiDIntResult; level = :agg)

Labels of the estimates returned by `coef(result; level = level)`.
"""
function StatsAPI.coefnames(result::DiDIntResult; level::Symbol = :agg)
    level == :agg && return copy(estimates_table(result, level).estimate)

    # A sub-aggregate is labelled by the values of its key columns
    keys = [string(col) for (col, _, _) in SUB_AGG_SPECS[result.specs.agg].keys]
    label(value) = value isa AbstractFloat && isinteger(value) ? string(Int(value)) : string(value)
    return [join((label(row[key]) for key in keys), "; ") for row in eachrow(estimates_table(result, level))]
end

# A DiDIntResult is a Tables.jl table: the wide table of all results. This is what makes
# DataFrame(result) and, for example, CSV.write(file, result) work.
Tables.istable(::Type{DiDIntResult}) = true
Tables.columnaccess(::Type{DiDIntResult}) = true
Tables.columns(result::DiDIntResult) = Tables.columns(result.table)
Tables.schema(result::DiDIntResult) = Tables.schema(result.table)
