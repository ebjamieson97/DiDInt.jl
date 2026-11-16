# Regression tests for construct_diff_df and edge_case_se.
# A small hand-built two-cohort design: A treated 2005, B treated 2007, C and D never
# treated, periods 2003-2008. Every lambda and every cell size is distinct, so a row
# computed at the wrong cell cannot match the hand-computed value by accident.

using Dates

diff_test_states = ["A", "B", "C", "D"]
diff_test_years = 2003:2008
diff_test_dates = [Date(y) for y in diff_test_years]
diff_test_time_to_index = Dict(d => i for (i, d) in enumerate(diff_test_dates))

# lambda_{s,t} = 100 s + t^2 and n_{s,t} = 2 s + t, with s and t the state and period positions
diff_test_lambda(s, t) = 100.0 * findfirst(==(s), diff_test_states) + diff_test_time_to_index[t]^2
diff_test_n(s, t) = 2 * findfirst(==(s), diff_test_states) + diff_test_time_to_index[t]

diff_test_lambda_df = DataFrame(state = String[], time = Date[], lambda = Float64[], lambda_index = Int[])
for s in diff_test_states, t in diff_test_dates
    push!(diff_test_lambda_df, (s, t, diff_test_lambda(s, t), nrow(diff_test_lambda_df) + 1))
end
diff_test_lambda_idx(s, t) = only(diff_test_lambda_df[(diff_test_lambda_df.state .== s) .& (diff_test_lambda_df.time .== t), :lambda_index])

# One row per observation, n_{s,t} rows per cell, in the columns construct_diff_df reads
diff_test_data = DataFrame(state_71X9yTx = String[], time_dmG5fpM = Date[])
for s in diff_test_states, t in diff_test_dates, _ in 1:diff_test_n(s, t)
    push!(diff_test_data, (s, t))
end

diff_test_treated = ["A", "B"]
diff_test_times = [Date(2005), Date(2007)]

function run_construct_diff_df(use_pre_controls)
    return DiDInt.construct_diff_df(diff_test_lambda_df, diff_test_treated, diff_test_times,
                                    diff_test_time_to_index, diff_test_dates, use_pre_controls,
                                    diff_test_data)
end

# Hand-computed long difference, cell sizes and lambda indices for state s in cell (t, r1)
function expected_diff_row(s, t, r1)
    return (diff = diff_test_lambda(s, t) - diff_test_lambda(s, r1),
            n = diff_test_n(s, t) + diff_test_n(s, r1),
            n_t = diff_test_n(s, t),
            lambda_idx = (diff_test_lambda_idx(s, r1), diff_test_lambda_idx(s, t)))
end

function row_matches_hand_computation(row)
    expected = expected_diff_row(row.state, row.t, row.r1)
    return row.diff ≈ expected.diff && row.n == expected.n && row.n_t == expected.n_t &&
           row.lambda_idx == expected.lambda_idx
end

@testset "ri_diff_df rows are computed at their own (t, r1) cell" begin
    diff_df, ri_diff_df = run_construct_diff_df(false)

    # A's placebo rows are B's cohort cells and vice versa
    a_rows = ri_diff_df[ri_diff_df.state .== "A", :]
    b_rows = ri_diff_df[ri_diff_df.state .== "B", :]
    @test Set(zip(a_rows.t, a_rows.r1)) == Set([(Date(2007), Date(2006)), (Date(2008), Date(2006))])
    @test Set(zip(b_rows.t, b_rows.r1)) == Set((Date(y), Date(2004)) for y in 2005:2008)
    @test all(a_rows.treated_time .== Date(2007))
    @test all(b_rows.treated_time .== Date(2005))

    # A few rows checked by hand: lambda_{A,2007} - lambda_{A,2006} = 25 - 16, and so on
    a_2007 = only(eachrow(a_rows[a_rows.t .== Date(2007), :]))
    @test a_2007.diff ≈ 9.0
    @test a_2007.n == (2 + 5) + (2 + 4)
    @test a_2007.n_t == 2 + 5
    b_2005 = only(eachrow(b_rows[b_rows.t .== Date(2005), :]))
    @test b_2005.diff ≈ 9.0 - 4.0
    @test b_2005.n == (4 + 3) + (4 + 2)
    b_2008 = only(eachrow(b_rows[b_rows.t .== Date(2008), :]))
    @test b_2008.diff ≈ 36.0 - 4.0

    # Every row of both tables matches the hand computation
    @test all(row_matches_hand_computation, eachrow(ri_diff_df))
    @test all(row_matches_hand_computation, eachrow(diff_df))
end

@testset "use_pre_controls moves correct not-yet-treated rows into diff_df" begin
    diff_df, ri_diff_df = run_construct_diff_df(true)

    # B is not yet treated in 2005 and 2006, so it is a control for A's cells at those t
    b_controls = diff_df[(diff_df.state .== "B") .& (diff_df.treat .== 0), :]
    @test Set(zip(b_controls.t, b_controls.r1)) == Set([(Date(2005), Date(2004)), (Date(2006), Date(2004))])
    @test all(row_matches_hand_computation, eachrow(b_controls))
    # B's rows at t >= 2007 for A's cells stay placebo-only
    @test Set(ri_diff_df[ri_diff_df.state .== "B", :t]) == Set([Date(2007), Date(2008)])
    @test all(row_matches_hand_computation, eachrow(diff_df))
    @test all(row_matches_hand_computation, eachrow(ri_diff_df))
end

@testset "lambda_idx order is (r1, t) for treated and control rows" begin
    diff_df, _ = run_construct_diff_df(false)
    for row in eachrow(diff_df)
        @test row.lambda_idx == (diff_test_lambda_idx(row.state, row.r1), diff_test_lambda_idx(row.state, row.t))
    end
end

@testset "edge_case_se equals the delta-method SE with nonzero cross-state covariances" begin
    # Var(ATT) for ATT = (lambda_{s,t} - lambda_{s,r1}) - (lambda_{c,t} - lambda_{c,r1}) is
    # a' V a with a = (-1, +1, +1, -1) on (s r1, s t, c r1, c t). A dense V makes every
    # cross-state covariance nonzero, as under ccc = "hom" or "time".
    M = [1.0 0.2 0.1 0.3; 0.4 1.5 0.2 0.1; 0.3 0.1 1.2 0.5; 0.2 0.3 0.4 0.9]
    V = M * M'
    a = [-1.0, 1.0, 1.0, -1.0]
    treated_idx = (1, 2)   # (r1, t) of the treated state
    control_idx = (3, 4)   # (r1, t) of the control state
    @test DiDInt.edge_case_se([treated_idx, control_idx], V) ≈ sqrt(a' * V * a)
    @test DiDInt.edge_case_se([control_idx, treated_idx], V) ≈ sqrt(a' * V * a)
end
