# Snapshot tests: rerun every case in snapshot_cases.jl and compare
# all result columns with the stored outputs in test/snapshots/. A later phase may change
# a snapshot only with a stated reason; regenerate with update_snapshots.jl.

include("snapshot_cases.jl")

# RI p-values depend on Julia's RNG stream, which is not guaranteed stable across Julia
# versions, so they are only compared on the minor version that recorded the snapshots.
recorded_version = strip(read(SNAPSHOT_VERSION_FILE, String))
compare_ri = recorded_version == "$(VERSION.major).$(VERSION.minor)"
if !compare_ri
    @info "Snapshots were recorded on Julia $recorded_version; skipping ri_pval_* columns on Julia $VERSION."
end

snapshot_data = all_snapshot_datasets()
snapshot_case_list = snapshot_cases()
stored = Dict(file => groupby(read_snapshot_file(file), :case)
              for file in unique(c.file for c in snapshot_case_list))

for case in snapshot_case_list
    @testset "$(case.id)" begin
        groups = stored[case.file]
        @test haskey(groups, (case = case.id,))
        if haskey(groups, (case = case.id,))
            expected = groups[(case = case.id,)]
            actual = flatten_result(case.id, run_snapshot_case(case, snapshot_data))
            problems = compare_snapshot(expected, actual; compare_ri = compare_ri)
            isempty(problems) || println("Snapshot mismatch in $(case.id):\n  ", join(first(problems, 10), "\n  "))
            @test isempty(problems)
        end
    end
end
