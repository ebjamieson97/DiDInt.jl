# Snapshot tests: rerun every case in snapshot_cases.jl and compare all result columns with
# the stored outputs in test/snapshots/.
#
# The stored outputs are large and are not part of the repository, so these tests only run
# where the files exist (they are skipped on CI and for users of the package). They guard
# against unintended changes in results while developing: create the files once with
#   julia --project -e 'using Pkg; Pkg.test(test_args=["--update-snapshots"])'
# and regenerate them, stating the reason, whenever a change in results is intended.

include("snapshot_cases.jl")

if !isfile(SNAPSHOT_VERSION_FILE)
    @info "No stored snapshots in $(SNAPSHOT_DIR): skipping the snapshot tests."
else
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
end
