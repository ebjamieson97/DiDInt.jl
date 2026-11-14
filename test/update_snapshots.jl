# Regenerate the snapshot files in test/snapshots/ (update spec §12.1).
# Run through the test environment so CSV is available:
#   julia --project -e 'using Pkg; Pkg.test(test_args=["--update-snapshots"])'
#   julia --project -e 'using Pkg; Pkg.test(test_args=["--update-snapshots", "--only=synth_stag|cov"])'
# With --only=<substring>, only cases whose id contains the substring are rerun; every
# other case keeps its stored rows. State the reason for every changed snapshot.

include("snapshot_cases.jl")

function update_snapshots(; only::Union{Nothing, AbstractString} = nothing)
    datasets = all_snapshot_datasets()
    cases = snapshot_cases()
    selected = isnothing(only) ? cases : filter(c -> occursin(only, c.id), cases)
    isempty(selected) && error("No snapshot case id contains \"$only\".")
    println("Regenerating $(length(selected)) of $(length(cases)) snapshot cases")

    for file in unique(c.file for c in cases)
        file_cases = filter(c -> c.file == file, selected)
        isempty(file_cases) && continue
        new_ids = Set(c.id for c in file_cases)
        # Keep stored rows of cases that are not being regenerated
        kept = isfile(snapshot_path(file)) && !isnothing(only) ?
               filter(r -> !(r.case in new_ids), read_snapshot_file(file)) :
               DataFrame(case = String[], row = Int[], column = String[], value = String[])
        fresh = [flatten_result(c.id, run_snapshot_case(c, datasets)) for c in file_cases]
        combined = vcat(kept, fresh...)
        # Store cases in grid order so diffs stay readable
        order = Dict(c.id => i for (i, c) in enumerate(cases))
        combined.case_order = [order[id] for id in combined.case]
        sort!(combined, [:case_order, :row])
        select!(combined, Not(:case_order))
        write_snapshot_file(file, combined)
        println("  wrote $(snapshot_path(file)) ($(length(file_cases)) cases regenerated)")
    end

    open(SNAPSHOT_VERSION_FILE, "w") do io
        println(io, "$(VERSION.major).$(VERSION.minor)")
    end
end

only_arg = findfirst(a -> startswith(a, "--only="), ARGS)
update_snapshots(only = isnothing(only_arg) ? nothing : ARGS[only_arg][length("--only=")+1:end])
