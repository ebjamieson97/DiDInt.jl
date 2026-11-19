# Draws the figures shown in README.md (docs/src/assets/readme_*.png) from the merit data.
# Run from the repository with: julia --project=docs docs/readme_figures.jl

using DiDInt
using CSV, DataFrames
using CairoMakie

treated_states = ["34", "57", "58", "59", "61", "64", "71", "72", "85", "88"]
treated_times = [2000, 1998, 1993, 1997, 1999, 1996, 1991, 1998, 1997, 2000]
test_data = CSV.read(joinpath(@__DIR__, "..", "test", "data", "merit.csv"), DataFrame)
test_data.state = string.(test_data.state)

plot_data = didint_plot("coll", "state", "year", test_data,
                        treatment_times = treated_times,
                        treated_states = treated_states,
                        covariates = ["asian", "black", "male"])

assets = joinpath(@__DIR__, "src", "assets")
figures = ["readme_event.png" => plot(plot_data; plot_type = "event"),
           "readme_conditional.png" => plot(plot_data; plot_type = "conditional"),
           "readme_conditional2.png" => plot(plot_data; plot_type = "conditional2", treated_states = "72"),
           "readme_integrated.png" => plot(plot_data; plot_type = "integrated"),
           "readme_integrated2.png" => plot(plot_data; plot_type = "integrated2")]
for (file, figure) in figures
    save(joinpath(assets, file), figure; px_per_unit = 1)
end
