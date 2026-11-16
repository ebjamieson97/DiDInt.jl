using Documenter
using DiDInt

# The Changelog page is the repository's CHANGELOG.md, copied into the docs at build time.
# docs/src/changelog.md is generated here and is not part of the repository; its "edit"
# link points to CHANGELOG.md.
changelog = read(joinpath(@__DIR__, "..", "CHANGELOG.md"), String)
write(joinpath(@__DIR__, "src", "changelog.md"),
      """
      ```@meta
      EditURL = "https://github.com/ebjamieson97/DiDInt.jl/blob/main/CHANGELOG.md"
      ```

      """ * changelog)

makedocs(
    sitename = "DiDInt.jl",
    modules = [DiDInt],
    format = Documenter.HTML(
        size_threshold      = 200_000_000,  # 200 MB
        size_threshold_warn = 100_000_000,
        assets = ["assets/headings.css"],),
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md",
        "Details" => "details.md",
        "Examples" => "examples.md",
        "Changelog" => "changelog.md"
    ]
)

deploydocs(
    repo = "github.com/ebjamieson97/DiDInt.jl.git",
    devbranch = "main"
)
