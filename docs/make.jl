using Documenter
using DiDInt

makedocs(
    sitename = "DiDInt.jl",
    modules = [DiDInt],
    format = Documenter.HTML(
        size_threshold      = 200_000_000,  # 200 MB
        size_threshold_warn = 100_000_000,),
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md",
        "Details" => "details.md",
        "Examples" => "examples.md"
    ]
)

deploydocs(
    repo = "github.com/ebjamieson97/DiDInt.jl.git",
    devbranch = "main"
)