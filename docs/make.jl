using Documenter
using DiDInt

makedocs(
    sitename = "DiDInt.jl",
    modules = [DiDInt],
    format = Documenter.HTML(),
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