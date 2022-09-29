using FuzzyCRegression
using Documenter

DocMeta.setdocmeta!(FuzzyCRegression, :DocTestSetup, :(using FuzzyCRegression); recursive=true)

makedocs(;
    modules=[FuzzyCRegression],
    authors="Aidan Toner-Rodgers",
    repo="https://github.com/aidantr/FuzzyCRegression.jl/blob/{commit}{path}#{line}",
    sitename="FuzzyCRegression.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aidantr.github.io/FuzzyCRegression.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aidantr/FuzzyCRegression.jl",
    devbranch="main",
)
