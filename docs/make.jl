using FuzzyCRegression
using Documenter

DocMeta.setdocmeta!(FuzzyCRegression, :DocTestSetup, :(using FuzzyCRegression); recursive=true)

makedocs(;
    modules=[FuzzyCRegression],
    authors="Aidan Toner-Rodgers",
    sitename="FuzzyCRegression.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aidantr.github.io/FuzzyCRegression.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "package.md",
        "API" => "API.md",
        
    ],
)

deploydocs(;
    repo="github.com/aidantr/FuzzyCRegression.jl.git",
    target = "build",
    branch = "gh-pages",
)


