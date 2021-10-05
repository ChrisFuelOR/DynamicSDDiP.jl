using DynamicSDDiP
using Documenter

DocMeta.setdocmeta!(DynamicSDDiP, :DocTestSetup, :(using DynamicSDDiP); recursive=true)

makedocs(;
    modules=[DynamicSDDiP],
    authors="Christian Fuellner",
    repo="https://github.com/ChrisFuelOR/DynamicSDDiP.jl/blob/{commit}{path}#{line}",
    sitename="DynamicSDDiP.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ChrisFuelOR.github.io/DynamicSDDiP.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ChrisFuelOR/DynamicSDDiP.jl",
)
