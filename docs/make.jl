using Documenter, CoupledODETools

makedocs(;
    modules=[CoupledODETools],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jamesjscully/CoupledODETools.jl/blob/{commit}{path}#L{line}",
    sitename="CoupledODETools.jl",
    authors="Jack Scully",
    assets=String[],
)

deploydocs(;
    repo="github.com/jamesjscully/CoupledODETools.jl",
)
