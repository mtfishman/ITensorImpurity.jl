using ITensorImpurity
using Documenter

DocMeta.setdocmeta!(ITensorImpurity, :DocTestSetup, :(using ITensorImpurity); recursive=true)

makedocs(;
    modules=[ITensorImpurity],
    authors="Matthew Fishman <mfishman@flatironinstitute.org> and contributors",
    repo="https://github.com/mtfishman/ITensorImpurity.jl/blob/{commit}{path}#{line}",
    sitename="ITensorImpurity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mtfishman.github.io/ITensorImpurity.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mtfishman/ITensorImpurity.jl",
    devbranch="main",
)
