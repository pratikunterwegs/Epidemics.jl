using Epidemics
using Documenter

DocMeta.setdocmeta!(Epidemics, :DocTestSetup, :(using Epidemics); recursive=true)

makedocs(;
    modules=[Epidemics],
    authors="Pratik Gupte",
    repo="https://github.com/pratikunterwegs/Epidemics.jl/blob/{commit}{path}#{line}",
    sitename="Epidemics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pratikunterwegs.github.io/Epidemics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pratikunterwegs/Epidemics.jl",
    devbranch="main",
)
