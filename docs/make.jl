using FDMT
using Documenter

DocMeta.setdocmeta!(FDMT, :DocTestSetup, :(using FDMT); recursive=true)

makedocs(;
    modules=[FDMT],
    authors="Kiran Shila <me@kiranshila.com> and contributors",
    repo="https://github.com/kiranshila/FDMT.jl/blob/{commit}{path}#{line}",
    sitename="FDMT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kiranshila.github.io/FDMT.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kiranshila/FDMT.jl",
    devbranch="main",
)
