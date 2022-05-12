using FastDMTransform
using Documenter

DocMeta.setdocmeta!(FastDMTransform, :DocTestSetup, :(using FastDMTransform); recursive=true)

makedocs(;
    modules=[FastDMTransform],
    authors="Kiran Shila <me@kiranshila.com> and contributors",
    repo="https://github.com/kiranshila/FastDMTransform.jl/blob/{commit}{path}#{line}",
    sitename="FastDMTransform.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kiranshila.github.io/FastDMTransform.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/kiranshila/FastDMTransform.jl",
    devbranch="main",
)
