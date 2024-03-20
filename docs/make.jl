using thesis_code
using Documenter

DocMeta.setdocmeta!(thesis_code, :DocTestSetup, :(using thesis_code); recursive=true)

makedocs(;
    modules=[thesis_code],
    authors="Bo Jacobs Strom <bojs.13@gmail.com> and contributors",
    sitename="thesis_code.jl",
    format=Documenter.HTML(;
        canonical="https://bo-js.github.io/thesis_code.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bo-js/thesis_code.jl",
    devbranch="main",
)
