using Documenter
using StochasticIntegrators

DocMeta.setdocmeta!(StochasticIntegrators, :DocTestSetup,
    :(using StochasticIntegrators); recursive = true)

makedocs(;
    modules = [StochasticIntegrators],
    warnonly = Documenter.except(:autodocs_block, :cross_references, :docs_block, :doctest,
        :eval_block, :example_block, :footnote, :linkcheck_remotes,
        :linkcheck, :meta_block, :parse_error, :setup_block),
    authors = "Michael Kraus, Tomasz M. Tyranowski",
    sitename = "StochasticIntegrators.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuliaGNI.github.io/StochasticIntegrators.jl",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Noise Processes" => "noise.md",
        "Methods" => "methods.md",
        "Implementation" => "implementation.md",
        "Index" => "reference.md"
    ]
)

deploydocs(;
    repo = "github.com/JuliaGNI/StochasticIntegrators.jl",
    devurl = "latest",
    devbranch = "main"
)
