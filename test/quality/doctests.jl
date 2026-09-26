using Documenter
using StochasticIntegrators

@eval Main import StochasticIntegrators

DocMeta.setdocmeta!(StochasticIntegrators, :DocTestSetup,
    :(using StochasticIntegrators); recursive = true)

doctest(StochasticIntegrators)
