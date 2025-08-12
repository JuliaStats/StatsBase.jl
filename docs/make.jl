using Documenter, StatsBase, StatsAPI, Statistics, Random, LinearAlgebra

# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

DocMeta.setdocmeta!(StatsBase, :DocTestSetup, :(using StatsBase))

makedocs(
    sitename = "StatsBase.jl",
    modules = [StatsBase, StatsAPI],
    format = Documenter.HTML(assets = ["assets/favicon.ico"]),
    pages = ["index.md",
             "weights.md",
             "scalarstats.md",
             "robust.md",
             "deviation.md",
             "cov.md",
             "counts.md",
             "ranking.md",
             "sampling.md",
             "empirical.md",
             "signalcorr.md",
             "multivariate.md",
             "misc.md",
             "statmodels.md",
             "transformations.md"],
    checkdocs=:exports
)

deploydocs(
    repo = "github.com/JuliaStats/StatsBase.jl.git"
)
