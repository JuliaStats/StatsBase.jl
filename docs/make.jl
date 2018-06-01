using Documenter, StatsBase

makedocs(
    format = :html,
    sitename = "StatsBase.jl",
    modules = [StatsBase],
    pages = ["index.md",
             "weights.md",
             "means.md",
             "scalarstats.md",
             "robust.md",
             "deviation.md",
             "cov.md",
             "counts.md",
             "ranking.md",
             "sampling.md",
             "empirical.md",
             "signalcorr.md",
             "misc.md",
             "statmodels.md"]
)

deploydocs(
    repo = "github.com/JuliaStats/StatsBase.jl.git",
    target = "build",
    julia  = "0.7",
    deps   = nothing,
    make   = nothing
)
