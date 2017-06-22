using Documenter, StatsBase

makedocs(
    format = :html,
    sitename = "StatsBase.jl",
    # options
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
    # options
    repo = "github.com/JuliaStats/StatsBase.jl.git",
    target = "build",
    julia  = "nightly",
    osname = "linux",
    deps = nothing,
    make = nothing
)
