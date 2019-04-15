using Documenter, StatsBase, Statistics, Random

# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

makedocs(
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
             "multivariate.md",
             "misc.md",
             "statmodels.md",
             "transformations.md"]
)

deploydocs(
    repo = "github.com/JuliaStats/StatsBase.jl.git"
)