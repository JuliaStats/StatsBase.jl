using Compat, StatsBase

@static if StatsBase.BASESTATS_IN_STATSBASE
    using Random
    using LinearAlgebra
    const Compatvar = var
    const Compatstd = std
    const Compatcov = cov
    const Compatcor = cor
else
    using Compat.LinearAlgebra
    using Compat.Random
    const Compatvar = Compat.var
    const Compatstd = Compat.std
    const Compatcov = Compat.cov
    const Compatcor = Compat.cor
end


tests = ["ambiguous",
         "weights",
         "moments",
         "scalarstats",
         "deviation",
         "cov",
         "counts",
         "ranking",
         "empirical",
         "hist",
         "rankcorr",
         "signalcorr",
         "misc",
         "robust",
         "sampling",
         "wsampling",
         "statmodels"]#,
         #"statquiz"]

if StatsBase.BASESTATS_IN_STATSBASE
    push!(tests, "base")
    include("dimensionful.jl")
end

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end
