using Compat, StatsBase

@static if StatsBase.BASESTATS_IN_STATSBASE
    using Random
    using LinearAlgebra
else
    using Compat.LinearAlgebra
    using Compat.Random
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
         "statmodels",
         "statscompat"]#,
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
