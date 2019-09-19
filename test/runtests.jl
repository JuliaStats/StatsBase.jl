using StatsBase
using Dates
using LinearAlgebra
using Random
using Statistics

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
         "partialcor",
         "transformations"]
         #"statquiz"]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end
