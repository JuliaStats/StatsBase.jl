using StatsBase

tests = ["mathfuns", 
         "weights",
         "means", 
         "scalarstats", 
         "counts", 
         "ranking", 
         "sampling", 
         "empirical", 
         "hist", 
         "rankcorr",
         "signalcorr", 
         "misc"]

println("Running tests:")

for t in tests
    tfile = joinpath("test", "$(t).jl")
    println(" * $(tfile) ...")
    include(tfile)
end

