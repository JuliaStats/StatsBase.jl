using StatsBase

tests = ["mathfuns", 
         "weights",
         "scalarstats", 
         "counts", 
         "ranking",  
         "empirical", 
         "hist", 
         "rankcorr",
         "signalcorr", 
         "misc", 
         "sampling"]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end

