using StatsBase

tests = ["mathfuns", 
         "weights",
         "moments",
         "scalarstats", 
         "cov",
         "counts", 
         "ranking",  
         "empirical", 
         "hist", 
         "rankcorr",
         "signalcorr", 
         "misc", 
         "sampling",
         "wsampling"]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end

