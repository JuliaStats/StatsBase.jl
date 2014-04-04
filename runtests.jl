using StatsBase

tests = ["means", 
         "scalarstats", 
         "counts", 
         "ranking", 
         "sampling", 
         "empirical", 
         "rankcorr",
         "signalcorr", 
         "canoncorr",
         "misc"]

println("Running tests:")

for t in tests
    tfile = joinpath("test", "$(t).jl")
    println(" * $(tfile) ...")
    include(tfile)
end

