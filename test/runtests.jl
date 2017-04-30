using StatsBase

opts = Base.JLOptions()
depwarns = isdefined(opts, :depwarn) ? opts.depwarn == "no" : true
test_deprecates = if haskey(ENV, "TEST_DEPRECATES")
    lowercase(ENV["TEST_DEPRECATES"]) == "true"
else
    false
end

tests = ["weights",
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
         "sampling",
         "wsampling",
         "statmodels"]#,
         #"statquiz"]

if !depwarns || test_deprecates
    push!(tests, "deprecates")
end

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end
