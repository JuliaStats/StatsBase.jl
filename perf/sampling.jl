# Benchmark on non-weighted sampling
# requires BenchmarkTools

using StatsBase
using BenchmarkTools
using Random

import StatsBase: direct_sample!, alias_sample!
import StatsBase: knuths_sample!, fisher_yates_sample!, self_avoid_sample!
import StatsBase: seqsample_a!, seqsample_c!

# set up a testing suite

suite = BenchmarkGroup()
suite["sampling"] = BenchmarkGroup(["sampling", "unweighted"])
suite["wsampling"] = BenchmarkGroup(["sampling", "weighted"])

# functions to choose params

function chooseParamsRep()
  n = 5 * (2 ^ rand(0:9))
  k = 2 ^ rand(1:16)
  return [1:n, zeros(Int, k)]
end

function chooseParamsNoRep()
  n = 5 * (2 ^ rand(0:9))
  k = rand(1:n)
  return [1:n, zeros(Int, k)]
end

function chooseParamsWRep()
  n = 5 * (2 ^ rand(0:9))
  k = 2 ^ rand(1:16)
  return [1:n, zeros(Int, k)]
end

# unweighted benchmarks

## No Replacement
suite["sampling"]["Direct"] = @benchmarkable  direct_sample!(p[1], p[2]) setup=(p=chooseParamsRep())
## With Replacement
suite["sampling"]["Knuths"] = @benchmarkable  knuths_sample!(p[1], p[2]) setup=(p=chooseParamsNoRep())
suite["sampling"]["Fisher_Yates"] = @benchmarkable  fisher_yates_sample!(p[1], p[2]) setup=(p=chooseParamsNoRep())
suite["sampling"]["Fisher_Yates"] = @benchmarkable  self_avoid_sample!(p[1], p[2]) setup=(p=chooseParamsNoRep())
suite["sampling"]["Seq_A"] = @benchmarkable  seqsample_a!(p[1], p[2]) setup=(p=chooseParamsNoRep())
suite["sampling"]["Seq_C"] = @benchmarkable  seqsample_c!(p[1], p[2]) setup=(p=chooseParamsNoRep())
suite["sampling"]["Sample_NoRep"] = @benchmarkable  sample!(p[1], p[2], replace=false, ordered=false) setup=(p=chooseParamsNoRep())
suite["sampling"]["Sample_NoRep_Ord"] = @benchmarkable  sample!(p[1], p[2], replace=false, ordered=true) setup=(p=chooseParamsNoRep())

# weighted benchmarks

## With Replacement
suite["wsampling"]["Direct"] = @benchmarkable  direct_sample!(Random.GLOBAL_RNG, p[1], p[2]) setup=(p=chooseParamsRep())
suite["wsampling"]["Alias"] = @benchmarkable  alias_sample!(Random.GLOBAL_RNG, p[1], p[2]) setup=(p=chooseParamsWRep())
suite["wsampling"]["Direct_S"] = @benchmarkable  sort!(direct_sample!(Random.GLOBAL_RNG, p[1], p[2])) setup=(p=chooseParamsWRep())
suite["wsampling"]["Sample_WRep"] = @benchmarkable  sample!(p[1], p[2], ordered=false) setup=(p=chooseParamsWRep())
suite["wsampling"]["Sample_WRep_Ord"] = @benchmarkable  sample!(p[1], p[2], ordered=true) setup=(p=chooseParamsWRep())

# run the tests
tune!(suite);
results = run(suite, verbose = true, seconds = 5)

# Show results
show(results)
