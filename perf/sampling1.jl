# Benchmark on non-weighted sampling 

# require the BenchmarkLite package
using BenchmarkLite

import StatsBase: direct_sample!, ordered_sample!

### generic sampling benchmarking

type SampleProc{Alg} <: Proc end

type DirectSample end
tsample!(s::DirectSample, a::AbstractArray, x::Array) = direct_sample!(a, x)

type OrderedSample end
tsample!(s::OrderedSample, a::AbstractArray, x::Array) = ordered_sample!(a, x)

Base.string(::SampleProc{DirectSample}) = "direct-sample"
Base.string(::SampleProc{OrderedSample}) = "ordered-sample"

# config is in the form of (n, k)

Base.length(p::SampleProc, cfg::(Int, Int)) = cfg[2]
Base.isvalid(p::SampleProc, cfg::(Int, Int)) = cfg[1] >= 1 && cfg[2] >= 1

Base.start(p::SampleProc, cfg::(Int, Int)) = Array(Int, cfg[2])
Base.run{Alg}(p::SampleProc{Alg}, cfg::(Int, Int), s::Vector{Int}) = tsample!(Alg(), 1:cfg[1], s)
Base.done(p::SampleProc, cfg, s) = nothing


### benchmarking

const ns = 5 * (2 .^ [0:9])
const ks = 2 .^ [1:16]

const procs = Proc[ SampleProc{DirectSample}(), 
                    SampleProc{OrderedSample}() ]

const cfgs = vec([(n, k) for k in ks, n in ns])

rtable = run(procs, cfgs; duration=0.3)
println()

show(rtable; unit=:mps)
println()
