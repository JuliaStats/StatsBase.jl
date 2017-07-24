# Benchmark on weighted sampling

using BenchmarkLite
using StatsBase

import StatsBase: direct_sample!, alias_sample!, xmultinom_sample!

### procedure definition

mutable struct WSampleProc{Alg} <: Proc end

abstract type WithRep end
abstract type NoRep end

mutable struct Direct <: WithRep end
tsample!(s::Direct, wv, x) = direct_sample!(1:length(wv), wv, x)

mutable struct Alias <: WithRep end
tsample!(s::Alias, wv, x) = alias_sample!(1:length(wv), wv, x)

mutable struct Xmultinom_S <: WithRep end
tsample!(s::Xmultinom_S, wv, x) = shuffle!(xmultinom_sample!(1:length(wv), wv, x))

mutable struct Xmultinom <: WithRep end
tsample!(s::Xmultinom, wv, x) = xmultinom_sample!(1:length(wv), wv, x)

mutable struct Direct_S <: WithRep end
tsample!(s::Direct_S, wv, x) = sort!(direct_sample!(1:length(wv), wv, x))

mutable struct Sample_WRep <: WithRep end
tsample!(s::Sample_WRep, wv, x) = sample!(1:length(wv), wv, x; ordered=false)

mutable struct Sample_WRep_Ord <: WithRep end
tsample!(s::Sample_WRep_Ord, wv, x) = sample!(1:length(wv), wv, x; ordered=true)


# config is in the form of (n, k)

Base.string(p::WSampleProc{Alg}) where {Alg} = lowercase(string(Alg))

Base.length(p::WSampleProc, cfg::Tuple{Int,Int}) = cfg[2]
Base.isvalid(p::WSampleProc{<:WithRep}, cfg::Tuple{Int,Int}) = ((n, k) = cfg; n >= 1 && k >= 1)
Base.isvalid(p::WSampleProc{<:NoRep}, cfg::Tuple{Int,Int}) = ((n, k) = cfg; n >= k >= 1)

function Base.start(p::WSampleProc, cfg::Tuple{Int,Int})
    n, k = cfg
    x = Vector{Int}(k)
    w = weights(fill(1.0/n, n))
    return (w, x)
end

Base.run(p::WSampleProc{Alg}, cfg::Tuple{Int,Int}, s) where {Alg} = tsample!(Alg(), s[1], s[2])
Base.done(p::WSampleProc, cfg, s) = nothing


### benchmarking

const ns = 5 * (2 .^ [0:9])
const ks = 2 .^ [1:16]

## with replacement

const procs1 = Proc[ WSampleProc{Direct}(),
                     WSampleProc{Alias}(),
                     WSampleProc{Xmultinom_S}(),
                     WSampleProc{Sample_WRep}(),
                     WSampleProc{Xmultinom}(),
                     WSampleProc{Direct_S}(),
                     WSampleProc{Sample_WRep_Ord}() ]

const cfgs1 = vec([(n, k) for k in ks, n in ns])

rtable1 = run(procs1, cfgs1; duration=0.2)
println()


## show results

println("Sampling With Replacement")
println("===================================")
show(rtable1; unit=:mps, cfghead="(n, k)")
println()


