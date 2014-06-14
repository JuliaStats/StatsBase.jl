# Benchmark on non-weighted sampling 

import StatsBase: direct_sample!, ordered_sample!

### generic sampling benchmarking

abstract SampleAlg

function benchmark_sampling1(s::SampleAlg, n::Int, k::Int, T::Int)
    a = [1:n]
    x = Array(Int, k)
    # warming
    for t=1:2; tsample!(s, a, x); end
    # benchmarking
    gc_disable()
    et = @elapsed for t=1:T
        tsample!(s, a, x)
    end 
    gc_enable()
    return float64(k) * T / (et * 1.0e6)  # seconds --> MPS
end

function benchmark_sampling(s::SampleAlg,               # the sampling algorithm
                            ns::AbstractVector{Int},    # a list of n-values to test
                            ks::AbstractVector{Int};    # a list of k-values to test
                            maxN::Int=2^18,             # the maximum number of total samples for each setting
                            verbose::Bool=true)         # whether to show progress information

    R = zeros(length(ns), length(ks))
    verbose && println("Benchmark sampling algorithm: $(typeof(s))")
    for (i, n) in enumerate(ns), (j, k) in enumerate(ks)
        T = div(maxN, k)
        mps = benchmark_sampling1(s, n, k, T)
        verbose && @printf("n = %5d, k = %5d: %10.4f MPS\n", n, k, mps)
        R[i,j] = mps
    end
    return R
end

function writeresults(io::IO, algname::String, 
                      ns::AbstractVector{Int}, 
                      ks::AbstractVector{Int},
                      dat::AbstractMatrix{Float64})

    for (i, n) in enumerate(ns), (j, k) in enumerate(ks)
        println(io, "$algname, $n, $k, $(dat[i,j])")
    end
end

### sampling with replacement

type Direct <: SampleAlg end
tsample!(s::Direct, a::Array, x::Array) = direct_sample!(a, x)

type Ordered <: SampleAlg end
tsample!(s::Ordered, a::Array, x::Array) = ordered_sample!(a, x)

const ns = 5 * (2 .^ [0:9])
const ks = 2 .^ [1:16]

direct_mps = benchmark_sampling(Direct(), ns, ks)
ordered_mps = benchmark_sampling(Ordered(), ns, ks)

const sample_wrep_path = joinpath(dirname(@__FILE__), "sample_wrep.csv")

open(sample_wrep_path, "w") do f
    println(f, "Algorithm, n, k, mps")
    writeresults(f, "direct", ns, ks, direct_mps)
    writeresults(f, "ordered", ns, ks, ordered_mps)
end
println("Results written to $(sample_wrep_path)")


