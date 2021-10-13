# Sampling from Population

## Sampling API

The package provides functions for sampling from a given population (with or without replacement).

```@docs
sample
sample!
wsample
wsample!
```

## Algorithms

Internally, this package implements multiple algorithms, and the `sample` (and `sample!`)
methods integrate them into a poly-algorithm, which chooses a specific algorithm based
on inputs.

Note that the choices made in `sample` are decided based on extensive benchmarking
(see `perf/sampling.jl` and `perf/wsampling.jl`). It performs reasonably fast for most cases.
That being said, if you know that a certain algorithm is particularly suitable for your context,
directly calling an internal algorithm function might be slightly more efficient.

Here are a list of algorithms implemented in the package. The functions below are not exported
(one can still import them from StatsBase via `using` though).

### Notations

- `a`: source array representing the population
- `x`: the destination array
- `wv`: the weight vector (of type `AbstractWeights`), for weighted sampling
- `n`: the length of `a`
- `k`: the length of `x`. For sampling without replacement, `k` must not exceed `n`.
- `rng`: optional random number generator (defaults to `Random.GLOBAL_RNG`)

All following functions write results to `x` (pre-allocated) and return `x`.


### Sampling Algorithms (Non-Weighted)

```@docs
StatsBase.direct_sample!(rng::Random.AbstractRNG, a::AbstractArray, x::AbstractArray)
samplepair
StatsBase.knuths_sample!
StatsBase.fisher_yates_sample!
StatsBase.self_avoid_sample!
StatsBase.seqsample_a!
StatsBase.seqsample_c!
StatsBase.seqsample_d!
```

### Weighted Sampling Algorithms

```@docs
StatsBase.direct_sample!(rng::Random.AbstractRNG, a::AbstractArray, wv::AbstractWeights, x::AbstractArray)
StatsBase.alias_sample!
StatsBase.naive_wsample_norep!
StatsBase.efraimidis_a_wsample_norep!
StatsBase.efraimidis_ares_wsample_norep!
```
