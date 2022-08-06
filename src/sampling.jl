
###########################################################
#
#   (non-weighted) sampling
#
###########################################################

using Random: Sampler, Random.GLOBAL_RNG

### Algorithms for sampling with replacement

function direct_sample!(rng::AbstractRNG, a::UnitRange, x::AbstractArray)
    s = Sampler(rng, 1:length(a))
    b = a[1] - 1
    if b == 0
        for i = 1:length(x)
            @inbounds x[i] = rand(rng, s)
        end
    else
        for i = 1:length(x)
            @inbounds x[i] = b + rand(rng, s)
        end
    end
    return x
end
direct_sample!(a::UnitRange, x::AbstractArray) = direct_sample!(Random.GLOBAL_RNG, a, x)

"""
    direct_sample!([rng], a::AbstractArray, x::AbstractArray)

Direct sampling: for each `j` in `1:k`, randomly pick `i` from `1:n`,
and set `x[j] = a[i]`, with `n=length(a)` and `k=length(x)`.

This algorithm consumes `k` random numbers.
"""
function direct_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    s = Sampler(rng, 1:length(a))
    for i = 1:length(x)
        @inbounds x[i] = a[rand(rng, s)]
    end
    return x
end
direct_sample!(a::AbstractArray, x::AbstractArray) = direct_sample!(Random.GLOBAL_RNG, a, x)

# check whether we can use T to store indices 1:n exactly, and
# use some heuristics to decide whether it is beneficial for k samples
# (true for a subset of hardware-supported numeric types)
_storeindices(n, k, ::Type{T}) where {T<:Integer} = n ≤ typemax(T)
_storeindices(n, k, ::Type{T}) where {T<:Union{Float32,Float64}} = k < 22 && n ≤ maxintfloat(T)
_storeindices(n, k, ::Type{Complex{T}}) where {T} = _storeindices(n, k, T)
_storeindices(n, k, ::Type{Rational{T}}) where {T} = k < 16 && _storeindices(n, k, T)
_storeindices(n, k, T) = false
storeindices(n, k, ::Type{T}) where {T<:Base.HWNumber}  = _storeindices(n, k, T)
storeindices(n, k, T) = false

# order results of a sampler that does not order automatically
function sample_ordered!(sampler!, rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n, k = length(a), length(x)
    # todo: if eltype(x) <: Real && eltype(a) <: Real,
    #       in some cases it might be faster to check
    #       issorted(a) to see if we can just sort x
    if storeindices(n, k, eltype(x))
        sort!(sampler!(rng, Base.OneTo(n), x), by=real, lt=<)
        @inbounds for i = 1:k
            x[i] = a[Int(x[i])]
        end
    else
        indices = Array{Int}(undef, k)
        sort!(sampler!(rng, Base.OneTo(n), indices))
        @inbounds for i = 1:k
            x[i] = a[indices[i]]
        end
    end
    return x
end

# special case of a range can be done more efficiently
sample_ordered!(sampler!, rng::AbstractRNG, a::AbstractRange, x::AbstractArray) =
    sort!(sampler!(rng, a, x), rev=step(a)<0)

# weighted case:
sample_ordered!(sampler!, rng::AbstractRNG, a::AbstractArray,
                wv::AbstractWeights, x::AbstractArray) =
    sample_ordered!(rng, a, x) do rng, a, x
        sampler!(rng, a, wv, x)
    end

### draw a pair of distinct integers in [1:n]

"""
    samplepair([rng], n)

Draw a pair of distinct integers between 1 and `n` without replacement.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function samplepair(rng::AbstractRNG, n::Int)
    i1 = rand(rng, 1:n)
    i2 = rand(rng, 1:n-1)
    return (i1, ifelse(i2 == i1, n, i2))
end
samplepair(n::Int) = samplepair(Random.GLOBAL_RNG, n)

"""
    samplepair([rng], a)

Draw a pair of distinct elements from the array `a` without replacement.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function samplepair(rng::AbstractRNG, a::AbstractArray)
    i1, i2 = samplepair(rng, length(a))
    return a[i1], a[i2]
end
samplepair(a::AbstractArray) = samplepair(Random.GLOBAL_RNG, a)

### Algorithm for sampling without replacement

"""
    knuths_sample!([rng], a, x)

*Knuth's Algorithm S* for random sampling without replacement.

Reference: D. Knuth. *The Art of Computer Programming*. Vol 2, 3.4.2, p.142.

This algorithm consumes `length(a)` random numbers. It requires no additional
memory space. Suitable for the case where memory is tight.
"""
function knuths_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray;
                        initshuffle::Bool=true)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    # initialize
    for i = 1:k
        @inbounds x[i] = a[i]
    end
    if initshuffle
        @inbounds for j = 1:k
            l = rand(rng, j:k)
            if l != j
                t = x[j]
                x[j] = x[l]
                x[l] = t
            end
        end
    end

    # scan remaining
    s = Sampler(rng, 1:k)
    for i = k+1:n
        if rand(rng) * i < k  # keep it with probability k / i
            @inbounds x[rand(rng, s)] = a[i]
        end
    end
    return x
end
knuths_sample!(a::AbstractArray, x::AbstractArray; initshuffle::Bool=true) =
    knuths_sample!(Random.GLOBAL_RNG, a, x; initshuffle=initshuffle)

"""
    fisher_yates_sample!([rng], a::AbstractArray, x::AbstractArray)

Fisher-Yates shuffling (with early termination).

Pseudo-code:
```
n = length(a)
k = length(x)

# Create an array of the indices
inds = collect(1:n)

for i = 1:k
    # swap element `i` with another random element in inds[i:n]
    # set element `i` in `x`
end
```

This algorithm consumes `k=length(x)` random numbers. It uses an integer array of
length `n=length(a)` internally to maintain the shuffled indices. It is considerably
faster than Knuth's algorithm especially when `n` is greater than `k`.
It is ``O(n)`` for initialization, plus ``O(k)`` for random shuffling
"""
function fisher_yates_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    inds = Vector{Int}(undef, n)
    for i = 1:n
        @inbounds inds[i] = i
    end

    @inbounds for i = 1:k
        j = rand(rng, i:n)
        t = inds[j]
        inds[j] = inds[i]
        inds[i] = t
        x[i] = a[t]
    end
    return x
end
fisher_yates_sample!(a::AbstractArray, x::AbstractArray) =
    fisher_yates_sample!(Random.GLOBAL_RNG, a, x)

"""
    self_avoid_sample!([rng], a::AbstractArray, x::AbstractArray)

Self-avoid sampling: use a set to maintain the index that has been sampled.
Each time draw a new index, if the index has already been sampled,
redraw until it draws an unsampled one.

This algorithm consumes about (or slightly more than) `k=length(x)` random numbers,
and requires ``O(k)`` memory to store the set of sampled indices.
Very fast when ``n >> k``, with `n=length(a)`.

However, if `k` is large and approaches ``n``, the rejection rate would increase
drastically, resulting in poorer performance.
"""
function self_avoid_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    s = Set{Int}()
    sizehint!(s, k)
    rgen = Sampler(rng, 1:n)

    # first one
    idx = rand(rng, rgen)
    x[1] = a[idx]
    push!(s, idx)

    # remaining
    for i = 2:k
        idx = rand(rng, rgen)
        while idx in s
            idx = rand(rng, rgen)
        end
        x[i] = a[idx]
        push!(s, idx)
    end
    return x
end
self_avoid_sample!(a::AbstractArray, x::AbstractArray) =
    self_avoid_sample!(Random.GLOBAL_RNG, a, x)

"""
    seqsample_a!([rng], a::AbstractArray, x::AbstractArray)

Random subsequence sampling using algorithm A described in the following paper (page 714):
Jeffrey Scott Vitter. "Faster Methods for Random Sampling". Communications of the ACM,
27 (7), July 1984.

This algorithm consumes ``O(n)`` random numbers, with `n=length(a)`.
The outputs are ordered.
"""
function seqsample_a!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    i = 0
    j = 0
    while k > 1
        u = rand(rng)
        q = (n - k) / n
        while q > u  # skip
            i += 1
            n -= 1
            q *= (n - k) / n
        end
        @inbounds x[j+=1] = a[i+=1]
        n -= 1
        k -= 1
    end

    if k > 0  # checking k > 0 is necessary: x can be empty
        s = trunc(Int, n * rand(rng))
        x[j+1] = a[i+(s+1)]
    end
    return x
end
seqsample_a!(a::AbstractArray, x::AbstractArray) = seqsample_a!(Random.GLOBAL_RNG, a, x)

"""
    seqsample_c!([rng], a::AbstractArray, x::AbstractArray)

Random subsequence sampling using algorithm C described in the following paper (page 715):
Jeffrey Scott Vitter. "Faster Methods for Random Sampling". Communications of the ACM,
27 (7), July 1984.

This algorithm consumes ``O(k^2)`` random numbers, with `k=length(x)`.
The outputs are ordered.
"""
function seqsample_c!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    i = 0
    j = 0
    while k > 1
        l = n - k + 1
        minv = l
        u = n
        while u >= l
            v = u * rand(rng)
            if v < minv
                minv = v
            end
            u -= 1
        end
        s = trunc(Int, minv) + 1
        x[j+=1] = a[i+=s]
        n -= s
        k -= 1
    end

    if k > 0
        s = trunc(Int, n * rand(rng))
        x[j+1] = a[i+(s+1)]
    end
    return x
end
seqsample_c!(a::AbstractArray, x::AbstractArray) = seqsample_c!(Random.GLOBAL_RNG, a, x)

"""
    seqsample_d!([rng], a::AbstractArray, x::AbstractArray)

Random subsequence sampling using algorithm D described in the following paper (page 716-17):
Jeffrey Scott Vitter. "Faster Methods for Random Sampling". Communications of the ACM,
27 (7), July 1984.

This algorithm consumes ``O(k)`` random numbers, with `k=length(x)`.
The outputs are ordered.
"""
function seqsample_d!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    N = length(a)
    n = length(x)
    n <= N || error("length(x) should not exceed length(a)")

    i = 0
    j = 0

    vprime = exp(-randexp(rng)/n)
    q1 = N - n + 1
    q2 = q1 / N
    alpha = 1 / 13 # choose alpha value
    threshold = alpha * n

    while n > 1 && threshold < N
        while true
           local X
            while true
                X = N * (1 - vprime)
                s = trunc(Int, X)
                if s < q1
                    break
                end
                vprime = exp(-randexp(rng)/n)
            end

            y = rand(rng) / q2
            lhs = exp(log(y) / (n - 1))
            rhs = ((q1 - s) / q1) * (N / (N - X))

            if lhs <= rhs
                vprime = lhs / rhs
                break
            end

            if n - 1 > s
                bottom = N - n
                limit = N - s
            else
                bottom = N - s - 1
                limit = q1
            end

            top = N - 1

            while top >= limit
                y = y * top / bottom
                bottom -= 1
                top -= 1
            end

            if log(y) < (n - 1)*(log(N) - log(N - X))
                vprime = exp(-randexp(rng) / (n-1))
                break
            end
            vprime = exp(-randexp(rng)/n)
        end

        j += 1
        i += s+1
        @inbounds x[j] = a[i]
        N = N - s - 1
        n -= 1
        q1 -= s
        q2 = q1 / N
        threshold -= alpha
    end

    if n > 1
        seqsample_a!(rng, a[i+1:end], @view x[j+1:end])
    else
        s = trunc(Int, N * vprime)
        @inbounds x[j+=1] = a[i+=s+1]
    end
end

seqsample_d!(a::AbstractArray, x::AbstractArray) = seqsample_d!(Random.GLOBAL_RNG, a, x)


### Interface functions (poly-algorithms)
"""
    sample([rng], a, [wv::AbstractWeights])

Select a single random element of `a`. Sampling probabilities are proportional to
the weights given in `wv`, if provided.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
sample(rng::AbstractRNG, a::AbstractArray) = a[rand(rng, 1:length(a))]
sample(a::AbstractArray) = sample(Random.GLOBAL_RNG, a)


"""
    sample!([rng], a, [wv::AbstractWeights], x; replace=true, ordered=false)

Draw a random sample of `length(x)` elements from an array `a`
and store the result in `x`. A polyalgorithm is used for sampling.
Sampling probabilities are proportional to the weights given in `wv`,
if provided. `replace` dictates whether sampling is performed with
replacement. `ordered` dictates whether
an ordered sample (also called a sequential sample, i.e. a sample where
items appear in the same order as in `a`) should be taken.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).

Output array `a` must not be the same object as `x` or `wv`
nor share memory with them, or the result may be incorrect.
"""
function sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray;
                 replace::Bool=true, ordered::Bool=false)
    Base.mightalias(a, x) &&
        throw(ArgumentError("output array a must not share memory with input array x"))
    1 == firstindex(a) == firstindex(x) ||
        throw(ArgumentError("non 1-based arrays are not supported"))
    n = length(a)
    k = length(x)
    k == 0 && return x

    if replace  # with replacement
        if ordered
            sample_ordered!(direct_sample!, rng, a, x)
        else
            direct_sample!(rng, a, x)
        end

    else  # without replacement
        k <= n || error("Cannot draw more samples without replacement.")

        if ordered
            if n > 10 * k * k
                seqsample_c!(rng, a, x)
            else
                seqsample_a!(rng, a, x)
            end
        else
            if k == 1
                @inbounds x[1] = sample(rng, a)
            elseif k == 2
                @inbounds (x[1], x[2]) = samplepair(rng, a)
            elseif n < k * 24
                fisher_yates_sample!(rng, a, x)
            else
                self_avoid_sample!(rng, a, x)
            end
        end
    end
    return x
end
sample!(a::AbstractArray, x::AbstractArray; replace::Bool=true, ordered::Bool=false) =
    sample!(Random.GLOBAL_RNG, a, x; replace=replace, ordered=ordered)


"""
    sample([rng], a, [wv::AbstractWeights], n::Integer; replace=true, ordered=false)

Select a random, optionally weighted sample of size `n` from an array `a`
using a polyalgorithm. Sampling probabilities are proportional to the weights
given in `wv`, if provided. `replace` dictates whether sampling is performed
with replacement. `ordered` dictates whether
an ordered sample (also called a sequential sample, i.e. a sample where
items appear in the same order as in `a`) should be taken.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function sample(rng::AbstractRNG, a::AbstractArray{T}, n::Integer;
                replace::Bool=true, ordered::Bool=false) where T
    sample!(rng, a, Vector{T}(undef, n); replace=replace, ordered=ordered)
end
sample(a::AbstractArray, n::Integer; replace::Bool=true, ordered::Bool=false) =
    sample(Random.GLOBAL_RNG, a, n; replace=replace, ordered=ordered)


"""
    sample([rng], a, [wv::AbstractWeights], dims::Dims; replace=true, ordered=false)

Select a random, optionally weighted sample from an array `a` specifying
the dimensions `dims` of the output array. Sampling probabilities are
proportional to the weights given in `wv`, if provided. `replace` dictates
whether sampling is performed with replacement. `ordered` dictates whether
an ordered sample (also called a sequential sample, i.e. a sample where
items appear in the same order as in `a`) should be taken.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function sample(rng::AbstractRNG, a::AbstractArray{T}, dims::Dims;
                replace::Bool=true, ordered::Bool=false) where T
    sample!(rng, a, Array{T}(undef, dims); replace=replace, ordered=ordered)
end
sample(a::AbstractArray, dims::Dims; replace::Bool=true, ordered::Bool=false) =
    sample(Random.GLOBAL_RNG, a, dims; replace=replace, ordered=ordered)

################################################################
#
#  Weighted sampling
#
################################################################

"""
    sample([rng], wv::AbstractWeights)

Select a single random integer in `1:length(wv)` with probabilities
proportional to the weights given in `wv`.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function sample(rng::AbstractRNG, wv::AbstractWeights)
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = 1
    cw = wv[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    return i
end
sample(wv::AbstractWeights) = sample(Random.GLOBAL_RNG, wv)

sample(rng::AbstractRNG, a::AbstractArray, wv::AbstractWeights) = a[sample(rng, wv)]
sample(a::AbstractArray, wv::AbstractWeights) = sample(Random.GLOBAL_RNG, a, wv)

"""
    direct_sample!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)

Direct sampling.

Draw each sample by scanning the weight vector.

Noting `k=length(x)` and `n=length(a)`, this algorithm:
* consumes `k` random numbers
* has time complexity ``O(n k)``, as scanning the weight vector each time takes ``O(n)``
* requires no additional memory space.
"""
function direct_sample!(rng::AbstractRNG, a::AbstractArray,
                        wv::AbstractWeights, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("Inconsistent lengths."))
    for i = 1:length(x)
        x[i] = a[sample(rng, wv)]
    end
    return x
end
direct_sample!(a::AbstractArray, wv::AbstractWeights, x::AbstractArray) =
    direct_sample!(Random.GLOBAL_RNG, a, wv, x)

function make_alias_table!(w::AbstractVector, wsum,
                           a::AbstractVector{Float64},
                           alias::AbstractVector{Int})
    # Arguments:
    #
    #   w [in]:         input weights
    #   wsum [in]:      pre-computed sum(w)
    #
    #   a [out]:        acceptance probabilities
    #   alias [out]:    alias table
    #
    # Note: a and w can be the same array, then that array will be
    #       overwritten inplace by acceptance probabilities
    #
    # Returns nothing
    #

    n = length(w)
    length(a) == length(alias) == n ||
        throw(DimensionMismatch("Inconsistent array lengths."))

    ac = n / wsum
    for i = 1:n
        @inbounds a[i] = w[i] * ac
    end

    larges = Vector{Int}(undef, n)
    smalls = Vector{Int}(undef, n)
    kl = 0  # actual number of larges
    ks = 0  # actual number of smalls

    for i = 1:n
        @inbounds ai = a[i]
        if ai > 1.0
            larges[kl+=1] = i  # push to larges
        elseif ai < 1.0
            smalls[ks+=1] = i  # push to smalls
        end
    end

    while kl > 0 && ks > 0
        s = smalls[ks]; ks -= 1  # pop from smalls
        l = larges[kl]; kl -= 1  # pop from larges
        @inbounds alias[s] = l
        @inbounds al = a[l] = (a[l] - 1.0) + a[s]
        if al > 1.0
            larges[kl+=1] = l  # push to larges
        else
            smalls[ks+=1] = l  # push to smalls
        end
    end

    # this loop should be redundant, except for rounding
    for i = 1:ks
        @inbounds a[smalls[i]] = 1.0
    end
    nothing
end

"""
    alias_sample!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)

Alias method.

Build an alias table, and sample therefrom.

Reference: Walker, A. J. "An Efficient Method for Generating Discrete Random Variables
with General Distributions." *ACM Transactions on Mathematical Software* 3 (3): 253, 1977.

Noting `k=length(x)` and `n=length(a)`, this algorithm takes ``O(n \\log n)`` time
for building the alias table, and then ``O(1)`` to draw each sample. It consumes ``2 k`` random numbers.
"""
function alias_sample!(rng::AbstractRNG, a::AbstractArray, wv::AbstractWeights, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("Inconsistent lengths."))

    # create alias table
    ap = Vector{Float64}(undef, n)
    alias = Vector{Int}(undef, n)
    make_alias_table!(wv, sum(wv), ap, alias)

    # sampling
    s = Sampler(rng, 1:n)
    for i = 1:length(x)
        j = rand(rng, s)
        x[i] = rand(rng) < ap[j] ? a[j] : a[alias[j]]
    end
    return x
end
alias_sample!(a::AbstractArray, wv::AbstractWeights, x::AbstractArray) =
    alias_sample!(Random.GLOBAL_RNG, a, wv, x)

"""
    naive_wsample_norep!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)

Naive implementation of weighted sampling without replacement.

It makes a copy of the weight vector at initialization, and sets the weight to zero
when the corresponding sample is picked.

Noting `k=length(x)` and `n=length(a)`, this algorithm consumes ``O(k)`` random numbers,
and has overall time complexity ``O(n k)``.
"""
function naive_wsample_norep!(rng::AbstractRNG, a::AbstractArray,
                              wv::AbstractWeights, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("Inconsistent lengths."))
    k = length(x)

    w = Vector{Float64}(undef, n)
    copyto!(w, wv)
    wsum = sum(wv)

    for i = 1:k
        u = rand(rng) * wsum
        j = 1
        c = w[1]
        while c < u && j < n
            @inbounds c += w[j+=1]
        end
        @inbounds x[i] = a[j]

        @inbounds wsum -= w[j]
        @inbounds w[j] = 0.0
    end
    return x
end
naive_wsample_norep!(a::AbstractArray, wv::AbstractWeights, x::AbstractArray) =
    naive_wsample_norep!(Random.GLOBAL_RNG, a, wv, x)

# Weighted sampling without replacement
# Instead of keys u^(1/w) where u = random(0,1) keys w/v where v = randexp(1) are used.
"""
    efraimidis_a_wsample_norep!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)

Weighted sampling without replacement using Efraimidis-Spirakis A algorithm.

Reference: Efraimidis, P. S., Spirakis, P. G. "Weighted random sampling with a reservoir."
*Information Processing Letters*, 97 (5), 181-185, 2006. doi:10.1016/j.ipl.2005.11.003.

Noting `k=length(x)` and `n=length(a)`, this algorithm takes ``O(n + k \\log k)``
processing time to draw ``k`` elements. It consumes ``n`` random numbers.
"""
function efraimidis_a_wsample_norep!(rng::AbstractRNG, a::AbstractArray,
                                     wv::AbstractWeights, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("a and wv must be of same length (got $n and $(length(wv)))."))
    k = length(x)

    # calculate keys for all items
    keys = randexp(rng, n)
    for i in 1:n
        @inbounds keys[i] = wv.values[i]/keys[i]
    end

    # return items with largest keys
    index = sortperm(keys; alg = PartialQuickSort(k), rev = true)
    for i in 1:k
        @inbounds x[i] = a[index[i]]
    end
    return x
end
efraimidis_a_wsample_norep!(a::AbstractArray, wv::AbstractWeights, x::AbstractArray) =
    efraimidis_a_wsample_norep!(Random.GLOBAL_RNG, a, wv, x)

# Weighted sampling without replacement
# Instead of keys u^(1/w) where u = random(0,1) keys w/v where v = randexp(1) are used.
"""
    efraimidis_ares_wsample_norep!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)

Implementation of weighted sampling without replacement using Efraimidis-Spirakis A-Res algorithm.

Reference: Efraimidis, P. S., Spirakis, P. G. "Weighted random sampling with a reservoir."
*Information Processing Letters*, 97 (5), 181-185, 2006. doi:10.1016/j.ipl.2005.11.003.

Noting `k=length(x)` and `n=length(a)`, this algorithm takes ``O(k \\log(k) \\log(n / k))``
processing time to draw ``k`` elements. It consumes ``n`` random numbers.
"""
function efraimidis_ares_wsample_norep!(rng::AbstractRNG, a::AbstractArray,
                                        wv::AbstractWeights, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("a and wv must be of same length (got $n and $(length(wv)))."))
    k = length(x)
    k > 0 || return x

    # initialize priority queue
    pq = Vector{Pair{Float64,Int}}(undef, k)
    i = 0
    s = 0
    @inbounds for _s in 1:n
        s = _s
        w = wv.values[s]
        w < 0 && error("Negative weight found in weight vector at index $s")
        if w > 0
            i += 1
            pq[i] = (w/randexp(rng) => s)
        end
        i >= k && break
    end
    i < k && throw(DimensionMismatch("wv must have at least $k strictly positive entries (got $i)"))
    heapify!(pq)

    # set threshold
    @inbounds threshold = pq[1].first

    @inbounds for i in s+1:n
        w = wv.values[i]
        w < 0 && error("Negative weight found in weight vector at index $i")
        w > 0 || continue
        key = w/randexp(rng)

        # if key is larger than the threshold
        if key > threshold
            # update priority queue
            pq[1] = (key => i)
            percolate_down!(pq, 1)

            # update threshold
            threshold = pq[1].first
        end
    end

    # fill output array with items in descending order
    @inbounds for i in k:-1:1
        x[i] = a[heappop!(pq).second]
    end
    return x
end
efraimidis_ares_wsample_norep!(a::AbstractArray, wv::AbstractWeights, x::AbstractArray) =
    efraimidis_ares_wsample_norep!(Random.GLOBAL_RNG, a, wv, x)

# Weighted sampling without replacement
# Instead of keys u^(1/w) where u = random(0,1) keys w/v where v = randexp(1) are used.
"""
    efraimidis_aexpj_wsample_norep!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)

Implementation of weighted sampling without replacement using Efraimidis-Spirakis A-ExpJ algorithm.

Reference: Efraimidis, P. S., Spirakis, P. G. "Weighted random sampling with a reservoir."
*Information Processing Letters*, 97 (5), 181-185, 2006. doi:10.1016/j.ipl.2005.11.003.

Noting `k=length(x)` and `n=length(a)`, this algorithm takes ``O(k \\log(k) \\log(n / k))``
processing time to draw ``k`` elements. It consumes ``O(k \\log(n / k))`` random numbers.
"""
function efraimidis_aexpj_wsample_norep!(rng::AbstractRNG, a::AbstractArray,
                                         wv::AbstractWeights, x::AbstractArray;
                                         ordered::Bool=false)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("a and wv must be of same length (got $n and $(length(wv)))."))
    k = length(x)
    k > 0 || return x

    # initialize priority queue
    pq = Vector{Pair{Float64,Int}}(undef, k)
    i = 0
    s = 0
    @inbounds for _s in 1:n
        s = _s
        w = wv.values[s]
        w < 0 && error("Negative weight found in weight vector at index $s")
        if w > 0
            i += 1
            pq[i] = (w/randexp(rng) => s)
        end
        i >= k && break
    end
    i < k && throw(DimensionMismatch("wv must have at least $k strictly positive entries (got $i)"))
    heapify!(pq)

    # set threshold
    @inbounds threshold = pq[1].first
    X = threshold*randexp(rng)

    @inbounds for i in s+1:n
        w = wv.values[i]
        w < 0 && error("Negative weight found in weight vector at index $i")
        w > 0 || continue
        X -= w
        X <= 0 || continue

        # update priority queue
        t = exp(-w/threshold)
        pq[1] = (-w/log(t+rand(rng)*(1-t)) => i)
        percolate_down!(pq, 1)

        # update threshold
        threshold = pq[1].first
        X = threshold * randexp(rng)
    end
    if ordered
        # fill output array with items sorted as in a
        sort!(pq, by=last)
        @inbounds for i in 1:k
            x[i] = a[pq[i].second]
        end
    else
        # fill output array with items in descending order
        @inbounds for i in k:-1:1
            x[i] = a[heappop!(pq).second]
        end
    end
    return x
end
efraimidis_aexpj_wsample_norep!(a::AbstractArray, wv::AbstractWeights, x::AbstractArray;
                                ordered::Bool=false) =
    efraimidis_aexpj_wsample_norep!(Random.GLOBAL_RNG, a, wv, x; ordered=ordered)

function sample!(rng::AbstractRNG, a::AbstractArray, wv::AbstractWeights, x::AbstractArray;
                 replace::Bool=true, ordered::Bool=false)
    Base.mightalias(a, x) &&
        throw(ArgumentError("output array a must not share memory with input array x"))
    Base.mightalias(a, wv) &&
        throw(ArgumentError("output array a must not share memory with weights array wv"))
    1 == firstindex(a) == firstindex(wv) == firstindex(x) ||
        throw(ArgumentError("non 1-based arrays are not supported"))
    n = length(a)
    k = length(x)

    if replace
        if ordered
            sample_ordered!(rng, a, wv, x) do rng, a, wv, x
                sample!(rng, a, wv, x; replace=true, ordered=false)
            end
        else
            if n < 40
                direct_sample!(rng, a, wv, x)
            else
                t = ifelse(n < 500, 64, 32)
                if k < t
                    direct_sample!(rng, a, wv, x)
                else
                    alias_sample!(rng, a, wv, x)
                end
            end
        end
    else
        k <= n || error("Cannot draw $k samples from $n samples without replacement.")
        efraimidis_aexpj_wsample_norep!(rng, a, wv, x; ordered=ordered)
    end
    return x
end
sample!(a::AbstractArray, wv::AbstractWeights, x::AbstractArray;
        replace::Bool=true, ordered::Bool=false) =
    sample!(Random.GLOBAL_RNG, a, wv, x; replace=replace, ordered=ordered)

sample(rng::AbstractRNG, a::AbstractArray{T}, wv::AbstractWeights, n::Integer;
       replace::Bool=true, ordered::Bool=false) where {T} =
    sample!(rng, a, wv, Vector{T}(undef, n); replace=replace, ordered=ordered)
sample(a::AbstractArray, wv::AbstractWeights, n::Integer;
       replace::Bool=true, ordered::Bool=false) =
    sample(Random.GLOBAL_RNG, a, wv, n; replace=replace, ordered=ordered)

sample(rng::AbstractRNG, a::AbstractArray{T}, wv::AbstractWeights, dims::Dims;
       replace::Bool=true, ordered::Bool=false) where {T} =
    sample!(rng, a, wv, Array{T}(undef, dims); replace=replace, ordered=ordered)
sample(a::AbstractArray, wv::AbstractWeights, dims::Dims;
       replace::Bool=true, ordered::Bool=false) =
    sample(Random.GLOBAL_RNG, a, wv, dims; replace=replace, ordered=ordered)

# wsample interface

"""
    wsample!([rng], a, w, x; replace=true, ordered=false)

Select a weighted sample from an array `a` and store the result in `x`. Sampling
probabilities are proportional to the weights given in `w`. `replace` dictates
whether sampling is performed with replacement. `ordered` dictates whether
an ordered sample (also called a sequential sample, i.e. a sample where
items appear in the same order as in `a`) should be taken.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
wsample!(rng::AbstractRNG, a::AbstractArray, w::RealVector, x::AbstractArray;
         replace::Bool=true, ordered::Bool=false) =
    sample!(rng, a, weights(w), x; replace=replace, ordered=ordered)
wsample!(a::AbstractArray, w::RealVector, x::AbstractArray;
         replace::Bool=true, ordered::Bool=false) =
    sample!(Random.GLOBAL_RNG, a, weights(w), x; replace=replace, ordered=ordered)

"""
    wsample([rng], [a], w)

Select a weighted random sample of size 1 from `a` with probabilities proportional
to the weights given in `w`. If `a` is not present, select a random weight from `w`.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
wsample(rng::AbstractRNG, w::RealVector) = sample(rng, weights(w))
wsample(w::RealVector) = wsample(Random.GLOBAL_RNG, w)
wsample(rng::AbstractRNG, a::AbstractArray, w::RealVector) = sample(rng, a, weights(w))
wsample(a::AbstractArray, w::RealVector) = wsample(Random.GLOBAL_RNG, a, w)


"""
    wsample([rng], [a], w, n::Integer; replace=true, ordered=false)

Select a weighted random sample of size `n` from `a` with probabilities proportional
to the weights given in `w` if `a` is present, otherwise select a random sample of size
`n` of the weights given in `w`. `replace` dictates whether sampling is performed with
replacement. `ordered` dictates whether
an ordered sample (also called a sequential sample, i.e. a sample where
items appear in the same order as in `a`) should be taken.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
wsample(rng::AbstractRNG, a::AbstractArray{T}, w::RealVector, n::Integer;
        replace::Bool=true, ordered::Bool=false) where {T} =
    wsample!(rng, a, w, Vector{T}(undef, n); replace=replace, ordered=ordered)
wsample(a::AbstractArray, w::RealVector, n::Integer;
        replace::Bool=true, ordered::Bool=false) =
    wsample(Random.GLOBAL_RNG, a, w, n; replace=replace, ordered=ordered)

"""
    wsample([rng], [a], w, dims::Dims; replace=true, ordered=false)

Select a weighted random sample from `a` with probabilities proportional to the
weights given in `w` if `a` is present, otherwise select a random sample of size
`n` of the weights given in `w`. The dimensions of the output are given by `dims`.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
wsample(rng::AbstractRNG, a::AbstractArray{T}, w::RealVector, dims::Dims;
        replace::Bool=true, ordered::Bool=false) where {T} =
    wsample!(rng, a, w, Array{T}(undef, dims); replace=replace, ordered=ordered)
wsample(a::AbstractArray, w::RealVector, dims::Dims;
        replace::Bool=true, ordered::Bool=false) =
    wsample(Random.GLOBAL_RNG, a, w, dims; replace=replace, ordered=ordered)
