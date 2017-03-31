
###########################################################
#
#   (non-weighted) sampling
#
###########################################################

### Algorithms for sampling with replacement

# Direct sampling
#
#   for each i, draw an index from 1:length(a), and
#   pick the corresponding element from a.
#
#   x[i] and x[j] are independent, when i != j
#
function direct_sample!(rng::AbstractRNG, a::UnitRange, x::AbstractArray)
    s = RandIntSampler(length(a))
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
direct_sample!(a::UnitRange, x::AbstractArray) = direct_sample!(Base.GLOBAL_RNG, a, x)

function direct_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    s = RandIntSampler(length(a))
    for i = 1:length(x)
        @inbounds x[i] = a[rand(rng, s)]
    end
    return x
end
direct_sample!(a::AbstractArray, x::AbstractArray) = direct_sample!(Base.GLOBAL_RNG, a, x)

### draw a pair of distinct integers in [1:n]

"""
    samplepair([rng], n)

Draw a pair of distinct integers between 1 and `n` without replacement.
"""
function samplepair(rng::AbstractRNG, n::Int)
    i1 = randi(rng, n)
    i2 = randi(rng, n-1)
    return (i1, ifelse(i2 == i1, n, i2))
end
samplepair(n::Int) = samplepair(Base.GLOBAL_RNG, n)

"""
    samplepair([rng], a)

Draw a pair of distinct elements from the array `a` without replacement.
"""
function samplepair(rng::AbstractRNG, a::AbstractArray)
    i1, i2 = samplepair(rng, length(a))
    return a[i1], a[i2]
end
samplepair(a::AbstractArray) = samplepair(Base.GLOBAL_RNG, a)

### Algorithm for sampling without replacement

# Knuth's Algorithm S
#
#   Reference: The Art of Computer Programming, Vol 2, 3.4.2, p.142
#
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
            l = randi(rng, j, k)
            if l != j
                t = x[j]
                x[j] = x[l]
                x[l] = t
            end
        end
    end

    # scan remaining
    s = RandIntSampler(k)
    for i = k+1:n
        if rand(rng) * i < k  # keep it with probability k / i
            @inbounds x[rand(rng, s)] = a[i]
        end
    end
    return x
end
knuths_sample!(a::AbstractArray, x::AbstractArray; initshuffle::Bool=true) =
    knuths_sample!(Base.GLOBAL_RNG, a, x; initshuffle=initshuffle)

# Fisher-Yates sampling
#
#   create an array of index inds = [1:n]
#
#   for i = 1:k
#       swap inds[i] with a random one in inds[i:n]
#       set x[i] = a[inds[i]]
#   end
#
#   O(n) for initialization + O(k) for random shuffling
#
function fisher_yates_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    inds = Vector{Int}(n)
    for i = 1:n
        @inbounds inds[i] = i
    end

    @inbounds for i = 1:k
        j = randi(rng, i, n)
        t = inds[j]
        inds[j] = inds[i]
        inds[i] = t
        x[i] = a[t]
    end
    return x
end
fisher_yates_sample!(a::AbstractArray, x::AbstractArray) = fisher_yates_sample!(Base.GLOBAL_RNG, a, x)

# Self-avoid sampling
#
#   Use a set to maintain the index that has been sampled,
#
#   each time draw a new index, if the index has already
#   been sampled, redraw until it draws an unsampled one.
#
function self_avoid_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    s = Set{Int}()
    sizehint!(s, k)
    rgen = RandIntSampler(n)

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
self_avoid_sample!(a::AbstractArray, x::AbstractArray) = self_avoid_sample!(Base.GLOBAL_RNG, a, x)

# Random subsequence sampling
#
#  References:
#
#   Jeffrey Scott Vitter.
#   "Faster Methods for Random Sampling"
#   Communications of the ACM, 27 (7), July 1984
#
#
# This paper presents three algorithms, respectively
# named Algorithm A, C, and D.
#
# These algorithms are for sequential sampling, and the
# outputs are ordered. They are implemented below
#

## Algorithm A (page 714)
#
#  Require O(n) random numbers.
#
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
seqsample_a!(a::AbstractArray, x::AbstractArray) = seqsample_a!(Base.GLOBAL_RNG, a, x)

## Algorithm C (page 715)
#
#  Require O(k^2) random numbers
#
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
seqsample_c!(a::AbstractArray, x::AbstractArray) = seqsample_c!(Base.GLOBAL_RNG, a, x)

## TODO: implement Algorithm D (page 716 - 717)


### Interface functions (poly-algorithms)
"""
    sample([rng], a, [wv::WeightVec])

Select a single random element of `a`. Sampling probabilities are proportional to
the weights given in `wv`, if provided.
"""
sample(rng::AbstractRNG, a::AbstractArray) = a[randi(rng, length(a))]
sample(a::AbstractArray) = sample(Base.GLOBAL_RNG, a)


"""
    sample!([rng], a, [wv::WeightVec], x; replace=true, ordered=false)

Draw a random sample of `length(x)` elements from an array `a`
and store the result in `x`. A polyalgorithm is used for sampling.
Sampling probabilities are proportional to the weights given in `wv`,
if provided. `replace` dictates whether sampling is performed with
replacement and `order` dictates whether an ordered sample, also called
a sequential sample, should be taken.
"""
function sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray;
                 replace::Bool=true, ordered::Bool=false)
    n = length(a)
    k = length(x)
    k == 0 && return x

    if replace  # with replacement
        if ordered
            sort!(direct_sample!(rng, a, x))
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
    sample!(Base.GLOBAL_RNG, a, x; replace=replace, ordered=ordered)


"""
    sample([rng], a, [wv::WeightVec], n::Integer; replace=true, ordered=false)

Select a random, optionally weighted sample of size `n` from an array `a`
using a polyalgorithm. Sampling probabilities are proportional to the weights
given in `wv`, if provided. `replace` dictates whether sampling is performed
with replacement and `order` dictates whether an ordered sample, also called
a sequential sample, should be taken.
"""
function sample{T}(rng::AbstractRNG, a::AbstractArray{T}, n::Integer;
                   replace::Bool=true, ordered::Bool=false)
    sample!(rng, a, Vector{T}(n); replace=replace, ordered=ordered)
end
sample(a::AbstractArray, n::Integer; replace::Bool=true, ordered::Bool=false) =
    sample(Base.GLOBAL_RNG, a, n; replace=replace, ordered=ordered)


"""
    sample([rng], a, [wv::WeightVec], dims::Dims; replace=true, ordered=false)

Select a random, optionally weighted sample from an array `a` specifying
the dimensions `dims` of the output array. Sampling probabilities are
proportional to the weights given in `wv`, if provided. `replace` dictates
whether sampling is performed with replacement and `order` dictates whether
an ordered sample, also called a sequential sample, should be taken.
"""
function sample{T}(rng::AbstractRNG, a::AbstractArray{T}, dims::Dims;
                   replace::Bool=true, ordered::Bool=false)
    sample!(rng, a, Array{T}(dims), rng; replace=replace, ordered=ordered)
end
sample(a::AbstractArray, dims::Dims; replace::Bool=true, ordered::Bool=false) =
    sample(Base.GLOBAL_RNG, a, dims; replace=replace, ordered=ordered)

################################################################
#
#  Weighted sampling
#
################################################################

"""
    sample([rng], wv::WeightVec)

Select a single random integer in `1:length(wv)` with probabilities
proportional to the weights given in `wv`.
"""
function sample(rng::AbstractRNG, wv::WeightVec)
    t = rand(rng) * sum(wv)
    w = values(wv)
    n = length(w)
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end
sample(wv::WeightVec) = sample(Base.GLOBAL_RNG, wv)

sample(rng::AbstractRNG, a::AbstractArray, wv::WeightVec) = a[sample(rng, wv)]
sample(a::AbstractArray, wv::WeightVec) = sample(Base.GLOBAL_RNG, a, wv)

function direct_sample!(rng::AbstractRNG, a::AbstractArray, wv::WeightVec, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("Inconsistent lengths."))
    for i = 1:length(x)
        x[i] = a[sample(rng, wv)]
    end
    return x
end
direct_sample!(a::AbstractArray, wv::WeightVec, x::AbstractArray) =
    direct_sample!(Base.GLOBAL_RNG, a, wv, x)

function make_alias_table!(w::AbstractVector{Float64}, wsum::Float64,
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
    # Note: a and w can be the same way, then that away will be
    #       overriden inplace by acceptance probabilities
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

    larges = Vector{Int}(n)
    smalls = Vector{Int}(n)
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

function alias_sample!(rng::AbstractRNG, a::AbstractArray, wv::WeightVec, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("Inconsistent lengths."))

    # create alias table
    ap = Vector{Float64}(n)
    alias = Vector{Int}(n)
    make_alias_table!(values(wv), sum(wv), ap, alias)

    # sampling
    s = RandIntSampler(n)
    for i = 1:length(x)
        j = rand(rng, s)
        x[i] = rand(rng) < ap[j] ? a[j] : a[alias[j]]
    end
    return x
end
alias_sample!(a::AbstractArray, wv::WeightVec, x::AbstractArray) =
    alias_sample!(Base.GLOBAL_RNG, a, wv, x)

function naive_wsample_norep!(rng::AbstractRNG, a::AbstractArray, wv::WeightVec, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("Inconsistent lengths."))
    k = length(x)

    w = Vector{Float64}(n)
    copy!(w, values(wv))
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
naive_wsample_norep!(a::AbstractArray, wv::WeightVec, x::AbstractArray) =
    naive_wsample_norep!(Base.GLOBAL_RNG, a, wv, x)

# Weighted sampling without replacement
#
# Algorithm A from:
#     Efraimidis PS, Spirakis PG (2006). "Weighted random sampling with a reservoir."
#     Information Processing  Letters, 97 (5), 181-185. ISSN 0020-0190.
#     doi:10.1016/j.ipl.2005.11.003.
#     URL http://www.sciencedirect.com/science/article/pii/S002001900500298X
#
# Instead of keys u^(1/w) where u = random(0,1) keys w/v where v = randexp(1) are used.
function efraimidis_a_wsample_norep!(rng::AbstractRNG, a::AbstractArray, wv::WeightVec, x::AbstractArray)
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
efraimidis_a_wsample_norep!(a::AbstractArray, wv::WeightVec, x::AbstractArray) =
    efraimidis_a_wsample_norep!(Base.GLOBAL_RNG, a, wv, x)

# Weighted sampling without replacement
#
# Algorithm A-Res from:
#     Efraimidis PS, Spirakis PG (2006). "Weighted random sampling with a reservoir."
#     Information Processing  Letters, 97 (5), 181-185. ISSN 0020-0190.
#     doi:10.1016/j.ipl.2005.11.003.
#     URL http://www.sciencedirect.com/science/article/pii/S002001900500298X
#
# Instead of keys u^(1/w) where u = random(0,1) keys w/v where v = randexp(1) are used.
function efraimidis_ares_wsample_norep!(rng::AbstractRNG, a::AbstractArray, wv::WeightVec, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("a and wv must be of same length (got $n and $(length(wv)))."))
    k = length(x)
    k > 0 || return x

    # initialize priority queue
    pq = Vector{Pair{Float64,Int}}(k)
    i = 0
    s = 0
    @inbounds for s in 1:n
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
efraimidis_ares_wsample_norep!(a::AbstractArray, wv::WeightVec, x::AbstractArray) =
    efraimidis_ares_wsample_norep!(Base.GLOBAL_RNG, a, wv, x)

# Weighted sampling without replacement
#
# Algorithm A-ExpJ from:
#     Efraimidis PS, Spirakis PG (2006). "Weighted random sampling with a reservoir."
#     Information Processing  Letters, 97 (5), 181-185. ISSN 0020-0190.
#     doi:10.1016/j.ipl.2005.11.003.
#     URL http://www.sciencedirect.com/science/article/pii/S002001900500298X
#
# Instead of keys u^(1/w) where u = random(0,1) keys w/v where v = randexp(1) are used.
function efraimidis_aexpj_wsample_norep!(rng::AbstractRNG, a::AbstractArray, wv::WeightVec, x::AbstractArray)
    n = length(a)
    length(wv) == n || throw(DimensionMismatch("a and wv must be of same length (got $n and $(length(wv)))."))
    k = length(x)
    k > 0 || return x

    # initialize priority queue
    pq = Vector{Pair{Float64,Int}}(k)
    i = 0
    s = 0
    @inbounds for s in 1:n
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

    # fill output array with items in descending order
    @inbounds for i in k:-1:1
        x[i] = a[heappop!(pq).second]
    end
    return x
end
efraimidis_aexpj_wsample_norep!(a::AbstractArray, wv::WeightVec, x::AbstractArray) =
    efraimidis_aexpj_wsample_norep!(Base.GLOBAL_RNG, a, wv, x)

function sample!(rng::AbstractRNG, a::AbstractArray, wv::WeightVec, x::AbstractArray;
                 replace::Bool=true, ordered::Bool=false)
    n = length(a)
    k = length(x)

    if replace
        if ordered
            sort!(direct_sample!(rng, a, wv, x))
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
        k <= n || error("Cannot draw $n samples from $k samples without replacement.")

        efraimidis_aexpj_wsample_norep!(rng, a, wv, x)
        if ordered
            sort!(x)
        end
    end
    return x
end
sample!(a::AbstractArray, wv::WeightVec, x::AbstractArray) = sample!(Base.GLOBAL_RNG, a, wv, x)

sample{T}(rng::AbstractRNG, a::AbstractArray{T}, wv::WeightVec, n::Integer;
          replace::Bool=true, ordered::Bool=false) =
    sample!(rng, a, wv, Vector{T}(n); replace=replace, ordered=ordered)
sample(a::AbstractArray, wv::WeightVec, n::Integer; replace::Bool=true, ordered::Bool=false) =
    sample(Base.GLOBAL_RNG, a, wv, n; replace=replace, ordered=ordered)

sample{T}(rng::AbstractRNG, a::AbstractArray{T}, wv::WeightVec, dims::Dims;
          replace::Bool=true, ordered::Bool=false) =
    sample!(rng, a, wv, Array{T}(dims); replace=replace, ordered=ordered)
sample(a::AbstractArray, wv::WeightVec, dims::Dims; replace::Bool=true, ordered::Bool=false) =
    sample!(Base.GLOBAL_RNG, a, wv, Array{T}(dims); replace=replace, ordered=ordered)

# wsample interface

"""
    wsample!([rng], a, w, x; replace=true, ordered=false)

Select a weighted sample from an array `a` and store the result in `x`. Sampling
probabilities are proportional to the weights given in `w`. `replace` dictates
whether sampling is performed with replacement and `order` dictates whether an
ordered sample, also called a sequential sample, should be taken.
"""
wsample!(rng::AbstractRNG, a::AbstractArray, w::RealVector, x::AbstractArray;
         replace::Bool=true, ordered::Bool=false) =
    sample!(rng, a, weights(w), x; replace=replace, ordered=ordered)
wsample!(a::AbstractArray, w::RealVector, x::AbstractArray;
         replace::Bool=true, ordered::Bool=false) =
    sample!(Base.GLOBAL_RNG, a, weights(w), x; replace=replace, ordered=ordered)

"""
    wsample([rng], [a], w)

Select a weighted random sample of size 1 from `a` with probabilities proportional
to the weights given in `w`. If `a` is not present, select a random weight from `w`.
"""
wsample(rng::AbstractRNG, w::RealVector) = sample(rng, weights(w))
wsample(w::RealVector) = wsample(Base.GLOBAL_RNG, w)
wsample(rng::AbstractRNG, a::AbstractArray, w::RealVector) = sample(rng, a, weights(w))
wsample(a::AbstractArray, w::RealVector) = wsample(Base.GLOBAL_RNG, a, w)


"""
    wsample([rng], [a], w, n::Integer; replace=true, ordered=false)

Select a weighted random sample of size `n` from `a` with probabilities proportional
to the weights given in `w` if `a` is present, otherwise select a random sample of size
`n` of the weights given in `w`. `replace` dictates whether sampling is performed with
replacement and `order` dictates whether an ordered sample, also called a sequential
sample, should be taken.
"""
wsample{T}(rng::AbstractRNG, a::AbstractArray{T}, w::RealVector, n::Integer;
           replace::Bool=true, ordered::Bool=false) =
    wsample!(rng, a, w, Vector{T}(n); replace=replace, ordered=ordered)
wsample(a::AbstractArray, w::RealVector, n::Integer; replace::Bool=true, ordered::Bool=false) =
    wsample(Base.GLOBAL_RNG, a, w, n; replace=replace, ordered=ordered)

"""
    wsample([rng], [a], w, dims::Dims; replace=true, ordered=false)

Select a weighted random sample from `a` with probabilities proportional to the
weights given in `w` if `a` is present, otherwise select a random sample of size
`n` of the weights given in `w`. The dimensions of the output are given by `dims`.
"""
wsample{T}(rng::AbstractRNG, a::AbstractArray{T}, w::RealVector, dims::Dims;
           replace::Bool=true, ordered::Bool=false) =
    wsample!(rng, a, w, Array{T}(dims); replace=replace, ordered=ordered)
wsample(a::AbstractArray, w::RealVector, dims::Dims; replace::Bool=true, ordered::Bool=false) =
    wsample(Base.GLOBAL_RNG, a, w, dims; replace=replace, ordered=ordered)
