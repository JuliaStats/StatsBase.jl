
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
function direct_sample!(a::UnitRange, x::AbstractArray)
    s = RandIntSampler(length(a))
    b = a[1] - 1
    if b == 0
        for i = 1:length(x)
            @inbounds x[i] = rand(s)
        end
    else
        for i = 1:length(x)
            @inbounds x[i] = b + rand(s)
        end
    end
    return x
end

function direct_sample!(a::AbstractArray, x::AbstractArray)
    s = RandIntSampler(length(a))
    for i = 1:length(x)
        @inbounds x[i] = a[rand(s)]
    end
    return x
end

# Expanded Multinomial sampling
#
#   for each element in a, we draw the number of its
#   occurrences, and fill this element to x for 
#   this number of times. 
#
function xmultinom_sample!(a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    offset = 0
    i = 1
    while offset < k
        rk = k - offset
        if i == n
            @inbounds ai = a[i]
            for j = 1:rk
                @inbounds x[offset + j] = ai
            end
            offset = k
        else
            m = rand_binom(rk, 1.0 / (n - i + 1))
            if m > 0
                @inbounds ai = a[i]
                for j = 1:m
                    @inbounds x[offset + j] = a[i]
                end
                offset += m
            end
            i += 1            
        end
    end
    x
end


### draw a pair of distinct integers in [1:n]

function samplepair(n::Int)
    i1 = randi(n)
    i2 = randi(n-1)
    return (i1, ifelse(i2 == i1, n, i2))
end

function samplepair(a::AbstractArray)
    i1, i2 = samplepair(length(a))
    return a[i1], a[i2]
end


### Algorithm for sampling without replacement

# Knuth's Algorithm S
#
#   Reference: The Art of Computer Programming, Vol 2, 3.4.2, p.142
#
function knuths_sample!(a::AbstractArray, x::AbstractArray; initshuffle::Bool=true) 
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    # initialize
    for i = 1:k
        x[i] = a[i]
    end
    if initshuffle
        for j = 1:k
            l = randi(j, k)
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
        if rand() * i < k  # keep it with probability k / i
            x[rand(s)] = a[i]            
        end
    end
    return x
end

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
function fisher_yates_sample!(a::AbstractArray, x::AbstractArray) 
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    inds = Array(Int, n)
    for i = 1:n
        @inbounds inds[i] = i
    end

    @inbounds for i = 1:k
        j = randi(i, n)
        t = inds[j]
        inds[j] = inds[i]
        inds[i] = t
        x[i] = a[t]
    end
    return x
end

# Self-avoid sampling
#
#   Use a set to maintain the index that has been sampled,
#
#   each time draw a new index, if the index has already
#   been sampled, redraw until it draws an unsampled one.  
#
function self_avoid_sample!(a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    s = Set{Int}()
    sizehint(s, k)
    rgen = RandIntSampler(n)

    # first one
    idx = rand(rgen)
    x[1] = a[idx]
    push!(s, idx)

    # remaining
    for i = 2:k
        idx = rand(rgen)
        while idx in s
            idx = rand(rgen)
        end
        x[i] = a[idx]
        push!(s, idx)
    end
    x
end



# Ordered sampling without replacement
# Original author: Mike Innes
function rand_first_index(n, k)
    r = rand()
    p = k/n
    i = 1
    while p < r
        i += 1
        p += (1-p)k/(n-(i-1))
    end
    return i
end

function ordered_sample_norep!(xs::AbstractArray, target::AbstractArray)
    n = length(xs)
    k = length(target)
    i = 0
    for j in 1:k
        step = rand_first_index(n, k)
        n -= step
        i += step
        target[j] = xs[i]
        k -= 1
    end
    return target
end


### Interface functions (poly-algorithms)

sample(a::AbstractArray) = a[randi(length(a))]

function sample!(a::AbstractArray, x::AbstractArray; replace::Bool=true, ordered::Bool=false)
    n = length(a)
    k = length(x)
    n == 0 && return x

    if replace  # with replacement
        if ordered
            if k > 10 * n
                xmultinom_sample!(a, x)
            else
                sort!(direct_sample!(a, x))
            end
        else
            direct_sample!(a, x)
        end

    else  # without replacement
        k <= n || error("Cannot draw more samples without replacement.")

        if ordered
            if k * 20 > n
                ordered_sample_norep!(a, x)
            else
                sort!(sample!(a, x; replace=false, ordered=false))
            end
        else
            if k == 1
                @inbounds x[1] = sample(a)
            elseif k == 2
                @inbounds (x[1], x[2]) = samplepair(a)
            elseif n < k * max(k, 100)
                fisher_yates_sample!(a, x)
            else
                self_avoid_sample!(a, x)
            end
        end
    end
    return x
end

function sample{T}(a::AbstractArray{T}, n::Integer; replace::Bool=true, ordered::Bool=false)
    sample!(a, Array(T, n); replace=replace, ordered=ordered)
end

function sample{T}(a::AbstractArray{T}, dims::Dims; replace::Bool=true, ordered::Bool=false)
    sample!(a, Array(T, dims); replace=replace, ordered=ordered)
end


################################################################
#
#  Weighted sampling
#
################################################################

function sample(wv::WeightVec)
    t = rand() * sum(wv)
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

sample(a::AbstractArray, wv::WeightVec) = a[sample(wv)]

# Original author: Mike Innes
function ordered_sample!(a::AbstractArray, wv::WeightVec, x::AbstractArray)
    n = length(a)
    k = length(x)
    offset = 0
    i = 1
    wsum = sum(wv)
    w = values(wv)

    while offset < k
        rk = k - offset
        wi = w[i]

        if i == n || wi >= wsum
            for j = 1 : rk
                @inbounds x[offset + j] = a[i]
            end
            offset = k
        else
            m = rand_binom(rk, wi / wsum)
            for j = 1 : m
                @inbounds x[offset + j] = a[i]
            end
            i += 1
            wsum -= wi
            offset += m
        end
    end
    x
end


function sample!(a::AbstractArray, wv::WeightVec, x::AbstractArray; ordered::Bool=false)
    n = length(a)
    k = length(x)

    if ordered
        ordered_sample!(a, wv, x)
    else
        if k > 100 * n
            shuffle!(ordered_sample!(a, wv, x))
        else
            for i = 1 : k
                @inbounds x[i] = sample(a, wv)
            end
        end
    end
    return x
end

sample{T}(a::AbstractArray{T}, wv::WeightVec, n::Integer; ordered::Bool=false) =
    sample!(a, wv, Array(T, n); ordered=ordered)

sample{T}(a::AbstractArray{T}, wv::WeightVec, dims::Dims; ordered::Bool=false) =
    sample!(a, wv, Array(T, dims); ordered=ordered)

# wsample interface

wsample(w::RealVector) = sample(weights(w))
wsample(a::AbstractArray, w::RealVector) = sample(a, weights(w))

wsample!(a::AbstractArray, w::RealVector, x::AbstractArray; ordered::Bool=false) =
    sample!(a, weights(w), x; ordered=ordered)

wsample{T}(a::AbstractArray{T}, w::RealVector, n::Integer; ordered::Bool=false) =
    wsample!(a, w, Array(T, n); ordered=ordered)

wsample{T}(a::AbstractArray{T}, w::RealVector, dims::Dims; ordered::Bool=false) =
    wsample!(a, w, Array(T, dims); ordered=ordered)


