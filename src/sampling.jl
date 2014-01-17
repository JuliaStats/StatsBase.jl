# Functions for sampling from populations

### draw a pair of distinct integers in [1:n] 

function samplepair(n::Int)
    i1 = randi(n)
    i2 = randi(n-1)
    return (i1, i2 == i1 ? n : i2)
end 

function samplepair(a::AbstractArray)
    i1, i2 = samplepair(length(a))
    return a[i1], a[i2]
end


### internal sampling algorithms
 
# Fisher-Yates sampling
immutable FisherYatesSampler
    n::Int 
    seq::Vector{Int}   # Internal sequence for shuffling

    FisherYatesSampler(n::Int) = new(n, [1:n])
end

function rand!(s::FisherYatesSampler, a::AbstractArray, x::AbstractArray)
    # draw samples without-replacement to x

    n::Int = s.n
    k::Int = length(x)
    if k > n
        throw(ArgumentError("Cannot draw more than n samples without replacement."))
    end

    seq::Vector{Int} = s.seq
    for i = 1:k
        j = randi(i, n)
        sj = seq[j]
        x[i] = a[sj]
        seq[j] = seq[i]
        seq[i] = sj
    end
    x
end

fisher_yates_sample!(a::AbstractArray, x::AbstractArray) = rand!(FisherYatesSampler(length(a)), a, x)

# self-avoiding sampling
function self_avoid_sample!{T}(a::AbstractArray{T}, x::AbstractArray)
    # This algorithm is suitable when length(x) << length(a)

    s = Set{T}()
    # sizehint(s, length(x))
    rgen = RandIntSampler(length(a))

    # first one    
    idx = rand(rgen)
    x[1] = a[idx]
    push!(s, idx)

    # remaining
    for i = 2:length(x)
        idx = rand(rgen)
        while (idx in s)
            idx = rand(rgen)
        end
        x[i] = a[idx]
        push!(s, idx)
    end
    x
end

# Ordered sampling without replacement
# Original author: Mike Innes

function ordered_sample!(a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    offset = 0
    i = 1
    
    while offset < k
        rk = k - offset
        if i == n
            for j = 1 : rk
                @inbounds x[offset + j] = a[i]
            end
            offset = k
        else
            m = rand_binom(rk, 1.0 / (n - i + 1))
            for j = 1 : m
                @inbounds x[offset + j] = a[i]
            end
            i += 1
            offset += m
        end
    end
    x
end

###########################################################
#
#   Interface functions
#
###########################################################

sample(a::AbstractArray) = a[randi(length(a))]

function sample!(a::AbstractArray, x::AbstractArray; replace::Bool=true, ordered::Bool=false)
    n = length(a)
    k = length(x)

    if isempty(x) 
        return x
    end

    if replace
        if ordered
            ordered_sample!(a, x)
        else
            s = RandIntSampler(n)
            for i = 1:k
                @inbounds x[i] = a[rand(s)]
            end
        end

    else  # without replacement
        k <= n || error("Cannot draw more samples without replacement.")

        if k == 1
            @inbounds x[1] = sample(a)

        elseif k == 2
            @inbounds (x[1], x[2]) = samplepair(a)

        elseif n < k * max(k, 100)
            fisher_yates_sample!(a, x)

        else
            self_avoid_sample!(a, x)
        end

        if ordered
            sort!(x)
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
    isempty(wv) && error("Input weight vector is empty.")
    t = rand() * sum(wv)
    w = values(wv)
    n = length(w)
    i = 1
    @inbounds cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end

sample(a::AbstractArray, wv::WeightVec) = a[sample(wv)]

# Author: Mike Innes
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

function sample!(a::AbstractArray, wv::WeightVec, x::AbstractArray; 
    replace::Bool=true, ordered::Bool=false)

    n = length(a)
    k = length(x)

    if ordered && replace
        ordered_sample!(a, wv, x)
    else
        if k > 100 * n
            shuffle!(ordered_sample!(a, wv, x))
        else
            for i = 1 : k
                @inbounds x[i] = sample(a, wv)
            end
        end

        if ordered
            sort!(x)
        end
    end
    return x
end

sample{T}(a::AbstractArray{T}, wv::WeightVec, n::Integer; replace::Bool=true, ordered::Bool=false) =
    sample!(a, wv, Array(T, n); replace=replace, ordered=ordered)

sample{T}(a::AbstractArray{T}, wv::WeightVec, dims::Dims; replace::Bool=true, ordered::Bool=false) =
    sample!(a, wv, Array(T, dims); replace=replace, ordered=ordered)    

# wsample interface

wsample(w::RealVector) = sample(weights(w))
wsample(a::AbstractArray, w::RealVector) = sample(a, weights(w))

wsample!(a::AbstractArray, w::RealVector, x::AbstractArray; replace::Bool=true, ordered::Bool=false) = 
    sample!(a, weights(w), x; replace=replace, ordered=ordered)

wsample{T}(a::AbstractArray{T}, w::RealVector, n::Integer; replace::Bool=true, ordered::Bool=false) = 
    wsample!(a, w, Array(T, n); replace=replace, ordered=ordered)

wsample{T}(a::AbstractArray{T}, w::RealVector, dims::Dims; replace::Bool=true, ordered::Bool=false) = 
    wsample!(a, w, Array(T, dims); replace=replace, ordered=ordered)    


