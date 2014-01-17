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
        if n == 1
            for j = 1 : rk
                @inbounds x[offset + j] = a[i]
            end
            offset = k
        else
            m = rand_binom(rk, 1.0 / n)
            for j = 1 : m
                @inbounds x[offset + j] = a[i]
            end
            i += 1
            n -= 1
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

function wsample(ws::AbstractArray; wsum::Float64 = sum(ws))
    t = rand() * wsum
    i = 0
    p = 0.
    while p < t
        i += 1
        p += ws[i]
    end
    return i
end

wsample(xs::AbstractArray, ws::AbstractArray; wsum::Float64 = sum(ws)) =
  xs[wsample(ws, wsum=wsum)]

# Author: Mike Innes
function ordered_wsample!(xs::AbstractArray, ws::AbstractArray, target::AbstractArray; wsum::Float64 = sum(ws))
    n = length(xs)
    k = length(target)
    length(ws) == n || throw(ArgumentError("Inconsistent argument dimensions."))

    j = 0
    for i = 1:n
        k > 0 || break
        num = i == n ? k : rand_binom(k, ws[i] / wsum) 
        for _ = 1:num
            j += 1
            target[j] = xs[i]
        end
        k -= num
        wsum -= ws[i]
    end
    return target
end

function wsample!(xs::AbstractArray, ws::AbstractArray, target::AbstractArray; wsum::Float64 = sum(ws), ordered::Bool = false)
    k = length(target)
    ordered && return ordered_wsample!(xs, ws, target, wsum = wsum)
    k > 100 && return ordered_wsample!(xs, ws, target, wsum = wsum) |> shuffle!

    for i = 1:k
        target[i] = wsample(xs)
    end
    return target
end

wsample(xs::AbstractArray, ws::AbstractArray, k; wsum::Float64 = sum(ws), ordered::Bool = false) =
  wsample!(xs, ws, similar(xs, k), wsum = wsum, ordered = ordered)

