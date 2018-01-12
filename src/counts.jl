# Counts of discrete values

#################################################
#
#  counts on given levels
#
#################################################

const IntUnitRange{T<:Integer} = UnitRange{T}

if isdefined(Base, :ht_keyindex2)
    const ht_keyindex2! = Base.ht_keyindex2
else
    using Base: ht_keyindex2!
end

#### functions for counting a single list of integers (1D)
"""
    addcounts!(r, x, levels::UnitRange{<:Int}, [wv::AbstractWeights])

Add the number of occurrences in `x` of each value in `levels` to an existing
array `r`. If a weighting vector `wv` is specified, the sum of weights is used
rather than the raw counts.
"""
function addcounts!(r::AbstractArray, x::IntegerArray, levels::IntUnitRange)
    # add counts of integers from x to r

    k = length(levels)
    length(r) == k || throw(DimensionMismatch())

    m0 = levels[1]
    m1 = levels[end]
    b = m0 - 1

    @inbounds for i in 1 : length(x)
        xi = x[i]
        if m0 <= xi <= m1
            r[xi - b] += 1
        end
    end
    return r
end

function addcounts!(r::AbstractArray, x::IntegerArray, levels::IntUnitRange, wv::AbstractWeights)
    k = length(levels)
    length(r) == k || throw(DimensionMismatch())

    m0 = levels[1]
    m1 = levels[end]
    b = m0 - 1
    w = values(wv)

    @inbounds for i in 1 : length(x)
        xi = x[i]
        if m0 <= xi <= m1
            r[xi - b] += w[i]
        end
    end
    return r
end


"""
    counts(x, [wv::AbstractWeights])
    counts(x, levels::UnitRange{<:Integer}, [wv::AbstractWeights])
    counts(x, k::Integer, [wv::AbstractWeights])

Count the number of times each value in `x` occurs. If `levels` is provided, only values
falling in that range will be considered (the others will be ignored without
raising an error or a warning). If an integer `k` is provided, only values in the
range `1:k` will be considered.

If a weighting vector `wv` is specified, the sum of the weights is used rather than the
raw counts.

The output is a vector of length `length(levels)`.
"""
function counts end

counts(x::IntegerArray, levels::IntUnitRange) =
    addcounts!(zeros(Int, length(levels)), x, levels)
counts(x::IntegerArray, levels::IntUnitRange, wv::AbstractWeights) =
    addcounts!(zeros(eltype(wv), length(levels)), x, levels, wv)
counts(x::IntegerArray, k::Integer) = counts(x, 1:k)
counts(x::IntegerArray, k::Integer, wv::AbstractWeights) = counts(x, 1:k, wv)
counts(x::IntegerArray) = counts(x, span(x))
counts(x::IntegerArray, wv::AbstractWeights) = counts(x, span(x), wv)


"""
    proportions(x, levels=span(x), [wv::AbstractWeights])

Return the proportion of values in the range `levels` that occur in `x`.
Equivalent to `counts(x, levels) / length(x)`. If a weighting vector `wv`
is specified, the sum of the weights is used rather than the raw counts.
"""
proportions(x::IntegerArray, levels::IntUnitRange) = counts(x, levels) .* inv(length(x))
proportions(x::IntegerArray, levels::IntUnitRange, wv::AbstractWeights) =
    counts(x, levels, wv) .* inv(sum(wv))

"""
    proportions(x, k::Integer, [wv::AbstractWeights])

Return the proportion of integers in 1 to `k` that occur in `x`.
"""
proportions(x::IntegerArray, k::Integer) = proportions(x, 1:k)
proportions(x::IntegerArray, k::Integer, wv::AbstractWeights) = proportions(x, 1:k, wv)
proportions(x::IntegerArray) = proportions(x, span(x))
proportions(x::IntegerArray, wv::AbstractWeights) = proportions(x, span(x), wv)

#### functions for counting a single list of integers (2D)

function addcounts!(r::AbstractArray, x::IntegerArray, y::IntegerArray, levels::NTuple{2,IntUnitRange})
    # add counts of integers from x to r

    n = length(x)
    length(y) == n || throw(DimensionMismatch())

    xlevels, ylevels = levels

    kx = length(xlevels)
    ky = length(ylevels)
    size(r) == (kx, ky) || throw(DimensionMismatch())

    mx0 = xlevels[1]
    mx1 = xlevels[end]
    my0 = ylevels[1]
    my1 = ylevels[end]

    bx = mx0 - 1
    by = my0 - 1

    for i = 1:n
        xi = x[i]
        yi = y[i]
        if (mx0 <= xi <= mx1) && (my0 <= yi <= my1)
            r[xi - bx, yi - by] += 1
        end
    end
    return r
end

function addcounts!(r::AbstractArray, x::IntegerArray, y::IntegerArray,
                    levels::NTuple{2,IntUnitRange}, wv::AbstractWeights)
    # add counts of integers from x to r

    n = length(x)
    length(y) == length(wv) == n || throw(DimensionMismatch())

    xlevels, ylevels = levels

    kx = length(xlevels)
    ky = length(ylevels)
    size(r) == (kx, ky) || throw(DimensionMismatch())

    mx0 = xlevels[1]
    mx1 = xlevels[end]
    my0 = ylevels[1]
    my1 = ylevels[end]

    bx = mx0 - 1
    by = my0 - 1
    w = values(wv)

    for i = 1:n
        xi = x[i]
        yi = y[i]
        if (mx0 <= xi <= mx1) && (my0 <= yi <= my1)
            r[xi - bx, yi - by] += w[i]
        end
    end
    return r
end

# facet functions

function counts(x::IntegerArray, y::IntegerArray, levels::NTuple{2,IntUnitRange})
    addcounts!(zeros(Int, length(levels[1]), length(levels[2])), x, y, levels)
end

function counts(x::IntegerArray, y::IntegerArray, levels::NTuple{2,IntUnitRange}, wv::AbstractWeights)
    addcounts!(zeros(eltype(wv), length(levels[1]), length(levels[2])), x, y, levels, wv)
end

counts(x::IntegerArray, y::IntegerArray, levels::IntUnitRange) =
    counts(x, y, (levels, levels))
counts(x::IntegerArray, y::IntegerArray, levels::IntUnitRange, wv::AbstractWeights) =
    counts(x, y, (levels, levels), wv)

counts(x::IntegerArray, y::IntegerArray, ks::NTuple{2,Integer}) =
    counts(x, y, (1:ks[1], 1:ks[2]))
counts(x::IntegerArray, y::IntegerArray, ks::NTuple{2,Integer}, wv::AbstractWeights) =
    counts(x, y, (1:ks[1], 1:ks[2]), wv)
counts(x::IntegerArray, y::IntegerArray, k::Integer) = counts(x, y, (1:k, 1:k))
counts(x::IntegerArray, y::IntegerArray, k::Integer, wv::AbstractWeights) =
    counts(x, y, (1:k, 1:k), wv)
counts(x::IntegerArray, y::IntegerArray) = counts(x, y, (span(x), span(y)))
counts(x::IntegerArray, y::IntegerArray, wv::AbstractWeights) = counts(x, y, (span(x), span(y)), wv)

proportions(x::IntegerArray, y::IntegerArray, levels::NTuple{2,IntUnitRange}) =
    counts(x, y, levels) .* inv(length(x))
proportions(x::IntegerArray, y::IntegerArray, levels::NTuple{2,IntUnitRange}, wv::AbstractWeights) =
    counts(x, y, levels, wv) .* inv(sum(wv))

proportions(x::IntegerArray, y::IntegerArray, ks::NTuple{2,Integer}) =
    proportions(x, y, (1:ks[1], 1:ks[2]))
proportions(x::IntegerArray, y::IntegerArray, ks::NTuple{2,Integer}, wv::AbstractWeights) =
    proportions(x, y, (1:ks[1], 1:ks[2]), wv)
proportions(x::IntegerArray, y::IntegerArray, k::Integer) = proportions(x, y, (1:k, 1:k))
proportions(x::IntegerArray, y::IntegerArray, k::Integer, wv::AbstractWeights) =
    proportions(x, y, (1:k, 1:k), wv)
proportions(x::IntegerArray, y::IntegerArray) = proportions(x, y, (span(x), span(y)))
proportions(x::IntegerArray, y::IntegerArray, wv::AbstractWeights) =
    proportions(x, y, (span(x), span(y)), wv)


#################################################
#
#  countmap on unknown levels
#
#  These methods are based on dictionaries, and
#  can be used on any kind of hashable values.
#
#################################################

## auxiliary functions

function _normalize_countmap(cm::Dict{T}, s::Real) where T
    r = Dict{T,Float64}()
    for (k, c) in cm
        r[k] = c / s
    end
    return r
end

## 1D


"""
    addcounts!(dict, x[, wv]; alg = :auto)

Add counts based on `x` to a count map. New entries will be added if new values come up.
If a weighting vector `wv` is specified, the sum of the weights is used rather than the
raw counts.

`alg` can be one of:
- `:auto` (default): if `StatsBase.radixsort_safe(eltype(x)) == true` then use
                     `:radixsort`, otherwise use `:dict`.

- `:radixsort`:      if `radixsort_safe(eltype(x)) == true` then use the
                     [radix sort](https://en.wikipedia.org/wiki/Radix_sort)
                     algorithm to sort the input vector which will generally lead to
                     shorter running time. However the radix sort algorithm creates a
                     copy of the input vector and hence uses more RAM. Choose `:dict`
                     if the amount of available RAM is a limitation.

- `:dict`:           use `Dict`-based method which is generally slower but uses less
                     RAM and is safe for any data type.
"""
function addcounts!(cm::Dict{T}, x::AbstractArray{T}; alg = :auto) where T
    # if it's safe to be sorted using radixsort then it should be faster
    # albeit using more RAM
    if radixsort_safe(T) && (alg == :auto || alg == :radixsort)
        addcounts_radixsort!(cm, x)
    elseif alg == :radixsort
        throw(ArgumentError("`alg = :radixsort` is chosen but type `radixsort_safe($T)` did not return `true`; use `alg = :auto` or `alg = :dict` instead"))
    else
        addcounts_dict!(cm,x)
    end
    return cm
end

"""Dict-based addcounts method"""
function addcounts_dict!(cm::Dict{T}, x::AbstractArray{T}) where T
    for v in x
        index = ht_keyindex2!(cm, v)
        if index > 0
            @inbounds cm.vals[index] += 1
        else
            @inbounds Base._setindex!(cm, 1, v, -index)
        end
    end
    return cm
end

"""Specialist addcounts methods for small bits types"""
function addcounts!(cm::Dict{Bool}, x::AbstractArray{Bool})
    sumx = sum(x)
    cm[true] = sumx
    cm[false] =length(x) = sumx
    cm
end

"""toindex functions to convert an integer to an array index between 1 to typemax{T}"""
function toindex(x::Int16)
  Int(x) + 32769
end

function toindex(x::Int8)
  Int(x) + 129
end

function toindex(x::T) where T <: Union{UInt8, UInt16}
  Int(x) + 1
end

function addcounts!(cm::Dict{T}, x::Vector{T}) where T <: Union{UInt8, UInt16, Int8, Int16}
  arr = zeros(Int, 2^(sizeof(T)*8))

  for xi in x
    @inbounds arr[toindex(xi)] += 1
  end

  for (i,arr1) in zip(typemin(T):typemax(T),arr)
    arr1 == 0 || @inbounds cm[i] = arr1
  end
  cm
end

const BaseRadixSortSafeTypes = Union{Int8, Int16, Int32, Int64, Int128,
                                     UInt8, UInt16, UInt32, UInt64, UInt128,
                                     Float32, Float64}

"Can the type be safely sorted by radixsort"
radixsort_safe(::Type{T}) where {T<:BaseRadixSortSafeTypes} = true
radixsort_safe(::Type) = false

function addcounts_radixsort!(cm::Dict{T}, x::AbstractArray{T}) where T
    # sort the x using radixsort
    sx = sort(x, alg = RadixSort)::typeof(x)

    tmpcount = 1
    last_sx = sx[1]

    # now the data is sorted: can just run through and accumulate values before
    # adding into the Dict
    for i in 2:length(sx)
        sxi = sx[i]
        if last_sx == sxi
            tmpcount += 1
        else
            cm[last_sx] = tmpcount
            last_sx = sxi
            tmpcount = 1
        end
    end

    cm[sx[end]] = tmpcount

    return cm
end

function addcounts!(cm::Dict{T}, x::AbstractArray{T}, wv::AbstractVector{W}) where {T,W<:Real}
    n = length(x)
    length(wv) == n || throw(DimensionMismatch())
    w = values(wv)
    z = zero(W)

    for i = 1 : n
        @inbounds xi = x[i]
        @inbounds wi = w[i]
        cm[xi] = get(cm, xi, z) + wi
    end
    return cm
end


"""
    countmap(x; alg = :auto)

Return a dictionary mapping each unique value in `x` to its number
of occurrences.

- `:auto` (default): if `StatsBase.radixsort_safe(eltype(x)) == true` then use
                     `:radixsort`, otherwise use `:dict`.

- `:radixsort`:      if `radixsort_safe(eltype(x)) == true` then use the
                     [radix sort](https://en.wikipedia.org/wiki/Radix_sort)
                     algorithm to sort the input vector which will generally lead to
                     shorter running time. However the radix sort algorithm creates a
                     copy of the input vector and hence uses more RAM. Choose `:dict`
                     if the amount of available RAM is a limitation.

- `:dict`:           use `Dict`-based method which is generally slower but uses less
                     RAM and is safe for any data type.
"""
countmap(x::AbstractArray{T}; alg = :auto) where {T} = addcounts!(Dict{T,Int}(), x; alg = alg)
# specialist function for Bool
countmap(x::AbstractArray{T}) where T<:Union{Bool, UInt8, UInt16, Int8, Int16} = addcounts!(Dict{T,Int}(), x) 
    
countmap(x::AbstractArray{T}, wv::AbstractVector{W}) where {T,W<:Real} = addcounts!(Dict{T,W}(), x, wv)


"""
    proportionmap(x)

Return a dictionary mapping each unique value in `x` to its
proportion in `x`.
"""
proportionmap(x::AbstractArray) = _normalize_countmap(countmap(x), length(x))
proportionmap(x::AbstractArray, wv::AbstractWeights) = _normalize_countmap(countmap(x, wv), sum(wv))
