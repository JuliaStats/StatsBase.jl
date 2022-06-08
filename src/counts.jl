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
    addcounts!(r, x, levels::UnitRange{<:Integer}, [wv::AbstractWeights])

Add the number of occurrences in `x` of each value in `levels` to an existing
array `r`. For each `xi ∈ x`, if `xi == levels[j]`, then we increment `r[j]`.

If a weighting vector `wv` is specified, the sum of weights is used rather than the
raw counts.
"""
function addcounts!(r::AbstractArray, x::IntegerArray, levels::IntUnitRange)
    # add counts of integers from x that fall within levels to r

    checkbounds(r, axes(levels)...)

    m0 = first(levels)
    m1 = last(levels)
    b = m0 - firstindex(levels) # firstindex(levels) == 1 because levels::IntUnitRange

    @inbounds for xi in x
        if m0 <= xi <= m1
            r[xi - b] += 1
        end
    end
    return r
end

function addcounts!(r::AbstractArray, x::IntegerArray, levels::IntUnitRange, wv::AbstractWeights)
    # add wv weighted counts of integers from x that fall within levels to r

    length(x) == length(wv) ||
        throw(DimensionMismatch("x and wv must have the same length, got $(length(x)) and $(length(wv))"))

    xv = vec(x) # discard shape because weights() discards shape

    checkbounds(r, axes(levels)...)

    m0 = first(levels)
    m1 = last(levels)
    b = m0 - 1

    @inbounds for i in eachindex(xv, wv)
        xi = xv[i]
        if m0 <= xi <= m1
            r[xi - b] += wv[i]
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

If a vector of weights `wv` is provided, the proportion of weights is computed rather
than the proportion of raw counts.

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
Equivalent to `counts(x, levels) / length(x)`.

If a vector of weights `wv` is provided, the proportion of weights is computed rather
than the proportion of raw counts.
"""
proportions(x::IntegerArray, levels::IntUnitRange) = counts(x, levels) .* inv(length(x))
proportions(x::IntegerArray, levels::IntUnitRange, wv::AbstractWeights) =
    counts(x, levels, wv) .* inv(sum(wv))

"""
    proportions(x, k::Integer, [wv::AbstractWeights])

Return the proportion of integers in 1 to `k` that occur in `x`.

If a vector of weights `wv` is provided, the proportion of weights is computed rather
than the proportion of raw counts.
"""
proportions(x::IntegerArray, k::Integer) = proportions(x, 1:k)
proportions(x::IntegerArray, k::Integer, wv::AbstractWeights) = proportions(x, 1:k, wv)
proportions(x::IntegerArray) = proportions(x, span(x))
proportions(x::IntegerArray, wv::AbstractWeights) = proportions(x, span(x), wv)

#### functions for counting a single list of integers (2D)

function addcounts!(r::AbstractArray, x::IntegerArray, y::IntegerArray, levels::NTuple{2,IntUnitRange})
    # add counts of pairs from zip(x,y) to r

    xlevels, ylevels = levels


    checkbounds(r, axes(xlevels, 1), axes(ylevels, 1))

    mx0 = first(xlevels)
    mx1 = last(xlevels)
    my0 = first(ylevels)
    my1 = last(ylevels)

    bx = mx0 - 1
    by = my0 - 1

    for i in eachindex(vec(x), vec(y))
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
    # add counts of pairs from zip(x,y) to r

    length(x) == length(y) == length(wv) ||
        throw(DimensionMismatch("x, y, and wv must have the same length, but got $(length(x)), $(length(y)), and $(length(wv))"))

    axes(x) == axes(y) ||
        throw(DimensionMismatch("x and y must have the same axes, but got $(axes(x)) and $(axes(y))"))

    xv, yv = vec(x), vec(y) # discard shape because weights() discards shape

    xlevels, ylevels = levels

    checkbounds(r, axes(xlevels, 1), axes(ylevels, 1))

    mx0 = first(xlevels)
    mx1 = last(xlevels)
    my0 = first(ylevels)
    my1 = last(ylevels)

    bx = mx0 - 1
    by = my0 - 1

    for i in eachindex(xv, yv, wv)
        xi = xv[i]
        yi = yv[i]
        if (mx0 <= xi <= mx1) && (my0 <= yi <= my1)
            r[xi - bx, yi - by] += wv[i]
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
    addcounts!(dict, x; alg = :auto)
    addcounts!(dict, x, wv)

Add counts based on `x` to a count map. New entries will be added if new values come up.

If a weighting vector `wv` is specified, the sum of the weights is used rather than the
raw counts.

`alg` is only allowed for unweighted counting and can be one of:
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
addcounts!(cm::Dict, x; alg = :auto) = _addcounts!(eltype(x), cm, x, alg = alg)

function _addcounts!(::Type{T}, cm::Dict, x; alg = :auto) where T
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
function addcounts_dict!(cm::Dict{T}, x) where T
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

# If the bits type is of small size i.e. it can have up to 65536 distinct values
# then it is always better to apply a counting-sort like reduce algorithm for
# faster results and less memory usage. However we still wish to enable others
# to write generic algorithms, therefore the methods below still accept the
# `alg` argument but it is ignored.
function _addcounts!(::Type{Bool}, cm::Dict{Bool}, x::AbstractArray{Bool}; alg = :ignored)
    sumx = sum(x)
    cm[true] = get(cm, true, 0) + sumx
    cm[false] = get(cm, false, 0) + length(x) - sumx
    cm
end

# specialized for `Bool` iterator
function _addcounts!(::Type{Bool}, cm::Dict{Bool}, x; alg = :ignored)
    sumx = 0
    len = 0
    for i in x
        sumx += i
        len += 1
    end
    cm[true] = get(cm, true, 0) + sumx
    cm[false] = get(cm, false, 0) + len - sumx
    cm
end

function _addcounts!(::Type{T}, cm::Dict{T}, x; alg = :ignored) where T <: Union{UInt8, UInt16, Int8, Int16}
    counts = zeros(Int, 2^(8sizeof(T)))

    @inbounds for xi in x
        counts[Int(xi) - typemin(T) + 1] += 1
    end

    for (i, c) in zip(typemin(T):typemax(T), counts)
        if c != 0
            index = ht_keyindex2!(cm, i)
            if index > 0
                @inbounds cm.vals[index] += c
            else
                @inbounds Base._setindex!(cm, c, i, -index)
            end
        end
    end
    cm
end

const BaseRadixSortSafeTypes = Union{Int8, Int16, Int32, Int64, Int128,
                                     UInt8, UInt16, UInt32, UInt64, UInt128,
                                     Float32, Float64}

"Can the type be safely sorted by radixsort"
radixsort_safe(::Type{T}) where T = T<:BaseRadixSortSafeTypes

function _addcounts_radix_sort_loop!(cm::Dict{T}, sx::AbstractVector{T}) where T
    isempty(sx) && return cm
    last_sx = first(sx)
    start_i = firstindex(sx)

    # now the data is sorted: can just run through and accumulate values before
    # adding into the Dict
    @inbounds for i in start_i+1:lastindex(sx)
        sxi = sx[i]
        if last_sx != sxi
            cm[last_sx] = get(cm, last_sx, 0) + i - start_i
            last_sx = sxi
            start_i = i
        end
    end

    last_sx = last(sx)
    cm[last_sx] = get(cm, last_sx, 0) + lastindex(sx) + 1 - start_i

    return cm
end

function _alg(x::AbstractArray)
    @static if VERSION >= v"1.9.0-DEV"
        return Base.DEFAULT_UNSTABLE
    else
        firstindex(x) == 1 ||
            throw(ArgumentError("alg = :radixsort requires either one based indexing or Julia >= 1.9. " *
                                "Use `alg = :dict` as an alternative."))
        return SortingAlgorithms.RadixSort
    end
end

function addcounts_radixsort!(cm::Dict{T}, x::AbstractArray{T}) where T
    # sort the x using radixsort
    sx = sort(vec(x), alg=_alg(x))

    # Delegate the loop to a separate function since sort might not
    # be inferred in Julia 0.6 after SortingAlgorithms is loaded.
    # It seems that sort is inferred in Julia 0.7.
    return _addcounts_radix_sort_loop!(cm, sx)
end

# fall-back for `x` an iterator
function addcounts_radixsort!(cm::Dict{T}, x) where T
    cx = vec(collect(x))
    sx = sort!(cx, alg = _alg(cx))
    return _addcounts_radix_sort_loop!(cm, sx)
end

function addcounts!(cm::Dict{T}, x::AbstractArray{T}, wv::AbstractVector{W}) where {T,W<:Real}
    # add wv weighted counts of integers from x to cm

    length(x) == length(wv) ||
        throw(DimensionMismatch("x and wv must have the same length, got $(length(x)) and $(length(wv))"))

    xv = vec(x) # discard shape because weights() discards shape

    z = zero(W)

    for i in eachindex(xv, wv)
        @inbounds xi = xv[i]
        @inbounds wi = wv[i]
        cm[xi] = get(cm, xi, z) + wi
    end
    return cm
end


"""
    countmap(x; alg = :auto)
    countmap(x::AbstractVector, wv::AbstractVector{<:Real})

Return a dictionary mapping each unique value in `x` to its number of occurrences.

If a weighting vector `wv` is specified, the sum of weights is used rather than the
raw counts.

`alg` is only allowed for unweighted counting and can be one of:
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
countmap(x; alg = :auto) = addcounts!(Dict{eltype(x),Int}(), x; alg = alg)
countmap(x::AbstractArray{T}, wv::AbstractVector{W}) where {T,W<:Real} = addcounts!(Dict{T,W}(), x, wv)


"""
    proportionmap(x)
    proportionmap(x::AbstractVector, w::AbstractVector{<:Real})

Return a dictionary mapping each unique value in `x` to its proportion in `x`.

If a vector of weights `wv` is provided, the proportion of weights is computed rather
than the proportion of raw counts.
"""
proportionmap(x::AbstractArray) = _normalize_countmap(countmap(x), length(x))
proportionmap(x::AbstractArray, wv::AbstractWeights) = _normalize_countmap(countmap(x, wv), sum(wv))
