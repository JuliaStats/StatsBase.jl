# Computing deviation in a variety of ways

## count the number of equal/non-equal pairs

"""
    counteq(a, b)

Count the number of indices at which the elements of the arrays
`a` and `b` are equal.
"""
function counteq(a::AbstractArray, b::AbstractArray)
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent lengths."))
    c = 0
    for i = 1:n
        @inbounds if a[i] == b[i]
            c += 1
        end
    end
    return c
end


"""
    countne(a, b)

Count the number of indices at which the elements of the arrays
`a` and `b` are not equal.
"""
function countne(a::AbstractArray, b::AbstractArray)
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent lengths."))
    c = 0
    for i = 1:n
        @inbounds if a[i] != b[i]
            c += 1
        end
    end
    return c
end


"""
    sqL2dist(a, b)

Compute the squared L2 distance between two arrays: `sumabs2(a - b)`.
"""
function sqL2dist{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T})
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Input dimension mismatch"))
    r = 0.0
    for i = 1:n
        @inbounds r += abs2(a[i] - b[i])
    end
    return r
end


# L2 distance
"""
    L2dist(a, b)

Compute the L2 distance between two arrays: `sqrt(sumabs2(a - b))`.
"""
L2dist{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = sqrt(sqL2dist(a, b))


# L1 distance
"""
    L1dist(a, b)

Compute the L1 distance between two arrays: `sumabs(a - b)`.
"""
function L1dist{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T})
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Input dimension mismatch"))
    r = 0.0
    for i = 1:n
        @inbounds r += abs(a[i] - b[i])
    end
    return r
end


# Linf distance
"""
    Linfdist(a, b)

Compute the L∞ distance, also called the Chebyshev distance, between
two arrays: `maxabs(a - b)`.
"""
function Linfdist{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T})
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Input dimension mismatch"))
    r = 0.0
    for i = 1:n
        @inbounds v = abs(a[i] - b[i])
        if r < v
            r = v
        end
    end
    return r
end


# Generalized KL-divergence
"""
    gkldiv(a, b)

Compute the generalized Kullback-Leibler divergence between two arrays:
`sum(a*log(a/b)-a+b)`.
"""
function gkldiv{T<:AbstractFloat}(a::AbstractArray{T}, b::AbstractArray{T})
    n = length(a)
    r = 0.0
    for i = 1:n
        @inbounds ai = a[i]
        @inbounds bi = b[i]
        if ai > 0
            r += (ai * log(ai / bi) - ai + bi)
        else
            r += bi
        end
    end
    return r::Float64
end


# MeanAD: mean absolute deviation
"""
    meanad(a, b)

Return the mean absolute deviation between two arrays: `mean(abs(a - b))`.
"""
meanad{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = L1dist(a, b) / length(a)


# MaxAD: maximum absolute deviation
"""
    maxad(a, b)

Return the maximum absolute deviation between two arrays: `maxabs(a - b)`.
"""
maxad{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = Linfdist(a, b)


# MSD: mean squared deviation
"""
    msd(a, b)

Return the mean squared deviation between two arrays: `mean(abs2(a - b))`.
"""
msd{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = sqL2dist(a, b) / length(a)


# RMSD: root mean squared deviation
"""
    rmsd(a, b; normalize=false)

Return the root mean squared deviation between two optionally
normalized arrays. The root mean squared deviation is computed
as `sqrt(msd(a, b))`.
"""
function rmsd{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}; normalize::Bool=false)
    v = sqrt(msd(a, b))
    if normalize
        amin, amax = extrema(a)
        v /= (amax - amin)
    end
    return v
end


# PSNR: peak signal-to-noise ratio
"""
    psnr(a, b, maxv)

Compute the peak signal-to-noise ratio between two arrays `a` and `b`.
`maxv` is the maximum possible value either array can take. The PSNR
is computed as `10 * log10(maxv^2 / msd(a, b))`.
"""
function psnr{T<:Real}(a::AbstractArray{T}, b::AbstractArray{T}, maxv::Real)
    20. * log10(maxv) - 10. * log10(msd(a, b))
end
