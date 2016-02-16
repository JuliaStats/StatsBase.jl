# Computing deviation in a variety of ways

## count the number of equal/non-equal pairs
"""
    counteq(a::AbstractArray, b::AbstractArray)

Count the number of equal pairs of elements in a and b, i.e countnz(a .== b).
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
    countne(a::AbstractArray, b::AbstractArray)

Count the number of non-equal pairs of elements in a and b, i.e countnz(a .!= b).
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

# squared L2 distance
"""
    sqL2dist(a::AbstractArray{T}, b::AbstractArray{T})

Computes the Squared L2 distance between `a` and `b`.
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
    L2dist(a::AbstractArray{T}, b::AbstractArray{T})

Computes the L2 distance between `a` and `b`.
"""
L2dist{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = sqrt(sqL2dist(a, b))

# L1 distance
"""
    L1dist(a::AbstractArray{T}, b::AbstractArray{T})

Computes the L1 distance between `a` and `b`.
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
    Linfdist(a::AbstractArray{T}, b::AbstractArray{T})

Computes the Linf distance between `a` and `b`.
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
    gkldiv(a::AbstractArray{T}, b::AbstractArray{T})

Computes the generalized Kullback-Leibler divergence between two arrays `a` and `b`.
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
    meanad(a::AbstractArray{T}, b::AbstractArray{T})

Computes the mean absolute deviation between `a` and `b`, i.e. `mean(abs(a - b))`.
"""
meanad{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = L1dist(a, b) / length(a)

# MaxAD: maximum absolute deviation
"""
    maxad(a::AbstractArray{T}, b::AbstractArray{T})

Computes the maximum absolute deviation between `a` and `b`, i.e. `maximum(abs(a - b))`.
"""
maxad{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = Linfdist(a, b)

# MSD: mean squared deviation
"""
    maxad(a::AbstractArray{T}, b::AbstractArray{T})

Computes the mean squared deviation between a and b, i.e. mean(abs2(a - b)).
"""
msd{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = sqL2dist(a, b) / length(a)

# RMSD: root mean squared deviation
"""
    rmsd(a::AbstractArray{T}, b::AbstractArray{T}; normalize::Bool=false)

Computes the root mean squared deviation, i.e. sqrt(msd(a, b)). The keyword argument normalize is default to false. If it is set to true, the result is normalized by (maximum(a) - minimum(a).
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
    psnr(a::AbstractArray{T}, b::AbstractArray{T}, maxv::Real)

Computes peak signal-to-noise ratio, i.e. 10 * log10(maxv^2 / msd(a, b)).
"""
function psnr{T<:Real}(a::AbstractArray{T}, b::AbstractArray{T}, maxv::Real)
    20. * log10(maxv) - 10. * log10(msd(a, b))
end
