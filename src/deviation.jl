# Computing deviation in a variety of ways

## count the number of equal/non-equal pairs

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
L2dist{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = sqrt(sqL2dist(a, b))

# L1 distance
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
function gkldiv{T<:FloatingPoint}(a::AbstractArray{T}, b::AbstractArray{T})
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
meanad{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = L1dist(a, b) / length(a)

# MaxAD: maximum absolute deviation
maxad{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = Linfdist(a, b)

# MSD: mean squared deviation
msd{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}) = sqL2dist(a, b) / length(a)

# RMSD: root mean squared deviation
function rmsd{T<:Number}(a::AbstractArray{T}, b::AbstractArray{T}; normalize::Bool=false) 
	v = sqrt(msd(a, b))
	if normalize
		amin, amax = extrema(a)
		v /= (amax - amin)
	end
	return v
end

# PSNR: peak signal-to-noise ratio
function psnr{T<:Real}(a::AbstractArray{T}, b::AbstractArray{T}, maxv::Real)
    20. * log10(maxv) - 10. * log10(msd(a, b))
end
