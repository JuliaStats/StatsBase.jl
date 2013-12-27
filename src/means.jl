# A variety of means

# Geometric mean
function gmean(a::RealArray)
    s = 0.0
    n = length(a)
    for i = 1 : n
        @inbounds s += log(a[i])
    end
    return exp(s / n)
end

# Harmonic mean
function hmean(a::RealArray)
    s = 0.0
    n = length(a)
    for i in 1 : n
        @inbounds s += inv(a[i])
    end
    return n / s
end

# Weighted mean
function wmean{T<:Number,W<:Real}(v::AbstractArray{T}, w::AbstractArray{W}; wsum::W=NaN)
	# wsum is the pre-computed sum of weights

	n = length(v)
	length(w) == n || throw(ArgumentError("Inconsistent array lengths."))
    sv = zero(T)

    if isnan(wsum)
    	sw = zero(W)
    	for i = 1 : n
    		@inbounds wi = w[i]
        	@inbounds sv += v[i] * wi
        	sw += wi
    	end
    else
    	for i = 1 : n
    		@inbounds sv += v[i] * w[i]
    	end
    	sw = wsum
    end

    return sv / sw
end

# Weighted mean (fast version for contiguous arrays using BLAS)
function wmean{T<:BlasReal}(v::Array{T}, w::Array{T}; wsum::T=NaN) 
	sw = (isnan(wsum) ? sum(w) : wsum)
	dot(v, w) / sw
end

wmean{T<:Number}(v::AbstractArray{T}, w::WeightVec) = wmean(v, values(w), sum(w))

