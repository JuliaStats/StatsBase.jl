# A variety of means

# Geometric mean
function geomean(a::RealArray)
    s = 0.0
    n = length(a)
    for i = 1 : n
        @inbounds s += log(a[i])
    end
    return exp(s / n)
end


# Harmonic mean
function harmmean(a::RealArray)
    s = 0.0
    n = length(a)
    for i in 1 : n
        @inbounds s += inv(a[i])
    end
    return n / s
end


# Trimmed mean

function trimmean(x::RealArray, p::FloatingPoint)
    n = length(x)
    n > 0 || error("x can not be empty.")
    rn = min(iround(n * p), n-1)
    if rn == n - 1
        return median(x)
    else
        sx = sort(x)
        nl = rn >> 1
        nh = (rn - nl)
        s = 0.0
        @inbounds for i = (1+nl) : (n-nh)
            s += x[i]
        end
        return s / (n - rn)
    end
end


# Weighted mean
function wmean{T<:Number,W<:Real}(v::AbstractArray{T}, w::AbstractArray{W}; wsum::W=NaN)
	# wsum is the pre-computed sum of weights

	n = length(v)
	length(w) == n || raise_dimerror()
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

mean{T<:Number}(v::AbstractArray{T}, w::WeightVec) = wmean(v, values(w); wsum=sum(w))

