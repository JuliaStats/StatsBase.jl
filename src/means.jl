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

# Weighted means

mean{T<:Number}(v::AbstractArray{T}, w::WeightVec) = dot(v, values(w)) / sum(w)

function wmean{T<:Number}(v::AbstractArray{T}, w::AbstractArray)
    Base.depwarn("wmean is deprecated, use mean(v, weights(w)) instead.", :wmean)
    mean(v, weights(w))
end

function mean!{T<:Number,W<:Real}(r::AbstractVector, v::AbstractMatrix{T}, w::WeightVec{W}, dim::Int)
    m, n = size(v)
    if dim == 1
        (length(r) == n && length(w) == m) || throw(DimensionMismatch("Dimensions mismatch"))
        At_mul_B!(r, v, values(w))
    elseif dim == 2        
        (length(r) == m && length(w) == n) || throw(DimensionMismatch("Dimensions mismatch"))
        A_mul_B!(r, v, values(w))
    else
        error("Invalid value of dim.")
    end
    return scale!(r, inv(sum(w)))
end

function mean{T<:Number,W<:Real}(v::AbstractMatrix{T}, w::WeightVec{W}, dim::Int)
    R = typeof(float(one(T) * one(W)))
    m, n = size(v)
    dim == 1 ? reshape(mean!(Array(R, n), v, w, 1), 1, n) :
    dim == 2 ? reshape(mean!(Array(R, m), v, w, 2), m, 1) :
    error("Invalid value of dim.")
end






