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

function wmean{T<:Number}(v::AbstractArray{T}, w::AbstractArray)
    Base.depwarn("wmean is deprecated, use mean(v, weights(w)) instead.", :wmean)
    mean(v, weights(w))
end

Base.mean(v::AbstractArray, w::WeightVec) = sum(v, w) / sum(w)

Base.mean!(r::AbstractArray, v::AbstractArray, w::WeightVec, dim::Int) =
    scale!(Base.sum!(r, v, w, dim), inv(sum(w)))

Base.mean{T<:Number,W<:Real}(v::AbstractArray{T}, w::WeightVec{W}, dim::Int) =
    mean!(Array(typeof((zero(T)*zero(W) + zero(T)*zero(W)) / one(W)), Base.reduced_dims(size(v), dim)), v, w, dim)

