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

import Base.Cartesian: @ngenerate, @nloops, @nref
@ngenerate N typeof(r) function _wsum!{T<:Number,N,W<:Real}(r::Union(Array, AbstractVector), v::AbstractArray{T,N},
                                                            values::AbstractVector{W}, dim::Int)
    for i = 1:N
        i == dim || size(r, i) == size(v, i) || error(DimensionMismatch(""))
    end
    fill!(r, 0)
    weight = zero(W)
    @nloops N i v d->(if d == dim
                           weight = values[i_d]
                           j_d = 1
                       else
                           j_d = i_d
                       end) @inbounds (@nref N r j) += (@nref N v i)*weight
    r
end

function mean!{T<:Number,W<:Real}(r::Union(Array, AbstractVector), v::AbstractMatrix{T}, w::WeightVec{W}, dim::Int)
    m, n = size(v)
    if dim == 1
        (length(r) == n && length(w) == m) || throw(DimensionMismatch("Dimensions mismatch"))
        At_mul_B!(vec(r), v, values(w))
    elseif dim == 2        
        (length(r) == m && length(w) == n) || throw(DimensionMismatch("Dimensions mismatch"))
        A_mul_B!(vec(r), v, values(w))
    else
        error("Invalid value of dim.")
    end
    return scale!(r, inv(sum(w)))
end

function mean!{S,N,T<:Number,W<:Real}(r::AbstractArray{S,N}, v::AbstractArray{T,N}, w::WeightVec{W}, dim::Int)
    m, n = size(v)
    if dim == 1 && isa(v, Array) && isa(r, Array)
        m = size(v, 1)
        n = div(length(v), m)
        (length(r) == n && length(w) == m) || throw(DimensionMismatch("Dimensions mismatch"))
        At_mul_B!(vec(r), reshape(v, m, n), values(w))
    elseif dim == ndims(v) && isa(v, Array) && isa(r, Array)
        n = size(v, ndims(v))
        m = div(length(v), n)
        (length(r) == m && length(w) == n) || throw(DimensionMismatch("Dimensions mismatch"))
        A_mul_B!(vec(r), reshape(v, m, n), values(w))
    elseif 1 <= dim <= ndims(v)
        n = size(v, dim)
        (length(r) == div(length(v), n) && length(w) == n) || throw(DimensionMismatch("Dimensions mismatch"))
        _wsum!(r, v, values(w), dim)
    else
        error("Invalid value of dim.")
    end
    return scale!(r, inv(sum(w)))
end

mean{T<:Number,W<:Real}(v::AbstractArray{T}, w::WeightVec{W}, dim::Int) =
    mean!(Array(typeof(one(T) / one(W)), Base.reduced_dims(size(v), dim)), v, w, dim)
