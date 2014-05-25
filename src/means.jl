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

sum(v::BitArray, w::WeightVec) = dot(vec(v), values(w))
sum(v::SparseMatrixCSC, w::WeightVec) = dot(vec(v), values(w))
sum(v::AbstractArray, w::WeightVec) = dot(vec(v), values(w))
mean(v::AbstractArray, w::WeightVec) = sum(v, w) / sum(w)

function wmean{T<:Number}(v::AbstractArray{T}, w::AbstractArray)
    Base.depwarn("wmean is deprecated, use mean(v, weights(w)) instead.", :wmean)
    mean(v, weights(w))
end

# General Cartesian-based weighted sum across dimensions
import Base.Cartesian: @ngenerate, @nloops, @nref
@ngenerate N typeof(r) function Base.sum!{T,N,S,W<:Real}(r::AbstractArray{T,N}, v::AbstractArray{S,N},
                                                         w::WeightVec{W}, dim::Int)
    vals = values(w)
    1 <= dim <= N || error("dim = $dim not in range (1,$N)")
    for i = 1:N
        (i == dim && size(r, i) == 1 && size(v, i) == length(vals)) || size(r, i) == size(v, i) || error(DimensionMismatch(""))
    end
    fill!(r, 0)
    weight = zero(W)
    @nloops N i v d->(if d == dim
                           weight = vals[i_d]
                           j_d = 1
                       else
                           j_d = i_d
                       end) @inbounds (@nref N r j) += (@nref N v i)*weight
    r
end

# Weighted sum via `A_mul_B!`/`At_mul_B!` for first and last
# dimensions of compatible arrays. `vec` and `reshape` are only
# guaranteed not to make a copy for Arrays, so only supports Arrays if
# these calls may be necessary.
function Base.sum!{W<:Real}(r::Union(Array, AbstractVector), v::Union(Array, AbstractMatrix), w::WeightVec{W}, dim::Int)
    if dim == 1
        m = size(v, 1)
        n = div(length(v), m)
        (length(r) == n && length(w) == m) || throw(DimensionMismatch(""))
        At_mul_B!(vec(r), isa(v, AbstractMatrix) ? v : reshape(v, m, n), values(w))
    elseif dim == ndims(v)
        n = size(v, ndims(v))
        m = div(length(v), n)
        (length(r) == m && length(w) == n) || throw(DimensionMismatch(""))
        A_mul_B!(vec(r), isa(v, AbstractMatrix) ? v : reshape(v, m, n), values(w))
    else
        invoke(sum!, (AbstractArray, AbstractArray, typeof(w), Int), r, v, w, dim)
    end
    r
end

Base.mean!(r::AbstractArray, v::AbstractArray, w::WeightVec, dim::Int) =
    scale!(Base.sum!(r, v, w, dim), inv(sum(w)))

Base.sum{T<:Number,W<:Real}(v::AbstractArray{T}, w::WeightVec{W}, dim::Int) =
    sum!(Array(typeof(zero(T)*zero(W) + zero(T)*zero(W)), Base.reduced_dims(size(v), dim)), v, w, dim)

Base.mean{T<:Number,W<:Real}(v::AbstractArray{T}, w::WeightVec{W}, dim::Int) =
    mean!(Array(typeof((zero(T)*zero(W) + zero(T)*zero(W)) / one(W)), Base.reduced_dims(size(v), dim)), v, w, dim)

