
###### Weight vector #####

immutable WeightVec{W,Vec<:RealVector}
    values::Vec
    sum::W
end

"""
    WeightVec(vs, [wsum])

Construct a `WeightVec` with weight values `vs` and sum of weights `wsum`.
If omitted, `wsum` is computed.
"""
WeightVec{Vec<:RealVector,W<:Real}(vs::Vec,wsum::W) = WeightVec{W,Vec}(vs, wsum)
WeightVec(vs::RealVector) = WeightVec(vs, sum(vs))

"""
    weights(vs)

Construct a `WeightVec` from a given array.
"""
weights(vs::RealVector) = WeightVec(vs)
weights(vs::RealArray) = WeightVec(vec(vs))

eltype(wv::WeightVec) = eltype(wv.values)
length(wv::WeightVec) = length(wv.values)
values(wv::WeightVec) = wv.values
sum(wv::WeightVec) = wv.sum
isempty(wv::WeightVec) = isempty(wv.values)

Base.getindex(wv::WeightVec, i) = getindex(wv.values, i)


##### Weighted sum #####

## weighted sum over vectors

"""
    wsum(v, w::AbstractVector, [dim])

Compute the weighted sum of an array `v` with weights `w`, optionally over the dimension `dim`.
"""
wsum(v::AbstractVector, w::AbstractVector) = dot(v, w)
wsum(v::AbstractArray, w::AbstractVector) = dot(vec(v), w)

# Note: the methods for BitArray and SparseMatrixCSC are to avoid ambiguities
Base.sum(v::BitArray, w::WeightVec) = wsum(v, values(w))
Base.sum(v::SparseMatrixCSC, w::WeightVec) = wsum(v, values(w))
Base.sum(v::AbstractArray, w::WeightVec) = dot(v, values(w))

## wsum along dimension
#
#  Brief explanation of the algorithm:
#  ------------------------------------
#
#  1. _wsum! provides the core implementation, which assumes that
#     the dimensions of all input arguments are consistent, and no
#     dimension checking is performed therein.
#
#     wsum and wsum! perform argument checking and call _wsum!
#     internally.
#
#  2. _wsum! adopt a Cartesian based implementation for general
#     sub types of AbstractArray. Particularly, a faster routine
#     that keeps a local accumulator will be used when dim = 1.
#
#     The internal function that implements this is _wsum_general!
#
#  3. _wsum! is specialized for following cases:
#     (a) A is a vector: we invoke the vector version wsum above.
#         The internal function that implements this is _wsum1!
#
#     (b) A is a dense matrix with eltype <: BlasReal: we call gemv!
#         The internal function that implements this is _wsum2_blas!
#
#     (c) A is a contiguous array with eltype <: BlasReal:
#         dim == 1: treat A like a matrix of size (d1, d2 x ... x dN)
#         dim == N: treat A like a matrix of size (d1 x ... x d(N-1), dN)
#         otherwise: decompose A into multiple pages, and apply _wsum2!
#         for each
#
#     (d) A is a general dense array with eltype <: BlasReal:
#         dim <= 2: delegate to (a) and (b)
#         otherwise, decompose A into multiple pages
#

function _wsum1!(R::AbstractArray, A::AbstractVector, w::AbstractVector, init::Bool)
    r = wsum(A, w)
    if init
        R[1] = r
    else
        R[1] += r
    end
    return R
end

function _wsum2_blas!{T<:BlasReal}(R::StridedVector{T}, A::StridedMatrix{T}, w::StridedVector{T}, dim::Int, init::Bool)
    beta = ifelse(init, zero(T), one(T))
    trans = dim == 1 ? 'T' : 'N'
    BLAS.gemv!(trans, one(T), A, w, beta, R)
    return R
end

function _wsumN!{T<:BlasReal,N}(R::ContiguousArray{T}, A::ContiguousArray{T,N}, w::StridedVector{T}, dim::Int, init::Bool)
    if dim == 1
        m = size(A, 1)
        n = div(length(A), m)
        _wsum2_blas!(view(R,:), reshape_view(A, (m, n)), w, 1, init)
    elseif dim == N
        n = size(A, N)
        m = div(length(A), n)
        _wsum2_blas!(view(R,:), reshape_view(A, (m, n)), w, 2, init)
    else # 1 < dim < N
        m = 1
        for i = 1:dim-1; m *= size(A, i); end
        n = size(A, dim)
        k = 1
        for i = dim+1:N; k *= size(A, i); end
        Av = reshape_view(A, (m, n, k))
        Rv = reshape_view(R, (m, k))
        for i = 1:k
            _wsum2_blas!(view(Rv,:,i), view(Av,:,:,i), w, 2, init)
        end
    end
    return R
end

function _wsumN!{T<:BlasReal,N}(R::ContiguousArray{T}, A::DenseArray{T,N}, w::StridedVector{T}, dim::Int, init::Bool)
    @assert N >= 3
    if dim <= 2
        m = size(A, 1)
        n = size(A, 2)
        npages = 1
        for i = 3:N
            npages *= size(A, i)
        end
        rlen = ifelse(dim == 1, n, m)
        Rv = reshape_view(R, (rlen, npages))
        for i = 1:npages
            _wsum2_blas!(view(Rv,:,i), view(A,:,:,i), w, dim, init)
        end
    else
        _wsum_general!(R, @functorize(identity), A, w, dim, init)
    end
    return R
end

# General Cartesian-based weighted sum across dimensions
@ngenerate N typeof(R) function _wsum_general!{T,RT,WT,N}(R::AbstractArray{RT}, f::supertype(typeof(@functorize(abs))),
                                                          A::AbstractArray{T,N}, w::AbstractVector{WT}, dim::Int, init::Bool)
    init && fill!(R, zero(RT))
    wi = zero(WT)
    if dim == 1
        @nextract N sizeR d->size(R,d)
        sizA1 = size(A, 1)
        @nloops N i d->(d>1? (1:size(A,d)) : (1:1)) d->(j_d = sizeR_d==1 ? 1 : i_d) begin
            @inbounds r = (@nref N R j)
            for i_1 = 1:sizA1
                @inbounds r += f(@nref N A i) * w[i_1]
            end
            @inbounds (@nref N R j) = r
        end
    else
        @nloops N i A d->(if d == dim
                               wi = w[i_d]
                               j_d = 1
                           else
                               j_d = i_d
                           end) @inbounds (@nref N R j) += f(@nref N A i) * wi
    end
    return R
end

@ngenerate N typeof(R) function _wsum_centralize!{T,RT,WT,N}(R::AbstractArray{RT}, f::supertype(typeof(@functorize(abs))),
                                                             A::AbstractArray{T,N}, w::AbstractVector{WT}, means,
                                                             dim::Int, init::Bool)
    init && fill!(R, zero(RT))
    wi = zero(WT)
    if dim == 1
        @nextract N sizeR d->size(R,d)
        sizA1 = size(A, 1)
        @nloops N i d->(d>1? (1:size(A,d)) : (1:1)) d->(j_d = sizeR_d==1 ? 1 : i_d) begin
            @inbounds r = (@nref N R j)
            @inbounds m = (@nref N means j)
            for i_1 = 1:sizA1
                @inbounds r += f((@nref N A i) - m) * w[i_1]
            end
            @inbounds (@nref N R j) = r
        end
    else
        @nloops N i A d->(if d == dim
                               wi = w[i_d]
                               j_d = 1
                           else
                               j_d = i_d
                           end) @inbounds (@nref N R j) += f((@nref N A i) - (@nref N means j)) * wi
    end
    return R
end


# N = 1
_wsum!{T<:BlasReal}(R::ContiguousArray{T}, A::DenseArray{T,1}, w::StridedVector{T}, dim::Int, init::Bool) =
    _wsum1!(R, A, w, init)

# N = 2
_wsum!{T<:BlasReal}(R::ContiguousArray{T}, A::DenseArray{T,2}, w::StridedVector{T}, dim::Int, init::Bool) =
    (_wsum2_blas!(view(R,:), A, w, dim, init); R)

# N >= 3
_wsum!{T<:BlasReal,N}(R::ContiguousArray{T}, A::DenseArray{T,N}, w::StridedVector{T}, dim::Int, init::Bool) =
    _wsumN!(R, A, w, dim, init)

_wsum!(R::AbstractArray, A::AbstractArray, w::AbstractVector, dim::Int, init::Bool) = _wsum_general!(R, @functorize(identity), A, w, dim, init)

## wsum! and wsum

wsumtype{T,W}(::Type{T}, ::Type{W}) = typeof(zero(T) * zero(W) + zero(T) * zero(W))
wsumtype{T<:BlasReal}(::Type{T}, ::Type{T}) = T


"""
    wsum!(R, A, w, dim; init=true)

Compute the weighted sum of `A` with weights `w` over the dimension `dim` and store
the result in `R`. If `init=false`, the sum is added to `R` rather than starting
from zero.
"""
function wsum!{T,N}(R::AbstractArray, A::AbstractArray{T,N}, w::AbstractVector, dim::Int; init::Bool=true)
    1 <= dim <= N || error("dim should be within [1, $N]")
    ndims(R) <= N || error("ndims(R) should not exceed $N")
    length(w) == size(A,dim) || throw(DimensionMismatch("Inconsistent array dimension."))
    # TODO: more careful examination of R's size
    _wsum!(R, A, w, dim, init)
end

function wsum{T<:Number,W<:Real}(A::AbstractArray{T}, w::AbstractVector{W}, dim::Int)
    length(w) == size(A,dim) || throw(DimensionMismatch("Inconsistent array dimension."))
    _wsum!(Array(wsumtype(T,W), Base.reduced_dims(size(A), dim)), A, w, dim, true)
end

# extended sum! and wsum

Base.sum!{W<:Real}(R::AbstractArray, A::AbstractArray, w::WeightVec{W}, dim::Int; init::Bool=true) =
    wsum!(R, A, values(w), dim; init=init)

Base.sum{T<:Number,W<:Real}(A::AbstractArray{T}, w::WeightVec{W}, dim::Int) = wsum(A, values(w), dim)


###### Weighted means #####

"""
    wmean(v, w::AbstractVector)

Compute the weighted mean of an array `v` with weights `w`.
"""
function wmean{T<:Number}(v::AbstractArray{T}, w::AbstractVector)
    Base.depwarn("wmean is deprecated, use mean(v, weights(w)) instead.", :wmean)
    mean(v, weights(w))
end

Base.mean(v::AbstractArray, w::WeightVec) = sum(v, w) / sum(w)

Base.mean!(R::AbstractArray, A::AbstractArray, w::WeightVec, dim::Int) =
    scale!(Base.sum!(R, A, w, dim), inv(sum(w)))

wmeantype{T,W}(::Type{T}, ::Type{W}) = typeof((zero(T)*zero(W) + zero(T)*zero(W)) / one(W))
wmeantype{T<:BlasReal}(::Type{T}, ::Type{T}) = T

Base.mean{T<:Number,W<:Real}(A::AbstractArray{T}, w::WeightVec{W}, dim::Int) =
    mean!(Array(wmeantype(T, W), Base.reduced_dims(size(A), dim)), A, w, dim)


###### Weighted median #####

function Base.median{W<:Real}(v::RealVector, w::WeightVec{W})
    isempty(v) && error("median of an empty array is undefined")
    if length(v) != length(w)
        error("data and weight vectors must be the same size")
    end
    @inbounds for x in w.values
        isnan(x) && error("weight vector cannot contain NaN entries")
    end
    @inbounds for x in v
        isnan(x) && return x
    end
    mask = w.values .!= 0
    if !any(mask)
        error("all weights are zero")
    end
    if all(w.values .<= 0)
        error("no positive weights found")
    end
    v = v[mask]
    wt = w[mask]
    midpoint = w.sum / 2
    maxval, maxind = findmax(wt)
    if maxval > midpoint
        v[maxind]
    else
        permute = sortperm(v)
        cumulative_weight = zero(eltype(wt))
        i = 0
        for (i, p) in enumerate(permute)
            if cumulative_weight == midpoint
                i += 1
                break
            elseif cumulative_weight > midpoint
                cumulative_weight -= wt[p]
                break
            end
            cumulative_weight += wt[p]
        end
        if cumulative_weight == midpoint
            middle(v[permute[i-2]], v[permute[i-1]])
        else
            middle(v[permute[i-1]])
        end
    end
end


"""
    wmedian(v, w)

Compute the weighted median of an array `v` with weights `w`, given as either a
vector or `WeightVec`.
"""
wmedian(v::RealVector, w::RealVector) = median(v, weights(w))
wmedian{W<:Real}(v::RealVector, w::WeightVec{W}) = median(v, w)

###### Weighted quantile #####

# http://stats.stackexchange.com/questions/13169/defining-quantiles-over-a-weighted-sample
# In the non weighted version, we compute a vector of index h(N, p)
# and take interpolation between floor and ceil of this index
# Here there is a supplementary function from index to weighted index k -> Sk

"""
    quantile(v, w::WeightVec, p)

Compute `p`th quantile(s) of `v` with weights `w`.
"""
function quantile{V, W <: Real}(v::RealVector{V}, w::WeightVec{W}, p::RealVector)

    # checks
    isempty(v) && error("quantile of an empty array is undefined")
    isempty(p) && throw(ArgumentError("empty quantile array"))

    w.sum == 0 && error("weight vector cannot sum to zero")
    length(v) == length(w) || error("data and weight vectors must be the same size, got $(length(v)) and $(length(w))")
    for x in w.values
        isnan(x) && error("weight vector cannot contain NaN entries")
        x < 0 && error("weight vector cannot contain negative entries")
    end

    # full sort
    vw = sort!(collect(zip(v, w.values)))

    wsum = w.sum

    # prepare percentiles
    ppermute = sortperm(p)
    p = p[ppermute]
    p = bound_quantiles(p)

    # prepare out vector
    N = length(vw)
    out = Array(typeof(zero(V)/1), length(p))
    fill!(out, vw[end][1])

    # start looping on quantiles
    cumulative_weight, Sk, Skold =  zero(W), zero(W), zero(W)
    vk, vkold = zero(V), zero(V)
    k = 1
    for i in 1:length(p)
        h = p[i] * (N - 1) * wsum
        if h == 0
            # happens when N or p or wsum equal zero
            out[ppermute[i]] = vw[1][1]
        else
            while Sk <= h
                # happens in particular when k == 1
                vk, wk = vw[k]
                cumulative_weight += wk
                if k >= N
                    # out was initialized with maximum v
                    return out
                end
                k += 1
                Skold, vkold = Sk, vk
                vk, wk = vw[k]
                Sk = (k - 1) * wk + (N - 1) * cumulative_weight
            end
            # in particular, Sk is different from Skold
            g = (h - Skold) / (Sk - Skold)
            out[ppermute[i]] = vkold + g * (vk - vkold)
        end
    end
    return out
end

# similarly to statistics.jl in Base
function bound_quantiles{T <: Real}(qs::AbstractVector{T})
    epsilon = 100 * eps()
    if (any(qs .< -epsilon) || any(qs .> 1+epsilon))
        throw(ArgumentError("quantiles out of [0,1] range"))
    end
    T[min(one(T), max(zero(T), q)) for q = qs]
end

quantile{W <: Real}(v::RealVector, w::WeightVec{W}, p::Number) = quantile(v, w, [p])[1]


"""
    wquantile(v, w, p)

Compute the `p`th quantile(s) of `v` with weights `w`, given as either a vector
or a `WeightVec`.
"""
wquantile{W <: Real}(v::RealVector, w::WeightVec{W}, p::RealVector) = quantile(v, w, p)
wquantile{W <: Real}(v::RealVector, w::WeightVec{W}, p::Number) = quantile(v, w, [p])[1]
wquantile(v::RealVector, w::RealVector, p::RealVector) = quantile(v, weights(w), p)
wquantile(v::RealVector, w::RealVector, p::Number) = quantile(v, weights(w), [p])[1]
