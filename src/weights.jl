
###### Weight vector #####

immutable WeightVec{W<:Real,Vec<:RealVector}
    values::Vec
    sum::W
end

WeightVec{Vec<:RealVector,W<:Real}(vs::Vec,wsum::W) = WeightVec{W,Vec}(vs, wsum)
WeightVec(vs::RealVector) = WeightVec(vs, sum(vs))

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
        _wsum_general!(R, A, w, dim, init)
    end
    return R
end

# General Cartesian-based weighted sum across dimensions
@ngenerate N typeof(R) function _wsum_general!{T,RT,WT,N}(R::AbstractArray{RT}, 
                                                          A::AbstractArray{T,N}, w::AbstractVector{WT}, dim::Int, init::Bool)
    init && fill!(R, zero(RT))
    wi = zero(WT)
    if dim == 1
        @nextract N sizeR d->size(R,d)
        sizA1 = size(A, 1)
        @nloops N i d->(d>1? (1:size(A,d)) : (1:1)) d->(j_d = sizeR_d==1 ? 1 : i_d) begin
            @inbounds r = (@nref N R j)
            for i_1 = 1:sizA1
                @inbounds r += (@nref N A i) * w[i_1]
            end
            @inbounds (@nref N R j) = r
        end 
    else
        @nloops N i A d->(if d == dim
                               wi = w[i_d]
                               j_d = 1
                           else
                               j_d = i_d
                           end) @inbounds (@nref N R j) += (@nref N A i) * wi
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

_wsum!(R::AbstractArray, A::AbstractArray, w::AbstractVector, dim::Int, init::Bool) = _wsum_general!(R, A, w, dim, init)

## wsum! and wsum

wsumtype{T,W}(::Type{T}, ::Type{W}) = typeof(zero(T) * zero(W) + zero(T) * zero(W))
wsumtype{T<:BlasReal}(::Type{T}, ::Type{T}) = T

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


