
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


##### Weighted sum #####

## weighted sum over vectors

wsum(v::AbstractVector, w::AbstractVector) = dot(v, w)
wsum(v::AbstractArray, w::AbstractVector) = dot(vec(v), w)

# Note: the methods for BitArray and SparseMatrixCSC are to avoid ambiguities
Base.sum(v::BitArray, w::WeightVec) = wsum(v, values(w))
Base.sum(v::SparseMatrixCSC, w::WeightVec) = wsum(v, values(w))
Base.sum(v::AbstractArray, w::WeightVec) = dot(v, values(w))

# General Cartesian-based weighted sum across dimensions
@ngenerate N typeof(r) function wsum!{T,N,S,W<:Real}(r::AbstractArray{T,N}, v::AbstractArray{S,N},
                                                     w::AbstractVector{W}, dim::Int)
    1 <= dim <= N || error("dim = $dim not in range [1,$N]")
    for i = 1:N
        (i == dim && size(r, i) == 1 && size(v, i) == length(w)) || size(r, i) == size(v, i) || error(DimensionMismatch(""))
    end
    fill!(r, 0)
    weight = zero(W)
    @nloops N i v d->(if d == dim
                           weight = w[i_d]
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
function wsum!{W<:Real}(r::Union(Array, AbstractVector), v::Union(Array, AbstractMatrix), w::AbstractVector{W}, dim::Int)
    if dim == 1
        m = size(v, 1)
        n = div(length(v), m)
        (length(r) == n && length(w) == m) || throw(DimensionMismatch(""))
        At_mul_B!(vec(r), isa(v, AbstractMatrix) ? v : reshape(v, m, n), w)
    elseif dim == ndims(v)
        n = size(v, ndims(v))
        m = div(length(v), n)
        (length(r) == m && length(w) == n) || throw(DimensionMismatch(""))
        A_mul_B!(vec(r), isa(v, AbstractMatrix) ? v : reshape(v, m, n), w)
    else
        invoke(wsum!, (AbstractArray, AbstractArray, typeof(w), Int), r, v, w, dim)
    end
    r
end

Base.sum!{W<:Real}(r::AbstractArray, v::AbstractArray, w::WeightVec{W}, dim::Int) =
    wsum!(r, v, values(w), dim)

wsum{T<:Number,W<:Real}(v::AbstractArray{T}, w::AbstractVector{W}, dim::Int) =
    wsum!(Array(typeof(zero(T)*zero(W) + zero(T)*zero(W)), Base.reduced_dims(size(v), dim)), v, w, dim)

Base.sum{T<:Number,W<:Real}(v::AbstractArray{T}, w::WeightVec{W}, dim::Int) = wsum(v, values(w), dim)

