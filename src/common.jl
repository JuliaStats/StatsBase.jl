# common utilities


## convenient type alias
#
#  These types signficantly reduces the need of using
#  type parameters in functions (which are often just
#  for the purpose of restricting the arrays to real)
#

typealias RealArray{T<:Real,N} AbstractArray{T,N}
typealias RealVector{T<:Real} AbstractArray{T,1}
typealias RealMatrix{T<:Real} AbstractArray{T,2}

typealias IntegerArray{T<:Integer,N} AbstractArray{T,N}
typealias IntegerVector{T<:Integer} AbstractArray{T,1}
typealias IntegerMatrix{T<:Integer} AbstractArray{T,2}

typealias RealFP Union(Float32, Float64)

## conversion from real to fp types

fptype(::Type{Float32}) = Float32
fptype(::Type{Float64}) = Float64
fptype(::Type{Bool}) = Float32
fptype(::Type{Int8}) = Float32
fptype(::Type{Uint8}) = Float32
fptype(::Type{Int16}) = Float32
fptype(::Type{Uint16}) = Float32
fptype(::Type{Int32}) = Float64
fptype(::Type{Uint32}) = Float64
fptype(::Type{Int64}) = Float64
fptype(::Type{Uint64}) = Float64
fptype(::Type{Int128}) = Float64
fptype(::Type{Uint128}) = Float64

## consistent way to raise error
#
# only change this, if we decide to change to way that
# dimension errors are raised
#

raise_dimerror(msg::String) = throw(DimensionMismatch(msg))
raise_dimerror() = raise_dimerror("")

# weight vector: 
#
#    a wrapper that indicator a vector representing a sequence of weights
#
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
