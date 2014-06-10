# common utilities

## convenient type alias
#
#  These types signficantly reduces the need of using
#  type parameters in functions (which are often just
#  for the purpose of restricting the arrays to real)
#
# These could be removed when the Base supports
# covariant type notation, i.e. AbstractVector{<:Real}
#

typealias RealArray{T<:Real,N} AbstractArray{T,N}
typealias RealVector{T<:Real} AbstractArray{T,1}
typealias RealMatrix{T<:Real} AbstractArray{T,2}

typealias IntegerArray{T<:Integer,N} AbstractArray{T,N}
typealias IntegerVector{T<:Integer} AbstractArray{T,1}
typealias IntegerMatrix{T<:Integer} AbstractArray{T,2}

typealias RealFP Union(Float32, Float64)

## conversion from real to fp types

fptype{T<:Union(Float32,Bool,Int8,Uint8,Int16,Uint16)}(::Type{T}) = Float32
fptype{T<:Union(Float64,Int64,Uint64,Int128,Uint128)}(::Type{T}) = Float64
fptype(::Type{Complex64}) = Complex64
fptype(::Type{Complex128}) = Complex128
