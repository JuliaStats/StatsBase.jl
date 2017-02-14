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

@compat RealArray{T<:Real,N} = AbstractArray{T,N}
@compat RealVector{T<:Real} = AbstractArray{T,1}
@compat RealMatrix{T<:Real} = AbstractArray{T,2}

@compat IntegerArray{T<:Integer,N} = AbstractArray{T,N}
@compat IntegerVector{T<:Integer} = AbstractArray{T,1}
@compat IntegerMatrix{T<:Integer} = AbstractArray{T,2}

@compat const RealFP = Union{Float32, Float64}

## conversion from real to fp types

@compat fptype{T<:Union{Float32,Bool,Int8,UInt8,Int16,UInt16}}(::Type{T}) = Float32
@compat fptype{T<:Union{Float64,Int32,UInt32,Int64,UInt64,Int128,UInt128}}(::Type{T}) = Float64
fptype(::Type{Complex64}) = Complex64
fptype(::Type{Complex128}) = Complex128
