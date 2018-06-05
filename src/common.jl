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

const RealArray{T<:Real,N} = AbstractArray{T,N}
const RealVector{T<:Real} = AbstractArray{T,1}
const RealMatrix{T<:Real} = AbstractArray{T,2}

const IntegerArray{T<:Integer,N} = AbstractArray{T,N}
const IntegerVector{T<:Integer} = AbstractArray{T,1}
const IntegerMatrix{T<:Integer} = AbstractArray{T,2}

const RealFP = Union{Float32, Float64}

## conversion from real to fp types

fptype(::Type{T}) where {T<:Union{Float32,Bool,Int8,UInt8,Int16,UInt16}} = Float32
fptype(::Type{T}) where {T<:Union{Float64,Int32,UInt32,Int64,UInt64,Int128,UInt128}} = Float64
fptype(::Type{Complex{Float32}}) = Complex{Float32}
fptype(::Type{Complex{Float64}}) = Complex{Float64}

# A convenient typealias for deprecating default corrected Bool
const DepBool = Union{Bool, Nothing}

function depcheck(fname::Symbol, b::DepBool)
    if b == nothing
        msg = "$fname will default to corrected=true in the future. Use corrected=false for previous behaviour."
        Base.depwarn(msg, fname)
        false
    else
        b
    end
end
