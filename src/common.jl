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

# A convenient typealias for deprecating default corrected Bool
const DepBool = Union{Bool, Nothing}

function depcheck(fname::Symbol, varname::Symbol, b::DepBool)
    if b === nothing
        msg = "$fname will default to $varname=true in the future. Use $varname=false for previous behaviour."
        Base.depwarn(msg, fname)
        false
    else
        b
    end
end
