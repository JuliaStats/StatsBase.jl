# Basic computational routines

import Base.BLAS: axpy!

## Inplace arithmetics

function negate!(x::AbstractArray)
    for i = 1:length(x)
        @inbounds x[i] = -x[i]
    end
    x
end

function add!(y::AbstractArray, x::Number)
    for i = 1:length(y)
        @inbounds y[i] += x
    end
    y
end

add!(y::AbstractArray, x::AbstractArray) = broadcast!(+, y, y, x)

function subtract!(y::AbstractArray, x::Number)
    for i = 1:length(y)
        @inbounds y[i] -= x
    end
    y
end

subtract!(y::AbstractArray, x::AbstractArray) = broadcast!(-, y, y, x)


## addscale!: y <- y + x * c

function generic_addscale!(y::AbstractArray, x::AbstractArray, c::Number)
    for i = 1:length(x)
        @inbounds y[i] += x[i] * c
    end
end

stride1(x::Vector) = 1
stride1(x::Array) = 1
stride1(x::StridedVector) = stride(x, 1)::Int

function blas_addscale!{T<:BlasFloat}(y::Union(Array{T},StridedVector{T}), 
                                      x::Union(Array{T},StridedVector{T}), c::T)
    n = length(x)
    n == 0 || axpy!(n, c, x, stride1(x), y, stride1(y))
end

function addscale!(y::AbstractArray, x::AbstractArray, c::Number)
    length(x) == length(y) || throw(DimensionMismatch("Inconsistent array lengths."))
    generic_addscale!(y, x, c)
    y
end


function addscale!{T<:BlasFloat}(y::Union(Array{T},StridedVector{T}), 
                                 x::Union(Array{T},StridedVector{T}), c::Number)
    n = length(x)
    length(y) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    if n < 256
        generic_addscale!(y, x, convert(T, c))
    else
        blas_addscale!(y, x, convert(T, c))
    end
    y
end

addscale(y::AbstractArray, x::AbstractArray, c::Number) = addscale!(copy(y), x, c)

