# Basic computational routines

import Base.BLAS: axpy!
import Base.LinAlg: vecnorm1, vecnorm2, vecnormInf

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


## some statistics computation

sumzero{T<:Number}(::Type{T}) = float(real(zero(T) + zero(T)))

sumabs{T<:Number}(x::AbstractArray{T}) = isempty(x) ? sumzero(T) : vecnorm1(x)
maxabs{T<:Number}(x::AbstractArray{T}) = isempty(x) ? sumzero(T) : vecnormInf(x)

function generic_sumabs2(x)
    s = start(x)
    (v, s) = next(x, s)
    av = float(abs2(v))
    T = typeof(av)
    sum::promote_type(Float64, T) = av
    while !done(x, s)
        (v, s) = next(x, s)
        sum += abs2(v)
    end
    return convert(T, sum)
end

sumabs2{T<:Number}(x::AbstractArray{T}) = isempty(x) ? sumzero(T) : generic_sumabs2(x)

function sumabs2{T<:BlasFloat}(x::Union(Array{T},StridedVector{T}))
    isempty(x) && return sumzero(T)
    if length(x) < 128
        generic_sumabs2(x)
    else
        incx = stride1(x)::Int
        BLAS.dot(length(x), x, incx, x, incx)
    end
end


function sumabsdiff{T<:Number}(x::AbstractArray{T}, y::AbstractArray{T})
    n = length(x)
    length(y) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    s = sumzero(T)
    for i = 1:n
        @inbounds s += abs(x[i] - y[i])
    end
    return s
end

function sumabsdiff{T<:Number}(x::AbstractArray{T}, y::Number)
    n = length(x)
    s = sumzero(T)
    yv = convert(T, y)::T
    for i = 1:n
        @inbounds s += abs(x[i] - yv)
    end
    return s
end

sumabsdiff{T<:Number}(x::Number, y::AbstractArray{T}) = sumabsdiff(y, x)


function sumabs2diff{T<:Number}(x::AbstractArray{T}, y::AbstractArray{T})
    n = length(x)
    length(y) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    s = sumzero(T)
    for i = 1:n
        @inbounds s += abs2(x[i] - y[i])
    end
    return s
end

function sumabs2diff{T<:Number}(x::AbstractArray{T}, y::Number)
    n = length(x)
    s = sumzero(T)
    yv = convert(T, y)::T
    for i = 1:n
        @inbounds s += abs2(x[i] - yv)
    end
    return s
end

sumabs2diff{T<:Number}(x::Number, y::AbstractArray{T}) = sumabs2diff(y, x)


function maxabsdiff{T<:Number}(x::AbstractArray{T}, y::AbstractArray{T})
    n = length(x)
    length(y) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    s = zero(T)
    for i = 1:n
        @inbounds av = abs(x[i] - y[i])
        if av > s  # this can handle NaN correctly
            s = av
        end
    end
    return s
end

function maxabsdiff{T<:Number}(x::AbstractArray{T}, y::Number)
    n = length(x)
    s = zero(T)
    yv = convert(T, y)::T
    for i = 1:n
        @inbounds av = abs(x[i] - yv)
        if av > s
            s = av
        end
    end
    return s
end

maxabsdiff{T<:Number}(x::Number, y::AbstractArray{T}) = maxabsdiff(y, x)


# scalar functions

xlogx(x::FloatingPoint) = x > zero(x) ? x * log(x) : zero(x)
xlogx(x::Real) = xlogx(float(x))

xlogy{T<:FloatingPoint}(x::T, y::T) = x > zero(T) ? x * log(y) : zero(x)
xlogy{T<:Real}(x::T, y::T) = xlogy(float(x), float(y))
xlogy(x::Real, y::Real) = xlogy(promote(x, y)...)

logistic(x::FloatingPoint) = one(x) / (one(x) + exp(-x))
logistic(x::Real) = logistic(float(x))

logit(x::FloatingPoint) = log(x / (one(x) - x))
logit(x::Real) = logit(float(x))

softplus(x::FloatingPoint) = x <= 0 ? log1p(exp(x)) : x + log1p(exp(-x))
softplus(x::Real) = softplus(float(x))

invsoftplus(x::FloatingPoint) = log(exp(x) - one(x))
invsoftplus(x::Real) = invsoftplus(float(x))

@vectorize_1arg Real xlogx
@vectorize_2arg Real xlogy
@vectorize_1arg Real logistic
@vectorize_1arg Real logit
@vectorize_1arg Real softplus
@vectorize_1arg Real invsoftplus

## logsumexp

logsumexp{T<:FloatingPoint}(x::T, y::T) = x > y ? x + log1p(exp(y - x)) : y + log1p(exp(x - y))
logsumexp{T<:Real}(x::T, y::T) = logsumexp(float(x), float(y))
logsumexp(x::Real, y::Real) = logsumexp(promote(x, y)...)

function logsumexp{T<:Real}(x::AbstractArray{T})
    isempty(x) && return -Inf
    u = maximum(x)
    s = 0.
    for i = 1:length(x)
        @inbounds s += exp(x[i] - u)
    end
    log(s) + u
end

## softmax

function softmax!{T<:Real}(r::AbstractArray, x::AbstractArray{T})
    n = length(x)
    length(r) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    u = maximum(x)
    s = 0.
    for i = 1:n
        @inbounds s += (r[i] = exp(x[i] - u))
    end
    for i = 1:n
        r[i] /= s
    end
    r
end

softmax!{T<:FloatingPoint}(x::AbstractArray{T}) = softmax!(x, x)
softmax{T<:Real}(x::AbstractArray{T}) = softmax!(Array(Float64, size(x)), x)

