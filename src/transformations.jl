### Transformations

abstract type DataTransform end

# apply the transform
transform!(t::S, x::AbstractArray{T,1}) where {T<:Real, S<:DataTransform} = transform!(x, t, x)
transform!(t::S, x::AbstractArray{T,2}) where {T<:Real, S<:DataTransform} = transform!(x, t, x)

transform(t::S, x::AbstractArray{T,1}) where {T<:Real, S<:DataTransform} = transform!(similar(x), t, x)
transform(t::S, x::AbstractArray{T,2}) where {T<:Real, S<:DataTransform} = transform!(similar(x), t, x)

# reconstruct the original data from transformed values
reconstruct!(t::S, x::AbstractArray{T,1}) where {T<:Real, S<:DataTransform} = reconstruct!(x, t, x)
reconstruct!(t::S, x::AbstractArray{T,2}) where {T<:Real, S<:DataTransform} = reconstruct!(x, t, x)

reconstruct(t::S, y::AbstractArray{T,1}) where {T<:Real, S<:DataTransform} = reconstruct!(similar(y), t, y)
reconstruct(t::S, y::AbstractArray{T,2}) where {T<:Real, S<:DataTransform} = reconstruct!(similar(y), t, y)

# Z-score transformation
struct ZScoreTransform{T<:Real} <: DataTransform
    dim::Int
    mean::Vector{T}
    scale::Vector{T}

    function ZScoreTransform(d::Int, m::Vector{T}, s::Vector{T}) where T
        lenm = length(m)
        lens = length(s)
        lenm == d || lenm == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        lens == d || lens == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        new{T}(d, m, s)
    end
end

function Base.getproperty(t::ZScoreTransform, p::Symbol)
    if p == :indim || p == :outdim
        return t.dim
    else
        return getfield(t, p)
    end
end

# fit a z-score transform
function fit(::Type{ZScoreTransform}, X::AbstractArray{T,2}; center::Bool=true, scale::Bool=true) where T<:Real
    d, n = size(X)
    n >= 2 || error("X must contain at least two columns.")

    m, s = mean_and_std(X, 2)

    return ZScoreTransform(d, (center ? vec(m) : zeros(T, 0)),
                              (scale ? vec(s) : zeros(T, 0)))
end

function transform!(y::AbstractVecOrMat{YT}, t::ZScoreTransform, x::AbstractVecOrMat{XT}) where {YT<:Real,XT<:Real}
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(y,2)
    size(x,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if x !== y
                copyto!(y, x)
            end
        else
            broadcast!(/, y, x, s)
        end
    else
        if isempty(s)
            broadcast!(-, y, x, m)
        else
            broadcast!((x,m,s)->(x-m)/s, y, x, m, s)
        end
    end
    return y
end

function reconstruct!(x::AbstractVecOrMat{YT}, t::ZScoreTransform, y::AbstractVecOrMat{XT}) where {YT<:Real,XT<:Real}
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(y,2)
    size(x,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if y !== x
                copyto!(x, y)
            end
        else
            broadcast!(*, x, y, s)
        end
    else
        if isempty(s)
            broadcast!(+, x, y, m)
        else
            broadcast!((y,m,s)->y*s+m, x, y, m, s)
        end
    end
    return x
end

# Unit transformation
struct UnitRangeTransform{T<:Real}  <: DataTransform
    dim::Int
    unit::Bool
    min::Vector{T}
    scale::Vector{T}

    function UnitRangeTransform(d::Int, unit::Bool, min::Vector{T}, max::Vector{T}) where T
        lenmin = length(min)
        lenmax = length(max)
        lenmin == d || lenmin == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        lenmax == d || lenmax == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        new{T}(d, unit, min, max)
    end
end

function Base.getproperty(t::UnitRangeTransform, p::Symbol)
    if p == :indim || p == :outdim
        return t.dim
    else
        return getfield(t, p)
    end
end

# fit a unit transform
function fit(::Type{UnitRangeTransform}, X::AbstractArray{T,2}; unit::Bool=true) where T<:Real
    d, n = size(X)

    tmin = zeros(T, d)
    tmax = zeros(T, d)
    copyto!(tmin, X[:, 1])
    copyto!(tmax, X[:, 1])
    for j = 2:n
        @inbounds for i = 1:d
            if X[i, j] < tmin[i]
                tmin[i] = X[i, j]
            elseif X[i, j] > tmax[i]
                tmax[i] = X[i, j]
            end
        end
    end
    for i = 1:d
        @inbounds tmax[i] = 1.0 / (tmax[i] - tmin[i])
    end
    return UnitRangeTransform(d, unit, tmin, tmax)
end

function transform!(y::AbstractVecOrMat{YT}, t::UnitRangeTransform, x::AbstractVecOrMat{XT}) where {YT<:Real,XT<:Real}
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(x,2)
    size(y,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        broadcast!((x,s,m) -> (x-m)*s, y, x, tscale, tmin)
    else
        broadcast!(*, y, x, tscale)
    end
    return y
end

function reconstruct!(x::AbstractVecOrMat{XT}, t::UnitRangeTransform, y::AbstractVecOrMat{YT}) where {YT<:Real,XT<:Real}
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(y,2)
    size(x,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        broadcast!((y,s,m) -> y/s + m, x, y, tscale, tmin)
    else
        broadcast!(/, x, y, tscale)
    end
    return x
end

"""
    standardize(DT, X; kwargs...)

Return a row-standardized matrix `X` using `DT` transformation.

# Example
```julia
julia> standardize(ZScoreTransform, [0.0 -0.5 0.5; 0.0 1.0 2.0])
2×3 Array{Float64,2}:
  0.0  -1.0  1.0
 -1.0   0.0  1.0

julia> standardize(UnitRangeTransform, [0.0 -0.5 0.5; 0.0 1.0 2.0])
2×3 Array{Float64,2}:
 0.5  0.0  1.0
 0.0  0.5  1.0
```
"""
function standardize(::Type{DT}, X::AbstractArray{T,2}; kwargs...) where {DT<:DataTransform, T<:Real}
    return transform(fit(DT, X; kwargs...), X)
end
