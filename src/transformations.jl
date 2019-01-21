### Transformations

abstract type AbstractDataTransform end

# apply the transform
"""
    transform!(t::AbstractDataTransform, x)

Apply transformation `t` to vector or matrix `x` in place.
"""
transform!(t::AbstractDataTransform, x::AbstractArray{<:Real,1}) = transform!(x, t, x)
transform!(t::AbstractDataTransform, x::AbstractArray{<:Real,2}) = transform!(x, t, x)

"""
    transform(t::AbstractDataTransform, x)

Return a row-standardized vector or matrix `x` using `t` transformation.
"""
transform(t::AbstractDataTransform, x::AbstractArray{<:Real,1}) = transform!(similar(x), t, x)
transform(t::AbstractDataTransform, x::AbstractArray{<:Real,2}) = transform!(similar(x), t, x)

# reconstruct the original data from transformed values
"""
    reconstruct!(t::AbstractDataTransform, y)

Perform an in-place reconstruction into an original data scale from a row-transformed
vector or matrix `y` using `t` transformation.
"""
reconstruct!(t::AbstractDataTransform, y::AbstractArray{<:Real,1}) = reconstruct!(y, t, y)
reconstruct!(t::AbstractDataTransform, y::AbstractArray{<:Real,2}) = reconstruct!(y, t, y)

"""
    reconstruct(t::AbstractDataTransform, y)

Return a reconstruction of an originally scaled data from a row-transformed vector
or matrix `y` using `t` transformation.
"""
reconstruct(t::AbstractDataTransform, y::AbstractArray{<:Real,1})  = reconstruct!(similar(y), t, y)
reconstruct(t::AbstractDataTransform, y::AbstractArray{<:Real,2})  = reconstruct!(similar(y), t, y)

"""
    Standardization (Z-score transformation)
"""
struct ZScoreTransform{T<:Real} <: AbstractDataTransform
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
    if p === :indim || p === :outdim
        return t.dim
    else
        return getfield(t, p)
    end
end

"""
    fit(ZScoreTransform, X; center=true, scale=true)

Fit standardization parameters to `X` and return a `ZScoreTransform` transformation object.

# Arguments

* `data`: matrix  of samples to fit transformation parameters.

# Keyword arguments

* `center`: if `true` (the default) center data so that its mean is zero.

* `scale`: if `true` (the default) scale the data so that its variance is equal to one.

# Examples

```jldoctest
julia> using StatsBase

julia> X = [0.0 -0.5 0.5; 0.0 1.0 2.0]
2×3 Array{Float64,2}:
 0.0  -0.5  0.5
 0.0   1.0  2.0

julia> dt = fit(ZScoreTransform, X)
ZScoreTransform{Float64}(2, [0.0, 1.0], [0.5, 1.0])

julia> StatsBase.transform(dt, X)
2×3 Array{Float64,2}:
  0.0  -1.0  1.0
 -1.0   0.0  1.0
```
"""
function fit(::Type{ZScoreTransform}, X::AbstractArray{<:Real,2}; center::Bool=true, scale::Bool=true)
    d, n = size(X)
    n >= 2 || error("X must contain at least two columns.")

    T = eltype(X)
    m, s = mean_and_std(X, 2)

    return ZScoreTransform(d, (center ? vec(m) : zeros(T, 0)),
                              (scale ? vec(s) : zeros(T, 0)))
end

function transform!(y::AbstractVecOrMat{<:Real}, t::ZScoreTransform, x::AbstractVecOrMat{<:Real})
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

function reconstruct!(x::AbstractVecOrMat{<:Real}, t::ZScoreTransform, y::AbstractVecOrMat{<:Real})
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

"""
    Unit range normalization
"""
struct UnitRangeTransform{T<:Real}  <: AbstractDataTransform
    dim::Int
    unit::Bool
    min::Vector{T}
    scale::Vector{T}

    function UnitRangeTransform(d::Int, unit::Bool, min::Vector{T}, max::Vector{T}) where {T}
        lenmin = length(min)
        lenmax = length(max)
        lenmin == d || lenmin == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        lenmax == d || lenmax == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        new{T}(d, unit, min, max)
    end
end

function Base.getproperty(t::UnitRangeTransform, p::Symbol)
    if p === :indim || p === :outdim
        return t.dim
    else
        return getfield(t, p)
    end
end

# fit a unit transform
"""
    fit(UnitRangeTransform, X; center=true, scale=true)

Fit a scaling parameters to `X` and return transformation description.

# Arguments

* `data`: matrix  of samples to fit transformation parameters.

# Keyword arguments

* `center`: if `true` (the default) centere data around zero.

* `scale`: if `true` (the default) perform variance scaling.

# Examples

```jldoctest
julia> using StatsBase

julia> X = [0.0 -0.5 0.5; 0.0 1.0 2.0]
2×3 Array{Float64,2}:
 0.0  -0.5  0.5
 0.0   1.0  2.0

julia> dt = fit(UnitRangeTransform, X)
UnitRangeTransform{Float64}(2, true, [-0.5, 0.0], [1.0, 0.5])

julia> StatsBase.transform(dt, X)
2×3 Array{Float64,2}:
 0.5  0.0  1.0
 0.0  0.5  1.0
```
"""
function fit(::Type{UnitRangeTransform}, X::AbstractArray{<:Real,2}; unit::Bool=true)
    d, n = size(X)

    tmin = X[:, 1]
    tmax = X[:, 1]
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
        @inbounds tmax[i] = 1 / (tmax[i] - tmin[i])
    end
    return UnitRangeTransform(d, unit, tmin, tmax)
end

function transform!(y::AbstractVecOrMat{<:Real}, t::UnitRangeTransform, x::AbstractVecOrMat{<:Real})
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

function reconstruct!(x::AbstractVecOrMat{<:Real}, t::UnitRangeTransform, y::AbstractVecOrMat{<:Real})
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

Return a row-standardized matrix `X` using `DT` transformation which is a subtype of `AbstractDataTransform`:

- `ZScoreTransform`
- `UnitRangeTransform`

# Example

```jldoctest
julia> using StatsBase

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
function standardize(::Type{DT}, X::AbstractArray{<:Real,2}; kwargs...) where {DT<:AbstractDataTransform}
    return transform(fit(DT, X; kwargs...), X)
end
