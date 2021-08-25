### Transformations

abstract type AbstractDataTransform end

# apply the transform
"""
    transform!(t::AbstractDataTransform, x)

Apply transformation `t` to vector or matrix `x` in place.
"""
transform!(t::AbstractDataTransform, x::AbstractMatrix{<:Real}) =
    transform!(x, t, x)
transform!(t::AbstractDataTransform, x::AbstractVector{<:Real}) =
    (transform!(t, reshape(x, :, 1)); x)

"""
    transform(t::AbstractDataTransform, x)

Return a standardized copy of vector or matrix `x` using transformation `t`.
"""
transform(t::AbstractDataTransform, x::AbstractMatrix{<:Real}) =
    transform!(similar(x), t, x)
transform(t::AbstractDataTransform, x::AbstractVector{<:Real}) =
    vec(transform(t, reshape(x, :, 1)))

# reconstruct the original data from transformed values
"""
    reconstruct!(t::AbstractDataTransform, y)

Perform an in-place reconstruction into an original data scale from a transformed
vector or matrix `y` using transformation `t`.
"""
reconstruct!(t::AbstractDataTransform, y::AbstractMatrix{<:Real}) =
    reconstruct!(y, t, y)
reconstruct!(t::AbstractDataTransform, y::AbstractVector{<:Real}) =
    (reconstruct!(t, reshape(y, :, 1)); y)

"""
    reconstruct(t::AbstractDataTransform, y)

Return a reconstruction of an originally scaled data from a transformed vector
or matrix `y` using transformation `t`.
"""
reconstruct(t::AbstractDataTransform, y::AbstractMatrix{<:Real}) =
    reconstruct!(similar(y), t, y)
reconstruct(t::AbstractDataTransform, y::AbstractVector{<:Real}) =
    vec(reconstruct(t, reshape(y, :, 1)))

"""
Standardization (Z-score transformation)
"""
struct ZScoreTransform{T<:Real, U<:AbstractVector{T}} <: AbstractDataTransform
    len::Int
    dims::Int
    mean::U
    scale::U

    function ZScoreTransform(l::Int, dims::Int, m::U, s::U) where {T<:Real, U<:AbstractVector{T}}
        lenm = length(m)
        lens = length(s)
        lenm == l || lenm == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        lens == l || lens == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        new{T, U}(l, dims, m, s)
    end
end

function Base.getproperty(t::ZScoreTransform, p::Symbol)
    if p === :indim || p === :outdim
        return t.len
    else
        return getfield(t, p)
    end
end

"""
    fit(ZScoreTransform, X; dims=nothing, center=true, scale=true)

Fit standardization parameters to vector or matrix `X`
and return a `ZScoreTransform` transformation object.

# Keyword arguments

* `dims`: if `1` fit standardization parameters in column-wise fashion;
  if `2` fit in row-wise fashion. The default is `nothing`, which is equivalent to `dims=2` with a deprecation warning.

* `center`: if `true` (the default) center data so that its mean is zero.

* `scale`: if `true` (the default) scale the data so that its variance is equal to one.

# Examples

```jldoctest
julia> using StatsBase

julia> X = [0.0 -0.5 0.5; 0.0 1.0 2.0]
2×3 Matrix{Float64}:
 0.0  -0.5  0.5
 0.0   1.0  2.0

julia> dt = fit(ZScoreTransform, X, dims=2)
ZScoreTransform{Float64, Vector{Float64}}(2, 2, [0.0, 1.0], [0.5, 1.0])

julia> StatsBase.transform(dt, X)
2×3 Matrix{Float64}:
  0.0  -1.0  1.0
 -1.0   0.0  1.0
```
"""
function fit(::Type{ZScoreTransform}, X::AbstractMatrix{<:Real};
             dims::Union{Integer,Nothing}=nothing, center::Bool=true, scale::Bool=true)
    if dims === nothing
        Base.depwarn("fit(t, x) is deprecated: use fit(t, x, dims=2) instead", :fit)
        dims = 2
    end
    if dims == 1
        n, l = size(X)
        n >= 2 || error("X must contain at least two rows.")
        m, s = mean_and_std(X, 1)
    elseif dims == 2
        l, n = size(X)
        n >= 2 || error("X must contain at least two columns.")
        m, s = mean_and_std(X, 2)
    else
        throw(DomainError(dims, "fit only accept dims to be 1 or 2."))
    end
    return ZScoreTransform(l, dims, (center ? vec(m) : similar(m, 0)),
                                    (scale ? vec(s) : similar(s, 0)))
end

function fit(::Type{ZScoreTransform}, X::AbstractVector{<:Real};
             dims::Integer=1, center::Bool=true, scale::Bool=true)
    if dims != 1
        throw(DomainError(dims, "fit only accepts dims=1 over a vector. Try fit(t, x, dims=1)."))
    end

    return fit(ZScoreTransform, reshape(X, :, 1); dims=dims, center=center, scale=scale)
end

function transform!(y::AbstractMatrix{<:Real}, t::ZScoreTransform, x::AbstractMatrix{<:Real})
    if t.dims == 1
        l = t.len
        size(x,2) == size(y,2) == l || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(y,1)
        size(x,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        m = t.mean
        s = t.scale

        if isempty(m)
            if isempty(s)
                if x !== y
                    copyto!(y, x)
                end
            else
                broadcast!(/, y, x, s')
            end
        else
            if isempty(s)
                broadcast!(-, y, x, m')
            else
                broadcast!((x,m,s)->(x-m)/s, y, x, m', s')
            end
        end
    elseif t.dims == 2
        t_ = ZScoreTransform(t.len, 1, t.mean, t.scale)
        transform!(y', t_, x')
    end
    return y
end

function reconstruct!(x::AbstractMatrix{<:Real}, t::ZScoreTransform, y::AbstractMatrix{<:Real})
    if t.dims == 1
        l = t.len
        size(x,2) == size(y,2) == l || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(y,1)
        size(x,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        m = t.mean
        s = t.scale

        if isempty(m)
            if isempty(s)
                if y !== x
                    copyto!(x, y)
                end
            else
                broadcast!(*, x, y, s')
            end
        else
            if isempty(s)
                broadcast!(+, x, y, m')
            else
                broadcast!((y,m,s)->y*s+m, x, y, m', s')
            end
        end
    elseif t.dims == 2
        t_ = ZScoreTransform(t.len, 1, t.mean, t.scale)
        reconstruct!(x', t_, y')
    end
    return x
end

"""
Unit range normalization
"""
struct UnitRangeTransform{T<:Real, U<:AbstractVector}  <: AbstractDataTransform
    len::Int
    dims::Int
    unit::Bool
    min::U
    scale::U

    function UnitRangeTransform(l::Int, dims::Int, unit::Bool, min::U, max::U) where {T, U<:AbstractVector{T}}
        lenmin = length(min)
        lenmax = length(max)
        lenmin == l || lenmin == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        lenmax == l || lenmax == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        new{T, U}(l, dims, unit, min, max)
    end
end

function Base.getproperty(t::UnitRangeTransform, p::Symbol)
    if p === :indim || p === :outdim
        return t.len
    else
        return getfield(t, p)
    end
end

# fit a unit transform
"""
    fit(UnitRangeTransform, X; dims=nothing, unit=true)

Fit a scaling parameters to vector or matrix `X`
and return a `UnitRangeTransform` transformation object.

# Keyword arguments

* `dims`: if `1` fit standardization parameters in column-wise fashion;
 if `2` fit in row-wise fashion. The default is `nothing`.

* `unit`: if `true` (the default) shift the minimum data to zero.

# Examples

```jldoctest
julia> using StatsBase

julia> X = [0.0 -0.5 0.5; 0.0 1.0 2.0]
2×3 Matrix{Float64}:
 0.0  -0.5  0.5
 0.0   1.0  2.0

julia> dt = fit(UnitRangeTransform, X, dims=2)
UnitRangeTransform{Float64, Vector{Float64}}(2, 2, true, [-0.5, 0.0], [1.0, 0.5])

julia> StatsBase.transform(dt, X)
2×3 Matrix{Float64}:
 0.5  0.0  1.0
 0.0  0.5  1.0
```
"""
function fit(::Type{UnitRangeTransform}, X::AbstractMatrix{<:Real};
             dims::Union{Integer,Nothing}=nothing, unit::Bool=true)
    if dims === nothing
        Base.depwarn("fit(t, x) is deprecated: use fit(t, x, dims=2) instead", :fit)
        dims = 2
    end
    dims ∈ (1, 2) || throw(DomainError(dims, "fit only accept dims to be 1 or 2."))
    tmin, tmax = _compute_extrema(X, dims)
    @. tmax = 1 / (tmax - tmin)
    l = length(tmin)
    return UnitRangeTransform(l, dims, unit, tmin, tmax)
end

function _compute_extrema(X::AbstractMatrix, dims::Integer)
    dims == 2 && return _compute_extrema(X', 1)
    l = size(X, 2)
    tmin = similar(X, l)
    tmax = similar(X, l)
    for i in 1:l
        @inbounds tmin[i], tmax[i] = extrema(@view(X[:, i]))
    end
    return tmin, tmax
end

function fit(::Type{UnitRangeTransform}, X::AbstractVector{<:Real};
             dims::Integer=1, unit::Bool=true)
    if dims != 1
        throw(DomainError(dims, "fit only accept dims=1 over a vector. Try fit(t, x, dims=1)."))
    end
    tmin, tmax = extrema(X)
    tmax = 1 / (tmax - tmin)
    return UnitRangeTransform(1, dims, unit, [tmin], [tmax])
end

function transform!(y::AbstractMatrix{<:Real}, t::UnitRangeTransform, x::AbstractMatrix{<:Real})
    if t.dims == 1
        l = t.len
        size(x,2) == size(y,2) == l || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(x,1)
        size(y,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        tmin = t.min
        tscale = t.scale

        if t.unit
            broadcast!((x,s,m)->(x-m)*s, y, x, tscale', tmin')
        else
            broadcast!(*, y, x, tscale')
        end
    elseif t.dims == 2
        t_ = UnitRangeTransform(t.len, 1, t.unit, t.min, t.scale)
        transform!(y', t_, x')
    end
    return y
end

function reconstruct!(x::AbstractMatrix{<:Real}, t::UnitRangeTransform, y::AbstractMatrix{<:Real})
    if t.dims == 1
        l = t.len
        size(x,2) == size(y,2) == l || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(y,1)
        size(x,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        tmin = t.min
        tscale = t.scale

        if t.unit
            broadcast!((y,s,m)->y/s+m, x, y, tscale', tmin')
        else
            broadcast!(/, x, y, tscale')
        end
    elseif t.dims == 2
        t_ = UnitRangeTransform(t.len, 1, t.unit, t.min, t.scale)
        reconstruct!(x', t_, y')
    end
    return x
end

"""
    standardize(DT, X; dims=nothing, kwargs...)

 Return a standardized copy of vector or matrix `X` along dimensions `dims`
 using transformation `DT` which is a subtype of `AbstractDataTransform`:

- `ZScoreTransform`
- `UnitRangeTransform`

# Example

```jldoctest
julia> using StatsBase

julia> standardize(ZScoreTransform, [0.0 -0.5 0.5; 0.0 1.0 2.0], dims=2)
2×3 Matrix{Float64}:
  0.0  -1.0  1.0
 -1.0   0.0  1.0

julia> standardize(UnitRangeTransform, [0.0 -0.5 0.5; 0.0 1.0 2.0], dims=2)
2×3 Matrix{Float64}:
 0.5  0.0  1.0
 0.0  0.5  1.0
```
"""
function standardize(::Type{DT}, X::AbstractVecOrMat{<:Real}; kwargs...) where {DT <: AbstractDataTransform}
    return transform(fit(DT, X; kwargs...), X)
end
