### Transformations

abstract type AbstractDataTransform end

# apply the transform
"""
    transform!(t::AbstractDataTransform, x; dims=nothing)

Apply transformation `t` to vector or matrix `x` in place along dimension `dims`.
"""
function transform!(t::AbstractDataTransform, x::AbstractVecOrMat{<:Real}; dims::Union{Integer,Nothing}=nothing)
    if dims == 1
        return transform!(x, t, x)
    elseif dims == 2
        return transform!(x', t, x')
    elseif dims === nothing
        Base.depwarn("transform!(t, x) is deprecated: use transform!(t, x, dims=2) instead", :transform)
    end
end

"""
    transform(t::AbstractDataTransform, x; dims=nothing)

Return a standardized vector or matrix `x` using `t` transformation along dimension `dims`.
"""
function transform(t::AbstractDataTransform, x::AbstractVecOrMat{<:Real}; dims::Union{Integer,Nothing}=nothing)
    if dims == 1
        return transform!(similar(x), t, x)
    elseif dims == 2
        return transform!(similar(x)', t, x')
    elseif dims === nothing
        Base.depwarn("transform!(t, x) is deprecated: use transform!(t, x, dims=2) instead", :transform)
    end
end

# reconstruct the original data from transformed values
"""
    reconstruct!(t::AbstractDataTransform, y)

Perform an in-place reconstruction into an original data scale from a row-transformed
vector or matrix `y` using `t` transformation.
"""
reconstruct!(t::AbstractDataTransform, y::AbstractVecOrMat{<:Real}) = reconstruct!(y, t, y)

"""
    reconstruct(t::AbstractDataTransform, y)

Return a reconstruction of an originally scaled data from a row-transformed vector
or matrix `y` using `t` transformation.
"""
reconstruct(t::AbstractDataTransform, y::AbstractVecOrMat{<:Real}) = reconstruct!(similar(y), t, y)

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
    fit(ZScoreTransform, X; dims=nothing, center=true, scale=true)

Fit standardization parameters to `X` and return a `ZScoreTransform` transformation object.

# Arguments

* `data`: matrix  of samples to fit transformation parameters.

# Keyword arguments

* `dims`: if `1` (the default) fit standardization parameters in column-wise fashion;
  if `2` fit in row-wise fashion.

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
function fit(::Type{ZScoreTransform}, X::AbstractMatrix{<:Real};
             dims::Union{Integer,Nothing}=nothing, center::Bool=true, scale::Bool=true)
    if dims == 1
        n, d = size(X)
        n >= 2 || error("X must contain at least two rows.")
        m, s = mean_and_std(X, 1)
    elseif dims == 2
        d, n = size(X)
        n >= 2 || error("X must contain at least two columns.")
        m, s = mean_and_std(X, 2)
    elseif dims === nothing
        Base.depwarn("fit(t, x) is deprecated: use fit(t, x, dims=2) instead", :fit)
    end
    T = eltype(X)
    return ZScoreTransform(d, (center ? vec(m) : zeros(T, 0)),
                              (scale ? vec(s) : zeros(T, 0)))
end

function fit(::Type{ZScoreTransform}, X::AbstractVector{<:Real}; center::Bool=true, scale::Bool=true)
    T = eltype(X)
    m, s = mean_and_std(X)
    return ZScoreTransform(1, (center ? [m] : zeros(T, 0)),
                              (scale ? [s] : zeros(T, 0)))
end

function transform!(y::AbstractMatrix{<:Real}, t::ZScoreTransform, x::AbstractMatrix{<:Real})
    d = t.dim
    size(x,2) == size(y,2) == d || throw(DimensionMismatch("Inconsistent dimensions."))
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
    return y
end

function transform!(y::AbstractVector{<:Real}, t::ZScoreTransform, x::AbstractVector{<:Real})
    t.dim == 1 || throw(DimensionMismatch("Inconsistent dimensions."))
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
            broadcast!(/, y, x, s[1])
        end
    else
        if isempty(s)
            broadcast!(-, y, x, m[1])
        else
            broadcast!((x,m,s)->(x-m)/s, y, x, m[1], s[1])
        end
    end
    return y
end

function reconstruct!(x::AbstractMatrix{<:Real}, t::ZScoreTransform, y::AbstractMatrix{<:Real})
    d = t.dim
    size(x,2) == size(y,2) == d || throw(DimensionMismatch("Inconsistent dimensions."))
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
    return x
end

function reconstruct!(x::AbstractVector{<:Real}, t::ZScoreTransform, y::AbstractVector{<:Real})
    t.dim == 1 || throw(DimensionMismatch("Inconsistent dimensions."))
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
            broadcast!(*, x, y, s[1])
        end
    else
        if isempty(s)
            broadcast!(+, x, y, m[1])
        else
            broadcast!((y,m,s)->y*s+m, x, y, m[1], s[1])
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
    fit(UnitRangeTransform, X; dims=nothing, unit=true)

Fit a scaling parameters to `X` and return transformation description.

# Arguments

* `data`: matrix  of samples to fit transformation parameters.

# Keyword arguments

* `dims`: if `1` (the default) fit standardization parameters in column-wise fashion; if `2` fit in row-wise fashion.

* `unit`: if `true` (the default) shift the minimum data to zero.

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
function fit(::Type{UnitRangeTransform}, X::AbstractMatrix{<:Real};
             dims::Union{Integer,Nothing}=nothing, unit::Bool=true)
    if dims == 1
        d, tmin, tmax = _extract_info(X)
    elseif dims == 2
        d, tmin, tmax = _extract_info(X')
    elseif dims == nothing
        Base.depwarn("fit(t, x) is deprecated: use fit(t, x, dims=2) instead", :fit)
    end

    for i = 1:d
        @inbounds tmax[i] = 1 / (tmax[i] - tmin[i])
    end
    return UnitRangeTransform(d, unit, tmin, tmax)
end

function _extract_info(X::AbstractMatrix{<:Real})
    n, d = size(X)
    tmin = X[1, :]
    tmax = X[1, :]
    for j = 1:d
        @inbounds for i = 2:n
            if X[i, j] < tmin[j]
                tmin[j] = X[i, j]
            elseif X_[i, j] > tmax[j]
                tmax[j] = X[i, j]
            end
        end
    end
    return d, tmin, tmax
end

function fit(::Type{UnitRangeTransform}, X::AbstractVector{<:Real}; unit::Bool=true)
    n = length(X)
    tmin = X[1]
    tmax = X[1]
    @inbounds for i = 2:n
        if X[i] < tmin
            tmin = X[i]
        elseif X[i] > tmax
            tmax = X[i]
        end
    end
    tmax = 1 / (tmax - tmin)
    return UnitRangeTransform(1, unit, [tmin], [tmax])
end

function transform!(y::AbstractMatrix{<:Real}, t::UnitRangeTransform, x::AbstractMatrix{<:Real})
    d = t.dim
    size(x,2) == size(y,2) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(x,1)
    size(y,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        broadcast!((x,s,m)->(x-m)*s, y, x, tscale', tmin')
    else
        broadcast!(*, y, x, tscale')
    end
    return y
end

function transform!(y::AbstractVector{<:Real}, t::UnitRangeTransform, x::AbstractVector{<:Real})
    t.dim == 1 || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(x,1)
    size(y,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        broadcast!((x,s,m)->(x-m)*s, y, x, tscale[1], tmin[1])
    else
        broadcast!(*, y, x, tscale[1])
    end
    return y
end

function reconstruct!(x::AbstractMatrix{<:Real}, t::UnitRangeTransform, y::AbstractMatrix{<:Real})
    d = t.dim
    size(x,2) == size(y,2) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(y,1)
    size(x,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        broadcast!((y,s,m)->y/s+m, x, y, tscale', tmin')
    else
        broadcast!(/, x, y, tscale')
    end
    return x
end

function reconstruct!(x::AbstractVector{<:Real}, t::UnitRangeTransform, y::AbstractVector{<:Real})
    t.dim == 1 || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(y,1)
    size(x,1) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        broadcast!((y,s,m)->y/s+m, x, y, tscale[1], tmin[1])
    else
        broadcast!(/, x, y, tscale[1])
    end
    return x
end

"""
    standardize(DT, X; kwargs...)

 Return a column-standardized copy of vector or matrix `X` using transformation `DT`
 which is a subtype of `AbstractDataTransform`:

- `ZScoreTransform`
- `UnitRangeTransform`

Return a column-standardized matrix while the input of `X` is a `AbstractMatrix`.
Return a standardized vector while the input of `X` is a `AbstractVector`.

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
function standardize(::Type{DT}, X::AbstractVector{<:Real}; kwargs...) where {DT <: AbstractDataTransform}
    return transform(fit(DT, X; kwargs...), X)
end

function standardize(::Type{DT}, X::AbstractMatrix{<:Real}; kwargs...) where {DT <: AbstractDataTransform}
    return transform(fit(DT, X; kwargs...), X)
end
