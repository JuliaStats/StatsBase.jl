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

indim(t::ZScoreTransform) = t.dim
outdim(t::ZScoreTransform) = t.dim

# fit a z-score transform
function fit(::Type{ZScoreTransform}, X::AbstractArray{T,2}; center::Bool=true, scale::Bool=true) where T<:Real
    d, n = size(X)
    n >= 2 || error("X must contain at least two columns.")

    m, s = mean_and_std(X, 2)

    return ZScoreTransform(d, (center ? vec(m) : zeros(T, 0)),
                              (scale ? vec(s) : zeros(T, 0)))
end

function transform!(y::AbstractArray{YT,1}, t::ZScoreTransform, x::AbstractArray{XT,1}) where {YT<:Real,XT<:Real}
    d = t.dim
    length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if x !== y
                copyto!(y, x)
            end
        else
            for i = 1:d
                @inbounds y[i] = x[i] / s[i]
            end
        end
    else
        if isempty(s)
            for i = 1:d
                @inbounds y[i] = x[i] - m[i]
            end
        else
            for i = 1:d
                @inbounds y[i] = (x[i] - m[i]) / s[i]
            end
        end
    end
    return y
end

function transform!(y::AbstractArray{YT,2}, t::ZScoreTransform, x::AbstractArray{XT,2}) where {YT<:Real,XT<:Real}
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(x,2)
    size(y,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if x !== y
                copyto!(y, x)
            end
        else
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds yj[i] = xj[i] / s[i]
                end
            end
        end
    else
        if isempty(s)
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds yj[i] = xj[i] - m[i]
                end
            end
        else
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds yj[i] = (xj[i] - m[i]) / s[i]
                end
            end
        end
    end
    return y
end

function reconstruct!(x::AbstractArray{XT,1}, t::ZScoreTransform, y::AbstractArray{YT,1}) where {YT<:Real,XT<:Real}
    d = t.dim
    length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if y !== x
                copyto!(x, y)
            end
        else
            for i = 1:d
                @inbounds x[i] = y[i] * s[i]
            end
        end
    else
        if isempty(s)
            for i = 1:d
                @inbounds x[i] = y[i] + m[i]
            end
        else
            for i = 1:d
                @inbounds x[i] = y[i] * s[i] + m[i]
            end
        end
    end
    return x
end

function reconstruct!(x::AbstractArray{XT,2}, t::ZScoreTransform, y::AbstractArray{YT,2}) where {YT<:Real,XT<:Real}
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
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds xj[i] = yj[i] * s[i]
                end
            end
        end
    else
        if isempty(s)
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds xj[i] = yj[i] + m[i]
                end
            end
        else
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds xj[i] = yj[i] * s[i] + m[i]
                end
            end
        end
    end
    return x
end

# UnitRangeTransform normalization

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

indim(t::UnitRangeTransform) = t.dim
outdim(t::UnitRangeTransform) = t.dim

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

function transform!(y::AbstractArray{YT,1}, t::UnitRangeTransform, x::AbstractArray{XT,1}) where {YT<:Real,XT<:Real}
    d = t.dim
    length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        for i = 1:d
            @inbounds y[i] = (x[i] - tmin[i]) * tscale[i]
        end
    else
        for i = 1:d
            @inbounds y[i] = x[i] * tscale[i]
        end
    end
    return y
end

function transform!(y::AbstractArray{YT,2}, t::UnitRangeTransform, x::AbstractArray{XT,2}) where {YT<:Real,XT<:Real}
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(x,2)
    size(y,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        for j = 1:n
            xj = view(x, :, j)
            yj = view(y, :, j)
            for i = 1:d
                @inbounds yj[i] = (xj[i] - tmin[i]) * tscale[i]
            end
        end
    else
        for j = 1:n
            xj = view(x, :, j)
            yj = view(y, :, j)
            for i = 1:d
                @inbounds yj[i] = xj[i] * tscale[i]
            end
        end
    end
    return y
end

function reconstruct!(x::AbstractArray{XT,1}, t::UnitRangeTransform, y::AbstractArray{YT,1}) where {YT<:Real,XT<:Real}
    d = t.dim
    length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        for i = 1:d
            @inbounds x[i] = y[i] / tscale[i] +  tmin[i]
        end
    else
        for i = 1:d
            @inbounds x[i] = y[i] / tscale[i]
        end
    end
    return x
end

function reconstruct!(x::AbstractArray{XT,2}, t::UnitRangeTransform, y::AbstractArray{YT,2}) where {YT<:Real,XT<:Real}
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(y,2)
    size(x,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    tmin = t.min
    tscale = t.scale

    if t.unit
        for j = 1:n
            xj = view(x, :, j)
            yj = view(y, :, j)
            for i = 1:d
                @inbounds xj[i] = yj[i] / tscale[i] + tmin[i]
            end
        end
    else
        for j = 1:n
            xj = view(x, :, j)
            yj = view(y, :, j)
            for i = 1:d
                @inbounds xj[i] = yj[i] / tscale[i]
            end
        end
    end
    return x
end
