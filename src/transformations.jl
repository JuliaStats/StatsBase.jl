### Transformation

abstract DataTransform

# apply the transform
transform!(t::S, x::AbstractArray{T,1}) where {T<:Real, S<:DataTransform} = transform!(x, t, x)
transform!(t::S, x::AbstractArray{T,2}) where {T<:Real, S<:DataTransform} = transform!(x, t, x)

transform{T<:Real, S<:DataTransform}(t::S, x::DenseArray{T,1}) = transform!(Array(Float64, size(x)), t, x)
transform{T<:Real, S<:DataTransform}(t::S, x::DenseArray{T,2}) = transform!(Array(Float64, size(x)), t, x)

# reconstruct the original data from transformed values
reconstruct!{T<:FloatingPoint, S<:DataTransform}(t::S, x::DenseArray{T,1}) = reconstruct!(x, t, x)
reconstruct!{T<:FloatingPoint, S<:DataTransform}(t::S, x::DenseArray{T,2}) = reconstruct!(x, t, x)

reconstruct{T<:Real, S<:DataTransform}(t::S, y::DenseArray{T,1}) = reconstruct!(Array(Float64, size(y)), t, y)
reconstruct{T<:Real, S<:DataTransform}(t::S, y::DenseArray{T,2}) = reconstruct!(Array(Float64, size(y)), t, y)

# Z-score transformation
immutable ZScoreTransform <: DataTransform
    dim::Int
    mean::Vector{Float64}
    scale::Vector{Float64}

    function ZScoreTransform(d::Int, m::Vector{Float64}, s::Vector{Float64})
        lenm = length(m)
        lens = length(s)
        lenm == d || lenm == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        lens == d || lens == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        new(d, m, s)
    end
end

indim(t::ZScoreTransform) = t.dim
outdim(t::ZScoreTransform) = t.dim

# fit a z-score transform
function fit{T<:Real}(::Type{ZScoreTransform}, X::DenseArray{T,2}; center::Bool=true, scale::Bool=true)
    d, n = size(X)
    n >= 2 || error("X must contain at least two columns.")

    m, s = mean_and_std(X, 2)

    return ZScoreTransform(d, (center ? float64(vec(m)) : Array(Float64, 0)),
                              (scale ? float64(vec(s)) : Array(Float64, 0)))
end


function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,1}, t::ZScoreTransform, x::DenseArray{XT,1})
    d = t.dim
    length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if !is(x, y)
                copy!(y, x)
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

function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,2}, t::ZScoreTransform, x::DenseArray{XT,2})
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(x,2)
    size(y,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if !is(x, y)
                copy!(y, x)
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

function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,1}, t::ZScoreTransform, y::DenseArray{YT,1})
    d = t.dim
    length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if !is(y, x)
                copy!(x, y)
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

function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,2}, t::ZScoreTransform, y::DenseArray{YT,2})
    d = t.dim
    size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
    n = size(y,2)
    size(x,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

    m = t.mean
    s = t.scale

    if isempty(m)
        if isempty(s)
            if !is(y, x)
                copy!(x, y)
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

immutable UnitRangeTransform  <: DataTransform
    dim::Int
    unit::Bool
    min::Vector{Float64}
    scale::Vector{Float64}

    function UnitRangeTransform(d::Int, unit::Bool, min::Vector{Float64}, max::Vector{Float64})
        lenmin = length(min)
        lenmax = length(max)
        lenmin == d || lenmin == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        lenmax == d || lenmax == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
        new(d, unit, min, max)
    end
end

indim(t::UnitRangeTransform) = t.dim
outdim(t::UnitRangeTransform) = t.dim

function fit{T<:Real}(::Type{UnitRangeTransform}, X::DenseArray{T,2}; unit::Bool=true)
    d, n = size(X)

    tmin = Array(Float64, d)
    tmax = Array(Float64, d)
    copy!(tmin, X[:, 1])
    copy!(tmax, X[:, 1])
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

function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,1}, t::UnitRangeTransform, x::DenseArray{XT,1})
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

function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,2}, t::UnitRangeTransform, x::DenseArray{XT,2})
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

function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,1}, t::UnitRangeTransform, y::DenseArray{YT,1})
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

function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,2}, t::UnitRangeTransform, y::DenseArray{YT,2})
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