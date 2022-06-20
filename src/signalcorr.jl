# Correlation analysis of signals
#
#  autocorrelation
#  cross-correlation
#  partial autocorrelation
#

#######################################
#
#   Helper functions
#
#######################################

default_laglen(lx::Int) = min(lx-1, round(Int,10*log10(lx)))
check_lags(lx::Int, lags::AbstractVector) = (maximum(lags) < lx || error("lags must be less than the sample length."))

function demean_col!(z::RealVector, x::RealMatrix, j::Int, demean::Bool)
    T = eltype(z)
    m = size(x, 1)
    @assert m == length(z)
    b = m * (j-1)
    if demean
        s = zero(T)
        for i = 1 : m
            s += x[b + i]
        end
        mv = s / m
        for i = 1 : m
            z[i] = x[b + i] - mv
        end
    else
        copyto!(z, 1, x, b+1, m)
    end
    z
end


#######################################
#
#   Auto-correlations
#
#######################################

default_autolags(lx::Int) = 0 : default_laglen(lx)

_autodot(x::AbstractVector{<:RealFP}, lx::Int, l::Int) = dot(x, 1:(lx-l), x, (1+l):lx)
_autodot(x::RealVector, lx::Int, l::Int) = dot(view(x, 1:(lx-l)), view(x, (1+l):lx))


## autocov
"""
    autocov!(r, x, lags; demean=true)

Compute the autocovariance of a vector or matrix `x` at `lags` and store the result
in `r`. `demean` denotes whether the mean of `x` should be subtracted from `x`
before computing the autocovariance.

If `x` is a vector, `r` must be a vector of the same length as `lags`.
If `x` is a matrix, `r` must be a matrix of size `(length(lags), size(x,2))`, and
where each column in the result will correspond to a column in `x`.

The output is not normalized. See [`autocor!`](@ref) for a method with normalization.
"""
function autocov!(r::RealVector, x::RealVector, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    length(r) == m || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    z::Vector{T} = demean ? x .- mean(x) : x
    for k = 1 : m  # foreach lag value
        r[k] = _autodot(z, lx, lags[k]) / lx
    end
    return r
end

function autocov!(r::RealMatrix, x::RealMatrix, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    size(r) == (m, ns) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    z = Vector{T}(undef, lx)
    for j = 1 : ns
        demean_col!(z, x, j, demean)
        for k = 1 : m
            r[k,j] = _autodot(z, lx, lags[k]) / lx
        end
    end
    return r
end


"""
    autocov(x, [lags]; demean=true)

Compute the autocovariance of a vector or matrix `x`, optionally specifying
the `lags` at which to compute the autocovariance. `demean` denotes whether
the mean of `x` should be subtracted from `x` before computing the autocovariance.

If `x` is a vector, return a vector of the same length as `lags`.
If `x` is a matrix, return a matrix of size `(length(lags), size(x,2))`,
where each column in the result corresponds to a column in `x`.

When left unspecified, the lags used are the integers from 0 to
`min(size(x,1)-1, 10*log10(size(x,1)))`.

The output is not normalized. See [`autocor`](@ref) for a function with normalization.
"""
function autocov(x::RealVector, lags::IntegerVector; demean::Bool=true)
    out = Vector{float(eltype(x))}(undef, length(lags))
    autocov!(out, x, lags; demean=demean)
end

function autocov(x::RealMatrix, lags::IntegerVector; demean::Bool=true)
    out = Matrix{float(eltype(x))}(undef, length(lags), size(x,2))
    autocov!(out, x, lags; demean=demean)
end

autocov(x::AbstractVecOrMat{<:Real}; demean::Bool=true) =
    autocov(x, default_autolags(size(x,1)); demean=demean)

## autocor

"""
    autocor!(r, x, lags; demean=true)

Compute the autocorrelation function (ACF) of a vector or matrix `x` at `lags`
and store the result in `r`. `demean` denotes whether the mean of `x` should
be subtracted from `x` before computing the ACF.

If `x` is a vector, `r` must be a vector of the same length as `lags`.
If `x` is a matrix, `r` must be a matrix of size `(length(lags), size(x,2))`, and
where each column in the result will correspond to a column in `x`.

The output is normalized by the variance of `x`, i.e. so that the lag 0
autocorrelation is 1. See [`autocov!`](@ref) for the unnormalized form.
"""
function autocor!(r::RealVector, x::RealVector, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    length(r) == m || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    z::Vector{T} = demean ? x .- mean(x) : x
    zz = dot(z, z)
    for k = 1 : m  # foreach lag value
        r[k] = _autodot(z, lx, lags[k]) / zz
    end
    return r
end

function autocor!(r::RealMatrix, x::RealMatrix, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    size(r) == (m, ns) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    z = Vector{T}(undef, lx)
    for j = 1 : ns
        demean_col!(z, x, j, demean)
        zz = dot(z, z)
        for k = 1 : m
            r[k,j] = _autodot(z, lx, lags[k]) / zz
        end
    end
    return r
end


"""
    autocor(x, [lags]; demean=true)

Compute the autocorrelation function (ACF) of a vector or matrix `x`,
optionally specifying the `lags`. `demean` denotes whether the mean
of `x` should be subtracted from `x` before computing the ACF.

If `x` is a vector, return a vector of the same length as `lags`.
If `x` is a matrix, return a matrix of size `(length(lags), size(x,2))`,
where each column in the result corresponds to a column in `x`.

When left unspecified, the lags used are the integers from 0 to
`min(size(x,1)-1, 10*log10(size(x,1)))`.

The output is normalized by the variance of `x`, i.e. so that the lag 0
autocorrelation is 1. See [`autocov`](@ref) for the unnormalized form.
"""
function autocor(x::RealVector, lags::IntegerVector; demean::Bool=true)
    out = Vector{float(eltype(x))}(undef, length(lags))
    autocor!(out, x, lags; demean=demean)
end

function autocor(x::RealMatrix, lags::IntegerVector; demean::Bool=true)
    out = Matrix{float(eltype(x))}(undef, length(lags), size(x,2))
    autocor!(out, x, lags; demean=demean)
end

autocor(x::AbstractVecOrMat{<:Real}; demean::Bool=true) =
    autocor(x, default_autolags(size(x,1)); demean=demean)


#######################################
#
#   Cross-correlations
#
#######################################

default_crosslags(lx::Int) = (l=default_laglen(lx); -l:l)

function _crossdot(x::AbstractVector{T}, y::AbstractVector{T}, lx::Int, l::Int) where {T<:RealFP}
    if l >= 0
        dot(x, 1:(lx-l), y, (1+l):lx)
    else
        dot(x, (1-l):lx, y, 1:(lx+l))
    end
end
function _crossdot(x::RealVector, y::RealVector, lx::Int, l::Int)
    if l >= 0
        dot(view(x, 1:(lx-l)), view(y, (1+l):lx))
    else
        dot(view(x, (1-l):lx), view(y, 1:(lx+l)))
    end
end

## crosscov

"""
    crosscov!(r, x, y, lags; demean=true)

Compute the cross covariance function (CCF) between real-valued vectors or matrices
`x` and `y` at `lags` and store the result in `r`. `demean` specifies whether the
respective means of `x` and `y` should be subtracted from them before computing their
CCF.

If both `x` and `y` are vectors, `r` must be a vector of the same length as
`lags`. If either `x` is a matrix and `y` is a vector, `r` must be a matrix of size
`(length(lags), size(x, 2))`; if `x` is a vector and `y` is a matrix, `r` must be a matrix
of size `(length(lags), size(y, 2))`. If both `x` and `y` are matrices, `r` must be a
three-dimensional array of size `(length(lags), size(x, 2), size(y, 2))`.

The output is not normalized. See [`crosscor!`](@ref) for a function with normalization.
"""
function crosscov!(r::RealVector, x::RealVector, y::RealVector, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    (length(y) == lx && length(r) == m) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    zx::Vector{T} = demean ? x .- mean(x) : x
    S = typeof(zero(eltype(y)) / 1)
    zy::Vector{S} = demean ? y .- mean(y) : y
    for k = 1 : m  # foreach lag value
        r[k] = _crossdot(zx, zy, lx, lags[k]) / lx
    end
    return r
end

function crosscov!(r::RealMatrix, x::RealMatrix, y::RealVector, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    (length(y) == lx && size(r) == (m, ns)) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    zx = Vector{T}(undef, lx)
    S = typeof(zero(eltype(y)) / 1)
    zy::Vector{S} = demean ? y .- mean(y) : y
    for j = 1 : ns
        demean_col!(zx, x, j, demean)
        for k = 1 : m
            r[k,j] = _crossdot(zx, zy, lx, lags[k]) / lx
        end
    end
    return r
end

function crosscov!(r::RealMatrix, x::RealVector, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    ns = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, ns)) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    zx::Vector{T} = demean ? x .- mean(x) : x
    S = typeof(zero(eltype(y)) / 1)
    zy = Vector{S}(undef, lx)
    for j = 1 : ns
        demean_col!(zy, y, j, demean)
        for k = 1 : m
            r[k,j] = _crossdot(zx, zy, lx, lags[k]) / lx
        end
    end
    return r
end

function crosscov!(r::AbstractArray{<:Real,3}, x::RealMatrix, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    nx = size(x, 2)
    ny = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, nx, ny)) || throw(DimensionMismatch())
    check_lags(lx, lags)

    # cached (centered) columns of x
    T = typeof(zero(eltype(x)) / 1)
    zxs = Vector{T}[]
    sizehint!(zxs, nx)
    for j = 1 : nx
        xj = x[:,j]
        if demean
            mv = mean(xj)
            for i = 1 : lx
                xj[i] -= mv
            end
        end
        push!(zxs, xj)
    end

    S = typeof(zero(eltype(y)) / 1)
    zy = Vector{S}(undef, lx)
    for j = 1 : ny
        demean_col!(zy, y, j, demean)
        for i = 1 : nx
            zx = zxs[i]
            for k = 1 : m
                r[k,i,j] = _crossdot(zx, zy, lx, lags[k]) / lx
            end
        end
    end
    return r
end


"""
    crosscov(x, y, [lags]; demean=true)

Compute the cross covariance function (CCF) between real-valued vectors or
matrices `x` and `y`, optionally specifying the `lags`. `demean` specifies
whether the respective means of `x` and `y` should be subtracted from them
before computing their CCF.

If both `x` and `y` are vectors, return a vector of the same length as
`lags`. Otherwise, compute cross covariances between each pairs of columns in `x` and `y`.

When left unspecified, the lags used are the integers from
`-min(size(x,1)-1, 10*log10(size(x,1)))` to `min(size(x,1), 10*log10(size(x,1)))`.

The output is not normalized. See [`crosscor`](@ref) for a function with normalization.
"""
function crosscov(x::RealVector, y::RealVector, lags::IntegerVector; demean::Bool=true)
    out = Vector{float(Base.promote_eltype(x, y))}(undef, length(lags))
    crosscov!(out, x, y, lags; demean=demean)
end

function crosscov(x::RealMatrix, y::RealVector, lags::IntegerVector; demean::Bool=true)
    out = Matrix{float(Base.promote_eltype(x, y))}(undef, length(lags), size(x,2))
    crosscov!(out, x, y, lags; demean=demean)
end

function crosscov(x::RealVector, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    out = Matrix{float(Base.promote_eltype(x, y))}(undef, length(lags), size(y,2))
    crosscov!(out, x, y, lags; demean=demean)
end

function crosscov(x::RealMatrix, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    out = Array{float(Base.promote_eltype(x, y)),3}(undef, length(lags), size(x,2), size(y,2))
    crosscov!(out, x, y, lags; demean=demean)
end

crosscov(x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real}; demean::Bool=true) =
    crosscov(x, y, default_crosslags(size(x,1)); demean=demean)


## crosscor
"""
    crosscor!(r, x, y, lags; demean=true)

Compute the cross correlation between real-valued vectors or matrices `x` and `y` at
`lags` and store the result in `r`. `demean` specifies whether the respective means of
`x` and `y` should be subtracted from them before computing their cross correlation.

If both `x` and `y` are vectors, `r` must be a vector of the same length as
`lags`. If either `x` is a matrix and `y` is a vector, `r` must be a matrix of size
`(length(lags), size(x, 2))`; if `x` is a vector and `y` is a matrix, `r` must be a matrix
of size `(length(lags), size(y, 2))`. If both `x` and `y` are matrices, `r` must be a
three-dimensional array of size `(length(lags), size(x, 2), size(y, 2))`.

The output is normalized by `sqrt(var(x)*var(y))`. See [`crosscov!`](@ref) for the
unnormalized form.
"""
function crosscor!(r::RealVector, x::RealVector, y::RealVector, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    (length(y) == lx && length(r) == m) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    zx::Vector{T} = demean ? x .- mean(x) : x
    S = typeof(zero(eltype(y)) / 1)
    zy::Vector{S} = demean ? y .- mean(y) : y
    sc = sqrt(dot(zx, zx) * dot(zy, zy))
    for k = 1 : m  # foreach lag value
        r[k] = _crossdot(zx, zy, lx, lags[k]) / sc
    end
    return r
end

function crosscor!(r::RealMatrix, x::RealMatrix, y::RealVector, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    (length(y) == lx && size(r) == (m, ns)) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    zx = Vector{T}(undef, lx)
    S = typeof(zero(eltype(y)) / 1)
    zy::Vector{S} = demean ? y .- mean(y) : y
    yy = dot(zy, zy)
    for j = 1 : ns
        demean_col!(zx, x, j, demean)
        sc = sqrt(dot(zx, zx) * yy)
        for k = 1 : m
            r[k,j] = _crossdot(zx, zy, lx, lags[k]) / sc
        end
    end
    return r
end

function crosscor!(r::RealMatrix, x::RealVector, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    ns = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, ns)) || throw(DimensionMismatch())
    check_lags(lx, lags)

    T = typeof(zero(eltype(x)) / 1)
    zx::Vector{T} = demean ? x .- mean(x) : x
    S = typeof(zero(eltype(y)) / 1)
    zy = Vector{S}(undef, lx)
    xx = dot(zx, zx)
    for j = 1 : ns
        demean_col!(zy, y, j, demean)
        sc = sqrt(xx * dot(zy, zy))
        for k = 1 : m
            r[k,j] = _crossdot(zx, zy, lx, lags[k]) / sc
        end
    end
    return r
end

function crosscor!(r::AbstractArray{<:Real,3}, x::RealMatrix, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    nx = size(x, 2)
    ny = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, nx, ny)) || throw(DimensionMismatch())
    check_lags(lx, lags)

    # cached (centered) columns of x
    T = typeof(zero(eltype(x)) / 1)
    zxs = Vector{T}[]
    sizehint!(zxs, nx)
    xxs = Vector{T}(undef, nx)

    for j = 1 : nx
        xj = x[:,j]
        if demean
            mv = mean(xj)
            for i = 1 : lx
                xj[i] -= mv
            end
        end
        push!(zxs, xj)
        xxs[j] = dot(xj, xj)
    end

    S = typeof(zero(eltype(y)) / 1)
    zy = Vector{S}(undef, lx)
    for j = 1 : ny
        demean_col!(zy, y, j, demean)
        yy = dot(zy, zy)
        for i = 1 : nx
            zx = zxs[i]
            sc = sqrt(xxs[i] * yy)
            for k = 1 : m
                r[k,i,j] = _crossdot(zx, zy, lx, lags[k]) / sc
            end
        end
    end
    return r
end


"""
    crosscor(x, y, [lags]; demean=true)

Compute the cross correlation between real-valued vectors or matrices `x` and `y`,
optionally specifying the `lags`. `demean` specifies whether the respective means of
`x` and `y` should be subtracted from them before computing their cross correlation.

If both `x` and `y` are vectors, return a vector of the same length as
`lags`. Otherwise, compute cross covariances between each pairs of columns in `x` and `y`.

When left unspecified, the lags used are the integers from
`-min(size(x,1)-1, 10*log10(size(x,1)))` to `min(size(x,1), 10*log10(size(x,1)))`.

The output is normalized by `sqrt(var(x)*var(y))`. See [`crosscov`](@ref) for the
unnormalized form.
"""
function crosscor(x::RealVector, y::RealVector, lags::IntegerVector; demean::Bool=true)
    out = Vector{float(Base.promote_eltype(x, y))}(undef, length(lags))
    crosscor!(out, x, y, lags; demean=demean)
end

function crosscor(x::RealMatrix, y::RealVector, lags::IntegerVector; demean::Bool=true)
    out = Matrix{float(Base.promote_eltype(x, y))}(undef, length(lags), size(x,2))
    crosscor!(out, x, y, lags; demean=demean)
end

function crosscor(x::RealVector, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    out = Matrix{float(Base.promote_eltype(x, y))}(undef, length(lags), size(y,2))
    crosscor!(out, x, y, lags; demean=demean)
end

function crosscor(x::RealMatrix, y::RealMatrix, lags::IntegerVector; demean::Bool=true)
    out = Array{float(Base.promote_eltype(x, y)),3}(undef, length(lags), size(x,2), size(y,2))
    crosscor!(out, x, y, lags; demean=demean)
end

crosscor(x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real}; demean::Bool=true) =
    crosscor(x, y, default_crosslags(size(x,1)); demean=demean)


#######################################
#
#   Partial auto-correlations
#
#   TODO: the codes below need cleanup.
#
#######################################

function pacf_regress!(r::RealMatrix, X::RealMatrix, lags::IntegerVector, mk::Integer)
    lx = size(X, 1)
    tmpX = ones(eltype(X), lx, mk + 1)
    for j = 1 : size(X,2)
        for l = 1 : mk
            for i = 1+l:lx
                tmpX[i,l+1] = X[i-l,j]
            end
        end
        for i = 1 : length(lags)
            l = lags[i]
            sX = view(tmpX, 1+l:lx, 1:l+1)
            r[i,j] = l == 0 ? 1 : (cholesky!(sX'sX, Val(false)) \ (sX'view(X, 1+l:lx, j)))[end]
        end
    end
    r
end

function pacf_yulewalker!(r::RealMatrix, X::AbstractMatrix{T}, lags::IntegerVector, mk::Integer) where T<:RealFP
    tmp = Vector{T}(undef, mk)
    for j = 1 : size(X,2)
        acfs = autocor(X[:,j], 1:mk)
        for i = 1 : length(lags)
            l = lags[i]
            r[i,j] = l == 0 ? 1 : l == 1 ? acfs[i] : -durbin!(view(acfs, 1:l), tmp)[l]
        end
    end
end


"""
    pacf!(r, X, lags; method=:regression)

Compute the partial autocorrelation function (PACF) of a matrix `X` at `lags` and
store the result in `r`. `method` designates the estimation method. Recognized values
are `:regression`, which computes the partial autocorrelations via successive
regression models, and `:yulewalker`, which computes the partial autocorrelations
using the Yule-Walker equations.

`r` must be a matrix of size `(length(lags), size(x, 2))`.
"""
function pacf!(r::RealMatrix, X::AbstractMatrix{T}, lags::IntegerVector; method::Symbol=:regression) where T<:RealFP
    lx = size(X, 1)
    m = length(lags)
    minlag, maxlag = extrema(lags)
    (0 <= minlag && 2maxlag < lx) || error("Invalid lag value.")
    size(r) == (m, size(X,2)) || throw(DimensionMismatch())

    if method == :regression
        pacf_regress!(r, X, lags, maxlag)
    elseif method == :yulewalker
        pacf_yulewalker!(r, X, lags, maxlag)
    else
        error("Invalid method: $method")
    end
    return r
end


"""
    pacf(X, lags; method=:regression)

Compute the partial autocorrelation function (PACF) of a real-valued vector
or matrix `X` at `lags`. `method` designates the estimation method. Recognized
values are `:regression`, which computes the partial autocorrelations via successive
regression models, and `:yulewalker`, which computes the partial autocorrelations
using the Yule-Walker equations.

If `x` is a vector, return a vector of the same length as `lags`.
If `x` is a matrix, return a matrix of size `(length(lags), size(x, 2))`,
where each column in the result corresponds to a column in `x`.
"""
function pacf(X::RealMatrix, lags::IntegerVector; method::Symbol=:regression)
    out = Matrix{float(eltype(X))}(undef, length(lags), size(X,2))
    pacf!(out, float(X), lags; method=method)
end

function pacf(x::RealVector, lags::IntegerVector; method::Symbol=:regression)
    vec(pacf(reshape(x, length(x), 1), lags, method=method))
end
