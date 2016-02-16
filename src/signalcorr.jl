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

function demean_col!{T<:RealFP}(z::Vector{T}, x::Matrix{T}, j::Int, demean::Bool)
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
        copy!(z, 1, x, b+1, m)
    end
    z
end


#######################################
#
#   Auto-correlations
#
#######################################

default_autolags(lx::Int) = 0 : default_laglen(lx)

_autodot{T<:RealFP}(x::Vector{T}, lx::Int, l::Int) = dot(x, 1:lx-l, x, 1+l:lx)


## autocov

function autocov!{T<:RealFP}(r::RealVector, x::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    length(r) == m || raise_dimerror()
    check_lags(lx, lags)

    z::Vector{T} = demean ? x .- mean(x) : x
    for k = 1 : m  # foreach lag value
        r[k] = _autodot(z, lx, lags[k]) / lx
    end
    return r
end

function autocov!{T<:RealFP}(r::RealMatrix, x::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    size(r) == (m, ns) || raise_dimerror()
    check_lags(lx, lags)

    z = Array(T, lx)
    for j = 1 : ns
        demean_col!(z, x, j, demean)
        for k = 1 : m
            r[k,j] = _autodot(z, lx, lags[k]) / lx
        end
    end
    return r
end

function autocov{T<:Real}(x::Vector{T}, lags::IntegerVector; demean::Bool=true)
    autocov!(Array(fptype(T), length(lags)), float(x), lags; demean=demean)
end

function autocov{T<:Real}(x::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    autocov!(Array(fptype(T), length(lags), size(x,2)), float(x), lags; demean=demean)
end

autocov{T<:Real}(x::VecOrMat{T}; demean::Bool=true) = autocov(x, default_autolags(size(x,1)); demean=demean)

"""
    autocov(x, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the autocovariance of. 
* `lags`: Specifies the lags at which the autocovariance of `x` is copmuted.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute auto-covariance of `x` at specified lags. If `x` is a vector, it returns a vector of the same length of lags. If `x` is a matrix, it returns a matrix of size `(length(lags), size(x,2))`, where each column in the result corresponding to a column in `x`.

Here, demean is a keyword argument (default value is true), which means that the function will subtract each `x` from its mean before computing the results. Otherwise, `x` is considered as having been centered.
"""
function autocov end

"""
    autocov!(r, x, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the autocovariance of. 
* `lags`: Specifies the lags at which the autocovariance of `x` is copmuted.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute auto-covariance of `x` at specified lags. If `x` is a vector, it returns a vector of the same length of lags. If `x` is a matrix, it returns a matrix of size `(length(lags), size(x,2))`, where each column in the result corresponding to a column in `x`.

Here, demean is a keyword argument (default value is true), which means that the function will subtract each `x` from its mean before computing the results. Otherwise, `x` is considered as having been centered.

the computed autocovariance is written back to `r`.
"""
function autocov end

## autocor

function autocor!{T<:RealFP}(r::RealVector, x::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    length(r) == m || raise_dimerror()
    check_lags(lx, lags)

    z::Vector{T} = demean ? x .- mean(x) : x
    zz = dot(z, z)
    for k = 1 : m  # foreach lag value
        r[k] = _autodot(z, lx, lags[k]) / zz
    end
    return r
end

function autocor!{T<:RealFP}(r::RealMatrix, x::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    size(r) == (m, ns) || raise_dimerror()
    check_lags(lx, lags)

    z = Array(T, lx)
    for j = 1 : ns
        demean_col!(z, x, j, demean)
        zz = dot(z, z)
        for k = 1 : m
            r[k,j] = _autodot(z, lx, lags[k]) / zz
        end
    end
    return r
end

function autocor{T<:Real}(x::Vector{T}, lags::IntegerVector; demean::Bool=true)
    autocor!(Array(fptype(T), length(lags)), float(x), lags; demean=demean)
end

function autocor{T<:Real}(x::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    autocor!(Array(fptype(T), length(lags), size(x,2)), float(x), lags; demean=demean)
end

autocor{T<:Real}(x::VecOrMat{T}; demean::Bool=true) = autocor(x, default_autolags(size(x,1)); demean=demean)

"""
    autocor(x, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the auto-correlation of. 
* `lags`: Specifies the lags at which the auto-correlation of `x` is copmuted.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute auto-correlation of `x` at specified lags. If `x` is a vector, it returns a vector of the same length of lags. If `x` is a matrix, it returns a matrix of size `(length(lags), size(x,2))`, where each column in the result corresponding to a column in `x`.

Here, demean is a keyword argument (default value is true), which means that the function will subtract each `x` from its mean before computing the results. Otherwise, `x` is considered as having been centered.
"""
function autocor end

"""
    autocor!(r, x, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the autocorrelation of. 
* `lags`: Specifies the lags at which the autocorrelation of `x` is copmuted.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute auto-correlation of `x` at specified lags. If `x` is a vector, it returns a vector of the same length of lags. If `x` is a matrix, it returns a matrix of size `(length(lags), size(x,2))`, where each column in the result corresponding to a column in `x`.

Here, demean is a keyword argument (default value is true), which means that the function will subtract each `x` from its mean before computing the results. Otherwise, `x` is considered as having been centered.

the computed auto-correlation is written back to `r`.
"""
function autocor! end

#######################################
#
#   Cross-correlations
#
#######################################

default_crosslags(lx::Int) = (l=default_laglen(lx); -l:l)

_crossdot{T<:RealFP}(x::Vector{T}, y::Vector{T}, lx::Int, l::Int) = (l >= 0 ? dot(x, 1:lx-l, y, 1+l:lx) : dot(x, 1-l:lx, y, 1:lx+l))

## crosscov

function crosscov!{T<:RealFP}(r::RealVector, x::Vector{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    (length(y) == lx && length(r) == m) || raise_dimerror()
    check_lags(lx, lags)

    zx::Vector{T} = demean ? x .- mean(x) : x
    zy::Vector{T} = demean ? y .- mean(y) : y
    for k = 1 : m  # foreach lag value
        r[k] = _crossdot(zx, zy, lx, lags[k]) / lx
    end
    return r
end

function crosscov!{T<:RealFP}(r::RealMatrix, x::Matrix{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    (length(y) == lx && size(r) == (m, ns)) || raise_dimerror()
    check_lags(lx, lags)

    zx = Array(T, lx)
    zy::Vector{T} = demean ? y .- mean(y) : y
    for j = 1 : ns
        demean_col!(zx, x, j, demean)
        for k = 1 : m
            r[k,j] = _crossdot(zx, zy, lx, lags[k]) / lx
        end
    end
    return r
end

function crosscov!{T<:RealFP}(r::RealMatrix, x::Vector{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    ns = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, ns)) || raise_dimerror()
    check_lags(lx, lags)

    zx::Vector{T} = demean ? x .- mean(x) : x
    zy = Array(T, lx)
    for j = 1 : ns
        demean_col!(zy, y, j, demean)
        for k = 1 : m
            r[k,j] = _crossdot(zx, zy, lx, lags[k]) / lx
        end
    end
    return r
end

function crosscov!{T<:RealFP}(r::AbstractArray{T,3}, x::Matrix{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    nx = size(x, 2)
    ny = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, nx, ny)) || raise_dimerror()
    check_lags(lx, lags)

    # cached (centered) columns of x
    zxs = Array(Vector{T}, 0)
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

    zx = Array(T, lx)
    zy = Array(T, lx)
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

function crosscov{T<:Real}(x::Vector{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    crosscov!(Array(fptype(T), length(lags)), float(x), float(y), lags; demean=demean)
end

function crosscov{T<:Real}(x::Matrix{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    crosscov!(Array(fptype(T), length(lags), size(x,2)), float(x), float(y), lags; demean=demean)
end

function crosscov{T<:Real}(x::Vector{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    crosscov!(Array(fptype(T), length(lags), size(y,2)), float(x), float(y), lags; demean=demean)
end

function crosscov{T<:Real}(x::Matrix{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    crosscov!(Array(fptype(T), length(lags), size(x,2), size(y,2)), float(x), float(y), lags; demean=demean)
end

crosscov{T<:Real}(x::VecOrMat{T}, y::VecOrMat{T}; demean::Bool=true) = crosscov(x, y, default_crosslags(size(x,1)); demean=demean)

"""
    crosscov(x, y, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `y`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `lags`: Specifies the lags at which the cross-covariance is computed.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute cross-covariance of `x` at specified lags. If both `x` and `y` are vectors, it returns a vector of th same length of lags. Otherwise, it computes cross covariances between each pairs of columns in x and y.
"""
function crosscov end    

"""
    crosscov!(r, x, y, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `y`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `lags`: Specifies the lags at which the cross-covariance is copmuted.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute cross-covariance of `x` at specified lags. If both `x` and `y` are vectors, it returns a vector of th same length of lags. Otherwise, it computes cross covariances between each pairs of columns in `x` and `y`.

The computed cross-covariance is written back to `r`.
"""
function crosscov end

## crosscor

function crosscor!{T<:RealFP}(r::RealVector, x::Vector{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    (length(y) == lx && length(r) == m) || raise_dimerror()
    check_lags(lx, lags)

    zx::Vector{T} = demean ? x .- mean(x) : x
    zy::Vector{T} = demean ? y .- mean(y) : y
    sc = sqrt(dot(zx, zx) * dot(zy, zy))
    for k = 1 : m  # foreach lag value
        r[k] = _crossdot(zx, zy, lx, lags[k]) / sc
    end
    return r
end

function crosscor!{T<:RealFP}(r::RealMatrix, x::Matrix{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    ns = size(x, 2)
    m = length(lags)
    (length(y) == lx && size(r) == (m, ns)) || raise_dimerror()
    check_lags(lx, lags)

    zx = Array(T, lx)
    zy::Vector{T} = demean ? y .- mean(y) : y
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

function crosscor!{T<:RealFP}(r::RealMatrix, x::Vector{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    ns = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, ns)) || raise_dimerror()
    check_lags(lx, lags)

    zx::Vector{T} = demean ? x .- mean(x) : x
    zy = Array(T, lx)
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

function crosscor!{T<:RealFP}(r::AbstractArray{T,3}, x::Matrix{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    lx = size(x, 1)
    nx = size(x, 2)
    ny = size(y, 2)
    m = length(lags)
    (size(y, 1) == lx && size(r) == (m, nx, ny)) || raise_dimerror()
    check_lags(lx, lags)

    # cached (centered) columns of x
    zxs = Array(Vector{T}, 0)
    sizehint!(zxs, nx)
    xxs = Array(T, nx)

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

    zx = Array(T, lx)
    zy = Array(T, lx)
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

function crosscor{T<:Real}(x::Vector{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    crosscor!(Array(fptype(T), length(lags)), float(x), float(y), lags; demean=demean)
end

function crosscor{T<:Real}(x::Matrix{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    crosscor!(Array(fptype(T), length(lags), size(x,2)), float(x), float(y), lags; demean=demean)
end

function crosscor{T<:Real}(x::Vector{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    crosscor!(Array(fptype(T), length(lags), size(y,2)), float(x), float(y), lags; demean=demean)
end

function crosscor{T<:Real}(x::Matrix{T}, y::Matrix{T}, lags::IntegerVector; demean::Bool=true)
    crosscor!(Array(fptype(T), length(lags), size(x,2), size(y,2)), float(x), float(y), lags; demean=demean)
end

crosscor{T<:Real}(x::VecOrMat{T}, y::VecOrMat{T}; demean::Bool=true) = crosscor(x, y, default_crosslags(size(x,1)); demean=demean)

"""
    crosscov(x, y, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `y`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `lags`: Specifies the lags at which the cross-covariance is computed.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute cross-covariance of `x` at specified lags. If both `x` and `y` are vectors, it returns a vector of th same length of lags. Otherwise, it computes cross covariances between each pairs of columns in x and y.
"""
function crosscov end    
    
"""
    crosscov!(r, x, y, lags, demean)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `y`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `lags`: Specifies the lags at which the cross-covariance is copmuted.
* `demean`: A `Bool` type variable, this is a named paramter. 

Compute cross-covariance of `x` at specified lags. If both `x` and `y` are vectors, it returns a vector of th same length of lags. Otherwise, it computes cross covariances between each pairs of columns in `x` and `y`.

The computed cross-covariance is written back to `r`.
"""
function crosscov! end    
#######################################
#
#   Partial auto-correlations
#
#   TODO: the codes below need cleanup.
#
#######################################

function pacf_regress!{T<:RealFP}(r::RealMatrix, X::AbstractMatrix{T}, lags::IntegerVector, mk::Integer)
    lx = size(X, 1)
    tmpX = ones(T, lx, mk + 1)
    for j = 1 : size(X,2)
        for l = 1 : mk
            for i = 1+l:lx
                tmpX[i,l+1] = X[i-l,j]
            end
        end
        for i = 1 : length(lags)
            l = lags[i]
            sX = sub(tmpX, 1+l:lx, 1:l+1)
            r[i,j] = (cholfact!(sX'sX)\(sX'sub(X, 1+l:lx, j)))[end]
        end
    end 
    r   
end

function pacf_yulewalker!{T<:RealFP}(r::RealMatrix, X::Matrix{T}, lags::IntegerVector, mk::Integer)
    tmp = Array(T, mk)
    for j = 1 : size(X,2)
        acfs = autocor(X[:,j], 1:mk)
        for i = 1 : length(lags)
            l = lags[i]
            r[i,j] = l == 0 ? one(T) : l == 1 ? acfs[i] : -durbin!(sub(acfs, 1:l), tmp)[l]
        end
    end
end

function pacf!{T<:RealFP}(r::RealMatrix, X::Matrix{T}, lags::IntegerVector; method::Symbol=:regression)
    lx = size(X, 1)
    m = length(lags)
    minlag, maxlag = extrema(lags)
    (0 <= minlag && 2maxlag < lx) || error("Invalid lag value.")
    size(r) == (m, size(X,2)) || raise_dimerror()

    if method == :regression
        pacf_regress!(r, X, lags, maxlag)        
    elseif method == :yulewalker
        pacf_yulewalker!(r, X, lags, maxlag)
    else
        error("Invalid method: $method")
    end
    return r
end

function pacf{T<:Real}(X::Matrix{T}, lags::IntegerVector; method::Symbol=:regression)
    pacf!(Array(fptype(T), length(lags), size(X,2)), float(X), lags; method=method)
end

function pacf{T<:Real}(x::Vector{T}, lags::IntegerVector; method::Symbol=:regression)
    vec(pacf(reshape(x, length(x), 1), lags, method=method))
end

"""
    pacf(x, lags; method=:regression)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `lags`: Specifies the lags at which the cross-covariance is computed.
* `method`: A named parameter, which specifies the choice of algorithm, which can be either :regresion or :yulewalker. The default value is :regression.

Compute partial auto-correlation of `x` at specified lags. If `x` is a vector, it returns a vector of th same length of lags. If `x` is a matrix, it returns a matrix of size `(length(lags), size(x,2))`, where each column in the result corresponding to a column in `x`.
"""

"""
    pacf!(x, lags; method=:regression)

### Args:
* `x`: Input `Vector` or `Matrix` to compute the cross-covariance of.
* `lags`: Specifies the lags at which the cross-covariance is computed.
* `method`: A named parameter, which specifies the choice of algorithm, which can be either :regresion or :yulewalker. The default value is :regression.

Compute partial auto-correlation of `x` at specified lags. If `x` is a vector, it returns a vector of th same length of lags. If `x` is a matrix, it returns a matrix of size `(length(lags), size(x,2))`, where each column in the result corresponding to a column in `x`.

Write the computed partial auto-correlation to r.
"""

