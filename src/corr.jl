# Correlation 

#######################################
#
#   Helper functions
#
#######################################

default_laglen(lx::Int) = min(lx-1, int(10log10(lx)))
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

    z::Vector{T} = demean ? x - mean(x) : x
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

## autocor

function autocor!{T<:RealFP}(r::RealVector, x::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    length(r) == m || raise_dimerror()
    check_lags(lx, lags)

    z::Vector{T} = demean ? x - mean(x) : x
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

const acf = autocor



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

    zx::Vector{T} = demean ? x - mean(x) : x
    zy::Vector{T} = demean ? y - mean(y) : y
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
    zy::Vector{T} = demean ? y - mean(y) : y
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

    zx::Vector{T} = demean ? x - mean(x) : x
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
    sizehint(zxs, nx)
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


## crosscor

function crosscor!{T<:RealFP}(r::RealVector, x::Vector{T}, y::Vector{T}, lags::IntegerVector; demean::Bool=true)
    lx = length(x)
    m = length(lags)
    (length(y) == lx && length(r) == m) || raise_dimerror()
    check_lags(lx, lags)

    zx::Vector{T} = demean ? x - mean(x) : x
    zy::Vector{T} = demean ? y - mean(y) : y
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
    zy::Vector{T} = demean ? y - mean(y) : y
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

    zx::Vector{T} = demean ? x - mean(x) : x
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
    sizehint(zxs, nx)
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

ccf = crosscor



#######################################
#
#   Spearman correlation
#
#######################################

# spearman correlation between two vectors
function cor_spearman(x::AbstractVector, y::AbstractVector)
    if any(isnan(x)) || any(isnan(y)) return NaN end
    return cor(tiedrank(x), tiedrank(y))
end

# spearman correlation over all pairs of columns of two matrices
function cor_spearman(X::AbstractMatrix, Y::AbstractMatrix)
    return cor(mapslices(tiedrank, X, 1), mapslices(tiedrank, Y, 1))
end
function cor_spearman(X::AbstractMatrix, y::AbstractVector)
    return cor(mapslices(tiedrank, X, 1), tiedrank(y))
end
function cor_spearman(x::AbstractVector, Y::AbstractMatrix)
    return cor(tiedrank(x), mapslices(tiedrank, Y, 1))
end

# spearman correlation over all pairs of columns of a matrix
function cor_spearman(X::AbstractMatrix)
    csp = cor(mapslices(tiedrank, X, 1))
    nanindex = vec(mapslices(any, isnan(X), 1))
    csp[nanindex, :] = NaN
    csp[:, nanindex] = NaN
    return csp
end

#
# Kendall's rank correlation
#

# Knigh JASA (1966)
function cor_kendall!{T<:Real,S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
    if any(isnan(x)) || any(isnan(y)) return NaN end
    n = length(x)
    if n != length(y) error("Vectors must have same length") end
    
    # Initial sorting
    pm = sortperm(y)
    x[:] = x[pm]
    y[:] = y[pm]
    pm[:] = sortperm(x)
    x[:] = x[pm]

    # Counting ties in x and y
    iT = 1
    nT = 0
    iU = 1
    nU = 0
    for i = 2:n
        if x[i] == x[i-1] 
            iT += 1
        else
            nT += iT*(iT - 1)
            iT = 1
        end
        if y[i] == y[i-1]
            iU += 1
        else
            nU += iU*(iU - 1)
            iU = 1
        end
    end
    if iT > 1 nT += iT*(iT - 1) end
    nT = div(nT,2)
    if iU > 1 nU += iU*(iU - 1) end
    nU = div(nU,2)

    # Sort y after x
    y[:] = y[pm]

    # Calculate double ties
    iV = 1
    nV = 0
    jV = 1
    for i = 2:n
        if x[i] == x[i-1] && y[i] == y[i-1] 
            iV += 1
        else
            nV += iV*(iV - 1)
            iV = 1
        end
    end
    if iV > 1 nV += iV*(iV - 1) end
    nV = div(nV,2)

    nD = div(n*(n - 1),2)
    return (nD - nT - nU + nV - 2swaps!(y))/sqrt((nD - nT)*(nD - nU))
end
cor_kendall(x::AbstractVector, y::AbstractVector) = cor_kendall!(copy(x), copy(y))
cor_kendall(x::AbstractVector, Y::AbstractMatrix) = [cor_kendall(x, Y[:,i]) for i in 1:size(Y, 2)]
cor_kendall(X::AbstractMatrix, y::AbstractVector) = [cor_kendall(X[:,i], y) for i in 1:size(X, 2)]
cor_kendall(X::AbstractMatrix, Y::AbstractMatrix) = [cor_kendall(X[:,i], Y[:,j]) for i in 1:size(X, 2), j in 1:size(Y, 2)]
function cor_kendall(X::AbstractMatrix)
    n = size(X, 2)
    C = eye(n)
    for j = 2:n
        for i = 1:j-1
            C[i,j] = cor_kendall!(X[:,i],X[:,j])
            C[j,i] = C[i,j]
        end
    end
    return C
end

# Auxilliary functions for Kendall's rank correlation
function swaps!(x::AbstractVector)
    n = length(x)
    if n == 1 return 0 end
    n2 = div(n, 2)
    xl = sub(x, 1:n2)
    xr = sub(x, n2+1:n)
    nsl = swaps!(xl)
    nsr = swaps!(xr)
    sort!(xl)
    sort!(xr)
    return nsl + nsr + mswaps(xl,xr)
end

function mswaps(x::AbstractVector, y::AbstractVector)
    i = 1
    j = 1
    nSwaps = 0
    n = length(x)
    while i <= n && j <= length(y)
        if y[j] < x[i]
            nSwaps += n - i + 1
            j += 1
        else
            i += 1
        end
    end
    return nSwaps
end


# Partial autoroccelation
function pacf{T<:BlasReal}(X::AbstractMatrix{T}, lags::AbstractVector{Int} = 0:min(size(X,1)-1, int(10log10(size(X,1)))); 
    method::Symbol = :regression)
    n, p = size(X)
    nk = length(lags)
    mk = maximum(lags)
    if minimum(lags) < 0 error("Negative autoroccelations not allowed") end
    if 2mk >= n error("Can at most calculate pacf for $(div(n,2) - 1) lags, you requested $mk") end
    val = Array(T, nk, p)
    if method == :regression
        tmpX = ones(T, n, mk + 1)
        for j = 1:p
            for l = 1:mk
                for i = 1+l:n
                    tmpX[i,l+1] = X[i-l,j]
                end
            end
            i = 1
            for l in lags
                sX = sub(tmpX, 1+l:n, 1:l+1)
                val[i,j] = (cholfact!(sX'sX)\(sX'sub(X, 1+l:n, j)))[end]
                i += 1
            end
        end
    elseif method == :yulewalker
        tmp = Array(T, mk)
        for j = 1:p
            acfs = acf(sub(X,1:n,j),1:mk)
            i = 1
            for l in lags
                if l == 0
                    val[i,j] = one(T)
                elseif l == 1
                    val[i,j] = acfs[i]
                else
                    val[i,j] = -durbin!(sub(acfs, 1:l), tmp)[l]
                end
                i += 1
            end
        end
    else
        error("No such method")
    end
    return val
end
pacf{T<:Real}(X::AbstractMatrix{T}, args1...; args2...) = pacf(float(X), args1...; args2...)
pacf(x::AbstractVector, args1...; args2...) = squeeze(pacf(reshape(x, length(x), 1), args1...; args2...),2)
