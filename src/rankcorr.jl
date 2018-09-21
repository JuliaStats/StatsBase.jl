# Rank-based correlations
#
# - Spearman's correlation
# - Kendall's correlation
#

#######################################
#
#   Spearman correlation
#
#######################################

"""
    corspearman(x, y=x)

Compute Spearman's rank correlation coefficient. If `x` and `y` are vectors, the
output is a float, otherwise it's a matrix corresponding to the pairwise correlations
of the columns of `x` and `y`.
"""
corspearman(x::RealVector, y::RealVector) = cor(tiedrank(x), tiedrank(y))

corspearman(X::RealMatrix, Y::RealMatrix) =
    cor(mapslices(tiedrank, X, dims=1), mapslices(tiedrank, Y, dims=1))
corspearman(X::RealMatrix, y::RealVector) = cor(mapslices(tiedrank, X, dims=1), tiedrank(y))
corspearman(x::RealVector, Y::RealMatrix) = cor(tiedrank(x), mapslices(tiedrank, Y, dims=1))

corspearman(X::RealMatrix) = (Z = mapslices(tiedrank, X, dims=1); cor(Z, Z))


#######################################
#
#   Kendall correlation
#
#######################################

# Knight JASA (1966)

function corkendall!(x::RealVector, y::RealVector)
    if any(isnan, x) || any(isnan, y) return NaN end
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
    return (nD - nT - nU + nV - 2swaps!(y)) / (sqrt(nD - nT) * sqrt(nD - nU))
end


"""
    corkendall(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors.
"""
corkendall(x::RealVector, y::RealVector) = corkendall!(float(copy(x)), float(copy(y)))

corkendall(X::RealMatrix, y::RealVector) = Float64[corkendall!(float(X[:,i]), float(copy(y))) for i in 1:size(X, 2)]

corkendall(x::RealVector, Y::RealMatrix) = (n = size(Y,2); reshape(Float64[corkendall!(float(copy(x)), float(Y[:,i])) for i in 1:n], 1, n))

corkendall(X::RealMatrix, Y::RealMatrix) = Float64[corkendall!(float(X[:,i]), float(Y[:,j])) for i in 1:size(X, 2), j in 1:size(Y, 2)]

function corkendall(X::RealMatrix)
    n = size(X, 2)
    C = Matrix{eltype(X)}(I, n, n)
    for j = 2:n
        for i = 1:j-1
            C[i,j] = corkendall!(X[:,i],X[:,j])
            C[j,i] = C[i,j]
        end
    end
    return C
end

# Auxilliary functions for Kendall's rank correlation

function swaps!(x::RealVector)
    n = length(x)
    if n == 1 return 0 end
    n2 = div(n, 2)
    xl = view(x, 1:n2)
    xr = view(x, n2+1:n)
    nsl = swaps!(xl)
    nsr = swaps!(xr)
    sort!(xl)
    sort!(xr)
    return nsl + nsr + mswaps(xl,xr)
end

function mswaps(x::RealVector, y::RealVector)
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

