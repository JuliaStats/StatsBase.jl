# Various correlation 

#
# spearman correlation functions
#

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

# autocorrelation for range
function autocor(x::AbstractVector, lags::Ranges)
    lx = length(x)
    if max(lags) > lx error("Autocorrelation distance must be less than sample size") end
    mx = mean(x)
    sxinv = 1/stdm(x, mx)
    xs = Array(typeof(sxinv), lx)
    for i = 1:lx
        xs[i] = (x[i] - mx)*sxinv
    end
    acf = Array(typeof(sxinv), length(lags))
    for i in 1:length(lags)
        acf[i] = dot(xs[1:end - lags[i]], xs[lags[i] + 1:end])/(lx - 1)
    end
    return acf
end

# autocorrelation at a specific lag
autocor(x::AbstractVector, lags::Real) = autocor(x, lags:lags)[1]

# autocorrelation at a default of zero to 10log10(length(v)) lags
autocor(v::AbstractVector) = autocor(v, 0:min(length(v) - 1, 10log10(length(v))))

# crosscorrelation between two vectors for range
function autocor(x::AbstractVector, y::AbstractVector, lags::Ranges)
  lx, ly = length(x), length(y)

  if lx != ly error("Input vectors must have same length") end
  if max(lags) > lx error("Autocorrelation distance must be less than sample size") end

  mx, my = mean(x), mean(y)
  sxinv, syinv = 1/stdm(x, mx), 1/stdm(y, my)
  xs, ys = Array(typeof(sxinv), lx), Array(typeof(syinv), ly)
  for i = 1:lx
    xs[i], ys[i] = (x[i] - mx)*sxinv, (y[i] - my)*syinv
  end
  acf = Array(typeof(sxinv), length(lags))
  for i in 1:length(lags)
    acf[i] = dot(xs[1:end - lags[i]], ys[lags[i] + 1:end])/(lx - 1)
  end
  return acf
end

# crosscorrelation between two vectors at a specific lag
autocor(x::AbstractVector, y::AbstractVector, lags::Real) = autocor(x, y, lags:lags)[1]

# crosscorrelation between two vectors at a default of zero to 10log10(length(v)) lags
autocor(x::AbstractVector, y::AbstractVector) = autocor(x, y, 0:min(length(x) - 1, 10log10(length(x))))

# crosscorrelation between all pairs of columns of a matrix for range
function autocor(x::AbstractMatrix, lags::Ranges)
  ncols = size(x, 2)
  acf = Array(eltype(x), length(lags), ncols, ncols)

  for i = 1:ncols
    for j = 1:ncols
      acf[:, i, j] = autocor(x[:, i], x[:, j], lags)
    end
  end

  return acf
end

# crosscorrelation between all pairs of columns of a matrix at a specific lag
autocor(x::AbstractMatrix, lags::Real) = reshape(autocor(x, lags:lags), size(x, 2), size(x, 2))

# crosscorrelation between all pairs of columns of a matrix at a default of zero to 10log10(length(v)) lags
autocor(x::AbstractMatrix) = autocor(x, 0:min(length(x) - 1, 10log10(length(x))))

# Unlike autocor, compute only autocorrelations (not crosscorrelations) between matrix columns for range
function autocordiag(x::AbstractMatrix, lags::Ranges)
  ncols = size(x, 2)
  acf = Array(eltype(x), length(lags), ncols)

  for i = 1:ncols
    acf[:, i] = autocor(x[:, i], lags)
  end

  return acf
end

# Unlike autocor, compute only autocorrelations (not crosscorrelations) between matrix columns at a specific lag
autocordiag(x::AbstractMatrix, lags::Real) = autocordiag(x, lags:lags)

# Unlike autocor, compute only autocorrelations (not crosscorrelations) between matrix columns at a default of zero to 10log10(length(v)) lags
autocordiag(x::AbstractMatrix) = autocordiag(x, 0:min(length(x) - 1, 10log10(length(x))))
