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

# Autocorrelation for range
function acf{T<:Real}(x::AbstractVector{T},
  lags::AbstractVector{Int}=convert(Vector{Int}, 0:min(length(x)-1, int(10log10(length(x)))));
  correlation::Bool=true, demean::Bool=true)
  lx, llags = length(x), length(lags)
  if max(lags) >= lx; error("Autocovariance distance must be less than sample size"); end

  xs = Array(T, lx)
  demean ? (mx = mean(x); for i = 1:lx; xs[i] = x[i]-mx; end) : (for i = 1:lx; xs[i] = x[i]; end)

  autocov_sumterm = Array(T, llags)
  for i in 1:llags
    autocov_sumterm[i] = dot(xs[1:end-lags[i]], xs[lags[i]+1:end])
  end
  autocov_sumterm

  if correlation
    demean ?
      (return autocov_sumterm/(min(lags)==0 ? autocov_sumterm[1] : dot(xs,xs))) : (return autocov_sumterm/dot(x,x))
  else
    return autocov_sumterm/lx
  end
end

# Autocorrelation at a specific lag
acf{T<:Real}(x::AbstractVector{T}, lags::Integer; correlation::Bool=true, demean::Bool=true) =
  acf(x, lags:lags, correlation=correlation, demean=demean)[1]

# Cross-correlation for range
function acf{T<:Real}(x::AbstractVector{T}, y::AbstractVector{T},
  lags::AbstractVector{Int}=convert(Vector{Int}, 0:min(length(x)-1, int(10log10(length(x)))));
  correlation::Bool=true, demean::Bool=true)
  lx, ly, llags = length(x), length(y), length(lags)
  if lx != ly error("Input vectors must have same length") end
  if max(lags) > lx; error("Cross-covariance distance must be less than sample size"); end

  xs, ys = Array(T, lx), Array(T, ly)
  if demean
    mx, my = mean(x), mean(y); for i = 1:lx; xs[i], ys[i] = x[i]-mx, y[i]-my; end
  else
    for i = 1:lx; xs[i], ys[i] = x[i], y[i]; end
  end

  crosscov_sumterm = Array(T, llags)
  for i in 1:llags
    crosscov_sumterm[i] = dot(xs[1:end-lags[i]], ys[lags[i]+1:end])
  end
  crosscov_sumterm

  if correlation
    demean ? (return crosscov_sumterm/((lx-1)*std(x)*std(y))) : (return crosscov_sumterm/sqrt(sum(x.^2)*sum(y.^2)))
  else
    return crosscov_sumterm/lx
  end  
end

# Cross-correlation at a specific lag
acf{T<:Real}(x::AbstractVector{T}, y::AbstractVector{T}, lags::Integer; correlation::Bool=true, demean::Bool=true) =
  acf(x, y, lags:lags, correlation=correlation, demean=demean)[1]

# Cross-correlation between all pairs of columns of a matrix for range
function acfall{T<:Real}(x::AbstractMatrix{T},
  lags::AbstractVector{Int}=convert(Vector{Int}, 0:min(length(x)-1, int(10log10(length(x)))));
  correlation::Bool=true, demean::Bool=true)
  ncols = size(x, 2)

  crosscorr = Array(T, length(lags), ncols, ncols)
  for i = 1:ncols
    for j = 1:ncols
      crosscorr[:, i, j] = acf(x[:, i], x[:, j], lags, correlation=correlation, demean=demean)
    end
  end
  crosscorr
end

# Cross-correlation between all pairs of columns of a matrix at a specific lag
acfall{T<:Real}(x::AbstractMatrix{T}, lags::Integer; correlation::Bool=true, demean::Bool=true) =
  reshape(acfall(x, lags:lags, correlation=correlation, demean=demean), size(x, 2), size(x, 2))

# Unlike acfall, compute only autocorrelation (not cross-correlation) of matrix columns for range
function acfdiag{T<:Real}(x::AbstractMatrix{T},
  lags::AbstractVector{Int}=convert(Vector{Int}, 0:min(length(x)-1, int(10log10(length(x)))));
  correlation::Bool=true, demean::Bool=true)
  ncols = size(x, 2)

  autocorr = Array(T, length(lags), ncols)
  for i = 1:ncols
    autocorr[:, i] = acf(x[:, i], lags, correlation=correlation, demean=demean)
  end
  autocorr
end

# Unlike acfall, compute only autocorrelation (not cross-correlation) of matrix columns at a specific lag
acfdiag{T<:Real}(x::AbstractMatrix{T}, lags::Integer; correlation::Bool=true, demean::Bool=true) =
  acfdiag(x, lags:lags, correlation=correlation, demean=demean)

# acf wrapper (with matrix as input) for range
function acf{T<:Real}(x::AbstractMatrix{T},
  lags::AbstractVector{Int}=convert(Vector{Int}, 0:min(length(x)-1, int(10log10(length(x)))));
  correlation::Bool=true, demean::Bool=true, diag::Bool=false)
  diag ?
    acfdiag(x, lags, correlation=correlation, demean=demean) :
    acfall(x, lags, correlation=correlation, demean=demean)
end

# acf wrapper (with matrix as input) at a specific lag
function acf{T<:Real}(x::AbstractMatrix{T}, lags::Integer; correlation::Bool=true, demean::Bool=true, diag::Bool=false)
 diag ?
   acfdiag(x, lags, correlation=correlation, demean=demean) :
   acfall(x, lags, correlation=correlation, demean=demean)
end
