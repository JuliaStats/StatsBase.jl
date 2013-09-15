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

# Autocovariance sum term for range
function acv_sumterm(x::AbstractVector, lags::Ranges; demean::Bool=true)
  lx, llags, typex = length(x), length(lags), eltype(x)
  if max(lags) > lx; error("Autocovariance distance must be less than sample size"); end

  xs = Array(typex, lx)
  demean ? (mx = mean(x); for i = 1:lx; xs[i] = x[i]-mx; end) : (for i = 1:lx; xs[i] = x[i]; end)

  autocov = Array(typex, llags)
  for i in 1:llags
    autocov[i] = dot(xs[1:end-lags[i]], xs[lags[i]+1:end])
  end
  autocov
end

# Autocovariance sum term at a specific lag
acv_sumterm(x::AbstractVector, lags::Real; demean::Bool=true) = acv_sumterm(x, lags:lags, demean=demean)[1]

# Autocovariance sum term at a default of zero to 10log10(length(x)) lags
acv_sumterm(x::AbstractVector; demean::Bool=true) =
  acv_sumterm(x, 0:min(length(x)-1, 10log10(length(x))), demean=demean)

# Cross-covariance sum term for range
function acv_sumterm(x::AbstractVector, y::AbstractVector, lags::Ranges; demean::Bool=true)
  lx, ly, llags, typex = length(x), length(y), length(lags), eltype(x)
  if lx != ly error("Input vectors must have same length") end
  if max(lags) > lx error("Cross-covariance distance must be less than sample size") end

  xs, ys = Array(typex, lx), Array(eltype(y), ly)
  if demean
    mx, my = mean(x), mean(y)
    for i = 1:lx; xs[i], ys[i] = x[i]-mx, y[i]-my; end
  else
    for i = 1:lx; xs[i], ys[i] = x[i], y[i]; end
  end

  autocov = Array(typex, llags)
  for i in 1:llags
    autocov[i] = dot(xs[1:end-lags[i]], ys[lags[i]+1:end])
  end
  autocov
end

# Cross-covariance sum term at a specific lag
acv_sumterm(x::AbstractVector, y::AbstractVector, lags::Real; demean::Bool=true) =
  acv_sumterm(x, y, lags:lags, demean=demean)[1]

# Cross-covariance sum term at a default of zero to 10log10(length(x)) lags
acv_sumterm(x::AbstractVector, y::AbstractVector; demean::Bool=true) =
  acv_sumterm(x, y, 0:min(length(x)-1, 10log10(length(x))), demean=demean)

# Cross-covariance sum term between all pairs of columns of a matrix for range
function acvall_sumterm(x::AbstractMatrix, lags::Ranges; demean::Bool=true)
  ncols = size(x, 2)

  autocov = Array(eltype(x), length(lags), ncols, ncols)
  for i = 1:ncols
    for j = 1:ncols
      autocov[:, i, j] = acv_sumterm(x[:, i], x[:, j], lags, demean=demean)
    end
  end
  autocov
end

# Cross-covariance sum term between all pairs of columns of a matrix at a specific lag
acvall_sumterm(x::AbstractMatrix, lags::Real; demean::Bool=true) =
  reshape(acvall_sumterm(x, lags:lags, demean=demean), size(x, 2), size(x, 2))

# Cross-covariance sum term between all pairs of columns of a matrix at a default of zero to 10log10(length(x)) lags
acvall_sumterm(x::AbstractMatrix; demean::Bool=true) =
  acvall_sumterm(x, 0:min(length(x)-1, 10log10(length(x))), demean=demean)

# Unlike acvall_sumterm, compute only autocovariance (not cross-covariance) sum term of matrix columns for range
function acvdiag_sumterm(x::AbstractMatrix, lags::Ranges; demean::Bool=true)
  ncols = size(x, 2)

  autocov = Array(eltype(x), length(lags), ncols)
  for i = 1:ncols
    autocov[:, i] = acv_sumterm(x[:, i], lags, demean=demean)
  end
  autocov
end

# Unlike acvall_sumterm, compute only autocovariance (not cross-covariance) sum term of matrix columns at a specific lag
acvdiag_sumterm(x::AbstractMatrix, lags::Real; demean::Bool=true) =
  acvdiag_sumterm(x, lags:lags, pop=pop, biased=biased)

# Unlike acvall_sumterm, compute only autocovariance (not cross-covariance) sum term of matrix columns at a default of
# zero to 10log10(length(x)) lags
acvdiag_sumterm(x::AbstractMatrix; demean::Bool=true) =
  acvdiag_sumterm(x, 0:min(length(x)-1, 10log10(length(x))), demean=demean)

# acv_sumterm wrapper for range
function acv_sumterm(x::AbstractMatrix, lags::Ranges; demean::Bool=true, diag::Bool=false)
 diag ? acvdiag_sumterm(x, lags, demean=demean) : acvall_sumterm(x, lags, demean=demean)
end

# acv_sumterm wrapper at a specific lag
function acv_sumterm(x::AbstractMatrix, lags::Real; demean::Bool=true, diag::Bool=false)
 diag ? acvdiag_sumterm(x, lags, demean=demean) : acvall_sumterm(x, lags, demean=demean)
end

# acv_sumterm wrapper at a default of zero to 10log10(length(x)) lags
function acv_sumterm(x::AbstractMatrix; demean::Bool=true, diag::Bool=false)
 diag ? acvdiag_sumterm(x, demean=demean) : acvall_sumterm(x, demean=demean)
end

# Autocovariance for range
acv(x::AbstractVector, lags::Ranges; demean::Bool=true) = acv_sumterm(x, lags, demean=demean)/length(x)

# Autocovariance at a specific lag
acv(x::AbstractVector, lags::Real; demean::Bool=true) = acv_sumterm(x, lags, demean=demean)/length(x)

# Autocovariance at a default of zero to 10log10(length(x)) lags
acv(x::AbstractVector; demean::Bool=true) = acv_sumterm(x, demean=demean)/length(x)

# Cross-covariance for range
acv(x::AbstractVector, y::AbstractVector, lags::Ranges; demean::Bool=true) =
  acv_sumterm(x, y, lags, demean=demean)/length(x)

# Cross-covariance at a specific lag
acv(x::AbstractVector, y::AbstractVector, lags::Real; demean::Bool=true) =
  acv_sumterm(x, y, lags, demean=demean)/length(x)

# Cross-covariance at a default of zero to 10log10(length(x)) lags
acv(x::AbstractVector, y::AbstractVector; demean::Bool=true) = acv_sumterm(x, y, demean=demean)/length(x)

# acv wrapper for range
acv(x::AbstractMatrix, lags::Ranges; demean::Bool=true, diag::Bool=false) =
  acv_sumterm(x, lags, demean=demean, diag=diag)/length(x)

# acv wrapper at a specific lag
acv(x::AbstractMatrix, lags::Real; demean::Bool=true, diag::Bool=false) =
  acv_sumterm(x, lags, demean=demean, diag=diag)/length(x)

# acv wrapper at a default of zero to 10log10(length(x)) lags
acv(x::AbstractMatrix; demean::Bool=true, diag::Bool=false) = acv_sumterm(x, demean=demean, diag=diag)/length(x)

# Autocorrelation for range
acf(x::AbstractVector, lags::Ranges) = acv_sumterm(x, lags)/((length(x)-1)*var(x))

# Autocorrelation at a specific lag
acf(x::AbstractVector, lags::Real) = acv_sumterm(x, lags)/((length(x)-1)*var(x))

# Autocorrelation at a default of zero to 10log10(length(x)) lags
acf(x::AbstractVector) = acv_sumterm(x)/((length(x)-1)*var(x))

# Cross-correlation for range
acf(x::AbstractVector, y::AbstractVector, lags::Ranges) = acv_sumterm(x, y, lags)/((length(x)-1)*var(x))

# Cross-correlation at a specific lag
acf(x::AbstractVector, y::AbstractVector, lags::Real) = acv_sumterm(x, y, lags)/((length(x)-1)*var(x))

# Cross-correlation at a default of zero to 10log10(length(x)) lags
acf(x::AbstractVector, y::AbstractVector) = acv_sumterm(x, y)/((length(x)-1)*var(x))

# acf wrapper for range
acf(x::AbstractMatrix, lags::Ranges, diag::Bool=false) = acv_sumterm(x, lags, diag=diag)/((length(x)-1)*var(x))

# acf wrapper at a specific lag
acf(x::AbstractMatrix, lags::Real, diag::Bool=false) = acv_sumterm(x, lags, diag=diag)/((length(x)-1)*var(x))

# acf wrapper at a default of zero to 10log10(length(x)) lags
acf(x::AbstractMatrix, diag::Bool=false) = acv_sumterm(x, diag=diag)/((length(x)-1)*var(x))
