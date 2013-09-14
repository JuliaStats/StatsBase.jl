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

# autocovariance for range
function acv(x::AbstractVector, lags::Ranges; pop::Bool=true, biased::Bool=false)
  lx, llags = length(x), length(lags)
  if max(lags) > lx; error("Autocovariance distance must be less than sample size"); end
  if (pop == true) && (biased == true); error("Only unbiased estimator of population acv is available"); end

  mx = mean(x)
  xs = Array(typeof(mx), lx)
  for i = 1:lx
    xs[i] = x[i] - mx
  end

  autocov = Array(typeof(mx), llags)
  for i in 1:llags
    autocov[i] = dot(xs[1:end - lags[i]], xs[lags[i] + 1:end])
  end
  
  pop ? (return autocov/lx) : (biased ? (return autocov/(lx-1)) : (return autocov/llags))
end

# autocovariance at a specific lag
acv(x::AbstractVector, lags::Real; pop::Bool=true, biased::Bool=false) = acv(x, lags:lags, pop=pop, biased=biased)[1]

# autocovariance at a default of zero to 10log10(length(x)) lags
acv(x::AbstractVector; pop::Bool=true, biased::Bool=false) = 
  acv(x, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# crosscovariance between two vectors for range
function acv(x::AbstractVector, y::AbstractVector, lags::Ranges; pop::Bool=true, biased::Bool=false)
  lx, ly, llags = length(x), length(y), length(lags)
  if lx != ly error("Input vectors must have same length") end
  if max(lags) > lx error("Autocorrelation distance must be less than sample size") end
  if (pop == true) && (biased == true); error("Only unbiased estimator of population acv is available"); end

  mx, my = mean(x), mean(y)
  xs, ys = Array(typeof(mx), lx), Array(typeof(my), ly)
  for i = 1:lx
    xs[i], ys[i] = x[i] - mx, y[i] - my
  end

  autocov = Array(typeof(mx), llags)
  for i in 1:llags
    autocov[i] = dot(xs[1:end - lags[i]], ys[lags[i] + 1:end])
  end

  pop ? (return autocov/lx) : (biased ? (return autocov/(lx-1)) : (return autocov/llags))
end

# crosscovariance between two vectors at a specific lag
acv(x::AbstractVector, y::AbstractVector, lags::Real; pop::Bool=true, biased::Bool=false) =
  acv(x, y, lags:lags, pop=pop, biased=biased)[1]

# crosscovariance between two vectors at a default of zero to 10log10(length(x)) lags
acv(x::AbstractVector, y::AbstractVector; pop::Bool=true, biased::Bool=false) =
  acv(x, y, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# crosscovariance between all pairs of columns of a matrix for range
function acvall(x::AbstractMatrix, lags::Ranges; pop::Bool=true, biased::Bool=false)
  ncols = size(x, 2)
  autocov = Array(eltype(x), length(lags), ncols, ncols)

  for i = 1:ncols
    for j = 1:ncols
      autocov[:, i, j] = acv(x[:, i], x[:, j], lags, pop=pop, biased=biased)
    end
  end

  return autocov
end

# crosscovariance between all pairs of columns of a matrix at a specific lag
acvall(x::AbstractMatrix, lags::Real; pop::Bool=true, biased::Bool=false) =
  reshape(acvall(x, lags:lags, pop=pop, biased=biased), size(x, 2), size(x, 2))

# crosscovariance between all pairs of columns of a matrix at a default of zero to 10log10(length(x)) lags
acvall(x::AbstractMatrix; pop::Bool=true, biased::Bool=false) =
  acvall(x, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# Unlike acvall, compute only autocovariance (not crosscovariance) of matrix columns for range
function acvdiag(x::AbstractMatrix, lags::Ranges; pop::Bool=true, biased::Bool=false)
  ncols = size(x, 2)
  autocov = Array(eltype(x), length(lags), ncols)

  for i = 1:ncols
    autocov[:, i] = acv(x[:, i], lags, pop=pop, biased=biased)
  end

  return autocov
end

# Unlike acvall, compute only autocovariance (not crosscovariance) of matrix columns at a specific lag
acvdiag(x::AbstractMatrix, lags::Real; pop::Bool=true, biased::Bool=false) =
  acvdiag(x, lags:lags, pop=pop, biased=biased)

# Unlike acvall, compute only autocovariance (not crosscovariance) of matrix columns at a default of zero to 10log10(length(x)) lags
acvdiag(x::AbstractMatrix; pop::Bool=true, biased::Bool=false) =
  acvdiag(x, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# acv wrapper for range
function acv(x::AbstractMatrix, lags::Ranges; pop::Bool=true, biased::Bool=false, diag::Bool=false)
 diag ? acvdiag(x, lags, pop=pop, biased=biased) : acvall(x, lags, pop=pop, biased=biased)
end

# acv wrapper at a specific lag
function acv(x::AbstractMatrix, lags::Real; pop::Bool=true, biased::Bool=false, diag::Bool=false)
 diag ? acvdiag(x, lags, pop=pop, biased=biased) : acvall(x, lags, pop=pop, biased=biased)
end

# acv wrapper at a default of zero to 10log10(length(x)) lags
function acv(x::AbstractMatrix; pop::Bool=true, biased::Bool=false, diag::Bool=false)
 diag ? acvdiag(x, pop=pop, biased=biased) : acvall(x, pop=pop, biased=biased)
end

# autocorrelation for range
acf(x::AbstractVector, lags::Ranges; pop::Bool=false, biased::Bool=true) =
  acv(x, lags, pop=pop, biased=biased)/var(x)

# autocorrelation at a specific lag
acf(x::AbstractVector, lags::Real; pop::Bool=false, biased::Bool=true) = acf(x, lags:lags, pop=pop, biased=biased)[1]

# autocorrelation at a default of zero to 10log10(length(x)) lags
acf(x::AbstractVector; pop::Bool=false, biased::Bool=true) =
  acf(x, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# crosscorrelation for range
acf(x::AbstractVector, y::AbstractVector, lags::Ranges; pop::Bool=false, biased::Bool=true) =
  acv(x, y, lags, pop=pop, biased=biased)/(std(x)*std(y))

# crosscorrelation at a specific lag
acf(x::AbstractVector, y::AbstractVector, lags::Real; pop::Bool=false, biased::Bool=true) =
  acf(x, y, lags:lags, pop=pop, biased=biased)[1]

# crosscorrelation at a default of zero to 10log10(length(x)) lags
acf(x::AbstractVector, y::AbstractVector; pop::Bool=false, biased::Bool=true) =
  acf(x, y, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# crosscorrelation between all pairs of columns of a matrix for range
function acfall(x::AbstractMatrix, lags::Ranges; pop::Bool=false, biased::Bool=true)
  ncols = size(x, 2)
  autocorr = Array(eltype(x), length(lags), ncols, ncols)

  for i = 1:ncols
    for j = 1:ncols
      autocorr[:, i, j] = acf(x[:, i], x[:, j], lags, pop=pop, biased=biased)
    end
  end

  return autocorr
end

# crosscorrelation between all pairs of columns of a matrix at a specific lag
acfall(x::AbstractMatrix, lags::Real; pop::Bool=false, biased::Bool=true) =
  reshape(acfall(x, lags:lags, pop=pop, biased=biased), size(x, 2), size(x, 2))

# crosscorrelation between all pairs of columns of a matrix at a default of zero to 10log10(length(x)) lags
acfall(x::AbstractMatrix; pop::Bool=false, biased::Bool=true) =
  acfall(x, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# Unlike acfall, compute only autocorrelation (not crosscorrelation) of matrix columns for range
function acfdiag(x::AbstractMatrix, lags::Ranges; pop::Bool=false, biased::Bool=true)
  ncols = size(x, 2)
  autocorr = Array(eltype(x), length(lags), ncols)

  for i = 1:ncols
    autocorr[:, i] = acf(x[:, i], lags, pop=pop, biased=biased)
  end

  return autocorr
end

# Unlike acfall, compute only autocorrelation (not crosscorrelation) of matrix columns at a specific lag
acfdiag(x::AbstractMatrix, lags::Real; pop::Bool=false, biased::Bool=true) =
  acfdiag(x, lags:lags, pop=pop, biased=biased)

# Unlike acfall, compute only autocorrelation (not crosscorrelation) of matrix columns at a default of zero to 10log10(length(x)) lags
acfdiag(x::AbstractMatrix; pop::Bool=false, biased::Bool=true) =
  acfdiag(x, 0:min(length(x)-1, 10log10(length(x))), pop=pop, biased=biased)

# acf wrapper for range
function acf(x::AbstractMatrix, lags::Ranges; pop::Bool=false, biased::Bool=true, diag::Bool=false)
 diag ? acfdiag(x, lags, pop=pop, biased=biased) : acfall(x, lags, pop=pop, biased=biased)
end

# acf wrapper at a specific lag
function acf(x::AbstractMatrix, lags::Real; pop::Bool=false, biased::Bool=true, diag::Bool=false)
 diag ? acfdiag(x, lags, pop=pop, biased=biased) : acfall(x, lags, pop=pop, biased=biased)
end

# acf wrapper at a default of zero to 10log10(length(x)) lags
function acf(x::AbstractMatrix; pop::Bool=false, biased::Bool=true, diag::Bool=false)
 diag ? acfdiag(x, pop=pop, biased=biased) : acfall(x, pop=pop, biased=biased)
end
