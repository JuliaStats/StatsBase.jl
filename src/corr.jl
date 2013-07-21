# Various correlation 

#
# spearman covariance functions
#

# spearman covariance between two vectors
function cov_spearman(x::AbstractVector, y::AbstractVector)
    return cov(tiedrank(x), tiedrank(y))
end

# spearman covariance over all pairs of columns of two matrices
function cov_spearman(X::AbstractMatrix, Y::AbstractMatrix)
    return [cov_spearman(X[:,i], Y[:,j]) for i = 1:size(X, 2), j = 1:size(Y,2)]
end
function cov_spearman(x::AbstractVector, Y::AbstractMatrix)
    return [cov_spearman(x, Y[:,i]) for i = 1:size(Y, 2)]
end
function cov_spearman(X::AbstractMatrix, y::AbstractVector)
    return [cov_spearman(X[:,i], y) for i = 1:size(X, 2)]
end
# spearman covariance over all pairs of columns of a matrix
cov_spearman(X::AbstractMatrix) = cov(tiedrank(X, 1))

#
# spearman correlation functions
#

# spearman correlation between two vectors
function cor_spearman(x::AbstractVector, y::AbstractVector)
    return cor(tiedrank(x), tiedrank(y))
end

# spearman correlation over all pairs of columns of two matrices
function cor_spearman(X::AbstractMatrix, Y::AbstractMatrix)
    return cor(tiedrank(X, 1), tiedrank(Y, 1))
end
function cor_spearman(X::AbstractMatrix, y::AbstractVector)
    return cor(tiedrank(X, 1), tiedrank(y))
end
function cor_spearman(x::AbstractVector, Y::AbstractMatrix)
    return cor(tiedrank(x), tiedrank(Y, 1))
end

# spearman correlation over all pairs of columns of a matrix
function cor_spearman(X::AbstractMatrix)
    return cor(tiedrank(X, 1))
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
