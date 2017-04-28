## extended methods for computing covariance and scatter matrix

# auxiliary functions

function _symmetrize!(a::DenseMatrix)
    m, n = size(a)
    m == n || error("a must be a square matrix.")
    for j = 1:n
        @inbounds for i = j+1:n
            vl = a[i,j]
            vr = a[j,i]
            a[i,j] = a[j,i] = middle(vl, vr)
        end
    end
    return a
end

function _scalevars(x::DenseMatrix, s::DenseVector, vardim::Int)
    vardim == 1 ? Diagonal(s) * x :
    vardim == 2 ? x * Diagonal(s) :
    error("vardim should be either 1 or 2.")
end

## scatter matrix

scattermat_zm(x::DenseMatrix, vardim::Int) = Base.unscaled_covzm(x, vardim)

scattermat_zm(x::DenseMatrix, wv::AbstractWeights, vardim::Int) =
    _symmetrize!(Base.unscaled_covzm(x, _scalevars(x, values(wv), vardim), vardim))


"""
    scattermat(X, [wv::AbstractWeights]; mean=nothing, vardim=1)

Compute the scatter matrix, which is an unnormalized covariance matrix.
A weighting vector `wv` can be specified to weight
the estimate.

# Arguments
* `mean=nothing`: a known mean value. `nothing` indicates that the mean is
unknown, and the function will compute the mean. Specifying `mean=0` indicates
that the data are centered and hence there's no need to subtract the mean.
* `vardim=1`: the dimension along which the variables are organized.
When `vardim = 1`, the variables are considered columns with observations in rows;
when `vardim = 2`, variables are in rows with observations in columns.
"""
function scattermat end


"""
    cov(X, wv::AbstractWeights; mean=nothing, vardim=1)

Compute the weighted covariance matrix. By default, the covariance
matrix is normalized by the sum of the weights. That is, `cov(X, wv)`
is equivalent to `scattermat(X, wv) / sum(wv)`.
"""
cov


"""
    mean_and_cov(x, [wv::AbstractWeights]; vardim=1) -> (mean, cov)

Return the mean and covariance matrix as a tuple. A weighting
vector `wv` can be specified. `vardim` that designates whether
the variables are columns in the matrix (`1`) or rows (`2`).
"""
function mean_and_cov end

scattermatm(x::DenseMatrix, mean, vardim::Int=1) =
    scattermat_zm(x .- mean, vardim)

scattermatm(x::DenseMatrix, mean, wv::AbstractWeights, vardim::Int=1) =
    scattermat_zm(x .- mean, wv, vardim)

scattermat(x::DenseMatrix, vardim::Int=1) =
    scattermatm(x, Base.mean(x, vardim), vardim)

scattermat(x::DenseMatrix, wv::AbstractWeights, vardim::Int=1) =
    scattermatm(x, Base.mean(x, wv, vardim), wv, vardim)

## weighted cov
function Base.covm(x::DenseMatrix, mean, wv::AbstractWeights, vardim::Int=1, corrected::Bool=true)
    scale!(scattermatm(x, mean, wv, vardim), cfactor(wv, corrected))
end

function Base.cov(x::DenseMatrix, wv::AbstractWeights, vardim::Int=1; corrected=true)
    Base.covm(x, Base.mean(x, wv, vardim), wv, vardim, corrected)
end

function mean_and_cov(x::DenseMatrix, vardim::Int=1; corrected=true)
    m = mean(x, vardim)
    return m, Base.covm(x, m, vardim, corrected)
end

function mean_and_cov(x::DenseMatrix, wv::AbstractWeights, vardim::Int=1; corrected=true)
    m = mean(x, wv, vardim)
    return m, Base.cov(x, wv, vardim; corrected=corrected)
end
