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
scattermatm(x::DenseMatrix, mean, vardim::Int=1) =
    scattermat_zm(x .- mean, vardim)

scattermatm(x::DenseMatrix, mean, wv::AbstractWeights, vardim::Int=1) =
    scattermat_zm(x .- mean, wv, vardim)

scattermat(x::DenseMatrix, vardim::Int=1) =
    scattermatm(x, Base.mean(x, vardim), vardim)

scattermat(x::DenseMatrix, wv::AbstractWeights, vardim::Int=1) =
    scattermatm(x, Base.mean(x, wv, vardim), wv, vardim)

scattermat_zm(x::DenseMatrix, vardim::Int) = Base.unscaled_covzm(x, vardim)

scattermat_zm(x::DenseMatrix, wv::AbstractWeights, vardim::Int) =
    _symmetrize!(Base.unscaled_covzm(x, _scalevars(x, values(wv), vardim), vardim))

"""
    cov(X, wv::AbstractWeights, [vardim, corrected])

Compute the weighted covariance matrix. Similar to `var` and `std` the biased covariance
matrix (`corrected=false`) can be computed by multiplying `scattermat(X, wv)` by
``\frac{1}{\sum{w}}`` to normalize. However, the unbiased covariance matrix
(`corrected=true`) is dependent on the type of weights used:

* AnalyticWeights: ``\\frac{1}{\sum w - \sum {w^2} / \sum{w}^2}``
* FrequencyWeights: ``\\frac{1}{\sum{w} - 1}``
* ProbabilityWeights: ``\\frac{n}{(n - 1) \sum w}`` where `n = length(w)`
"""
Base.cov(x::DenseMatrix, wv::AbstractWeights, corrected::Bool) =
    Base.covm(x, Base.mean(x, wv, 1), wv, 1, corrected)

Base.cov(x::DenseMatrix, wv::AbstractWeights, vardim::Int, corrected::Bool) =
    Base.covm(x, Base.mean(x, wv, vardim), wv, vardim, corrected)

Base.covm(x::DenseMatrix, mean, wv::AbstractWeights, corrected::Bool) =
    scale!(scattermatm(x, mean, wv, 1), varcorrection(wv, corrected))

Base.covm(x::DenseMatrix, mean, wv::AbstractWeights, vardim::Int, corrected::Bool) =
    scale!(scattermatm(x, mean, wv, vardim), varcorrection(wv, corrected))

"""
    mean_and_cov(x, [wv::AbstractWeights, vardim, corrected]) -> (mean, cov)

Return the mean and covariance matrix as a tuple. A weighting
vector `wv` can be specified. `vardim` that designates whether
the variables are columns in the matrix (`1`) or rows (`2`).
Finally, bias correction can be applied to the covariance calculation if
`corrected=true`.
See `cov` documentation for more details.
"""
function mean_and_cov(x::DenseMatrix, corrected::Bool=true)
    m = mean(x, 1)
    return m, Base.covm(x, m, 1, corrected)
end

function mean_and_cov(x::DenseMatrix, vardim::Int, corrected::Bool=true)
    m = mean(x, vardim)
    return m, Base.covm(x, m, vardim, corrected)
end

function mean_and_cov(x::DenseMatrix, wv::AbstractWeights, corrected::Bool)
    m = mean(x, wv, 1)
    return m, Base.cov(x, wv, 1, corrected)
end

function mean_and_cov(x::DenseMatrix, wv::AbstractWeights, vardim::Int, corrected::Bool)
    m = mean(x, wv, vardim)
    return m, Base.cov(x, wv, vardim, corrected)
end
