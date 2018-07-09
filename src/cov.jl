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


scattermat_zm(x::DenseMatrix, vardim::Int) = unscaled_covzm(x, vardim)


scattermat_zm(x::DenseMatrix, wv::AbstractWeights, vardim::Int) =
    _symmetrize!(unscaled_covzm(x, _scalevars(x, values(wv), vardim), vardim))

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
    cov(X, w::AbstractWeights; mean=nothing, vardim=1, corrected=false)

Compute the weighted covariance matrix. Similar to `var` and `std` the biased covariance
matrix (`corrected=false`) is computed by multiplying `scattermat(X, w)` by
``\\frac{1}{\\sum{w}}`` to normalize. However, the unbiased covariance matrix
(`corrected=true`) is dependent on the type of weights used:
* `AnalyticWeights`: ``\\frac{1}{\\sum w - \\sum {w^2} / \\sum w}``
* `FrequencyWeights`: ``\\frac{1}{\\sum{w} - 1}``
* `ProbabilityWeights`: ``\\frac{n}{(n - 1) \\sum w}`` where ``n`` equals `count(!iszero, w)`
* `Weights`: `ArgumentError` (bias correction not supported)
"""
cov


"""
    mean_and_cov(x, [wv::AbstractWeights]; vardim=1, corrected=false) -> (mean, cov)

Return the mean and covariance matrix as a tuple. A weighting
vector `wv` can be specified. `vardim` that designates whether
the variables are columns in the matrix (`1`) or rows (`2`).
Finally, bias correction is applied to the covariance calculation if
`corrected=true`. See [`cov`](@ref) documentation for more details.
"""
function mean_and_cov end


scattermatm(x::DenseMatrix, mean, vardim::Int=1) =
    scattermat_zm(x .- mean, vardim)

scattermatm(x::DenseMatrix, mean, wv::AbstractWeights, vardim::Int=1) =
    scattermat_zm(x .- mean, wv, vardim)

scattermat(x::DenseMatrix, vardim::Int=1) =
    scattermatm(x, mean(x, dims = vardim), vardim)

scattermat(x::DenseMatrix, wv::AbstractWeights, vardim::Int=1) =
    scattermatm(x, Statistics.mean(x, wv, vardim), wv, vardim)

## weighted cov
covm(x::DenseMatrix, mean, w::AbstractWeights, vardim::Int=1;
     corrected::DepBool=nothing) =
    rmul!(scattermatm(x, mean, w, vardim), varcorrection(w, depcheck(:covm, corrected)))


cov(x::DenseMatrix, w::AbstractWeights, vardim::Int=1; corrected::DepBool=nothing) =
    covm(x, mean(x, w, vardim), w, vardim; corrected=depcheck(:cov, corrected))

function corm(x::DenseMatrix, mean, w::AbstractWeights, vardim::Int=1)
    c = covm(x, mean, w, vardim; corrected=false)
    s = stdm(x, w, mean, vardim; corrected=false)
    cov2cor!(c, s)
end

"""
    cor(X, w::AbstractWeights, vardim=1)

Compute the Pearson correlation matrix of `X` along the dimension
`vardim` with a weighting `w` .
"""
cor(x::DenseMatrix, w::AbstractWeights, vardim::Int=1) =
    corm(x, mean(x, w, vardim), w, vardim)

function mean_and_cov(x::DenseMatrix, vardim::Int=1; corrected::Bool=true)
    m = mean(x, dims = vardim)
    return m, covm(x, m, vardim, corrected=corrected)
end
function mean_and_cov(x::DenseMatrix, wv::AbstractWeights, vardim::Int=1;
                      corrected::DepBool=nothing)
    m = mean(x, wv, vardim)
    return m, cov(x, wv, vardim; corrected=depcheck(:mean_and_cov, corrected))
end

"""
    cov2cor(C, s)

Compute the correlation matrix from the covariance matrix `C` and a vector of standard
deviations `s`. Use `StatsBase.cov2cor!` for an in-place version.
"""
cov2cor(C::AbstractMatrix, s::AbstractArray) = cov2cor!(copy(C), s)

"""
    cor2cov(C, s)

Compute the covariance matrix from the correlation matrix `C` and a vector of standard
deviations `s`. Use `StatsBase.cor2cov!` for an in-place version.
"""
cor2cov(C::AbstractMatrix, s::AbstractArray) = cor2cov!(copy(C), s)

"""
    cor2cov!(C, s)

Converts the correlation matrix `C` to a covariance matrix in-place using a vector of
standard deviations `s`.
"""
function cor2cov!(C::AbstractMatrix, s::AbstractArray)
    n = length(s)
    size(C) == (n, n) || throw(DimensionMismatch("inconsistent dimensions"))
    for i in CartesianIndices(size(C))
        @inbounds C[i] *= s[i[1]] * s[i[2]]
    end
    return C
end
