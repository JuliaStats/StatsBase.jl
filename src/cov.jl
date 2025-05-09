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

function _scalevars(x::DenseMatrix, s::AbstractWeights, dims::Int)
    dims == 1 ? Diagonal(s) * x :
    dims == 2 ? x * Diagonal(s) :
    error("dims should be either 1 or 2.")
end

## scatter matrix

_unscaled_covzm(x::DenseMatrix, dims::Integer) = unscaled_covzm(x, dims)
_unscaled_covzm(x::DenseMatrix, wv::AbstractWeights, dims::Integer) =
    _symmetrize!(unscaled_covzm(x, _scalevars(x, wv, dims), dims))

"""
    scattermat(X, [wv::AbstractWeights]; mean=nothing, dims=1)

Compute the scatter matrix, which is an unnormalized covariance matrix.
A weighting vector `wv` can be specified to weight
the estimate.

# Arguments
* `mean=nothing`: a known mean value. `nothing` indicates that the mean is
  unknown, and the function will compute the mean. Specifying `mean=0` indicates
  that the data are centered and hence there's no need to subtract the mean.
* `dims=1`: the dimension along which the variables are organized.
  When `dims = 1`, the variables are considered columns with observations in rows;
  when `dims = 2`, variables are in rows with observations in columns.
"""
function scattermat end


"""
    cov(X, w::AbstractWeights, vardim=1; mean=nothing, corrected=false)

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
    mean_and_cov(x, [wv::AbstractWeights,] vardim=1; corrected=false) -> (mean, cov)

Return the mean and covariance matrix as a tuple. A weighting
vector `wv` can be specified. `vardim` that designates whether
the variables are columns in the matrix (`1`) or rows (`2`).
Finally, bias correction is applied to the covariance calculation if
`corrected=true`. See [`cov`](@ref) documentation for more details.
"""
function mean_and_cov end

scattermat(x::DenseMatrix; mean=nothing, dims::Int=1) =
    _scattermatm(x, mean, dims)
_scattermatm(x::DenseMatrix, ::Nothing, dims::Int) =
    _unscaled_covzm(x .- mean(x, dims=dims), dims)
_scattermatm(x::DenseMatrix, mean, dims::Int=1) =
    _unscaled_covzm(x .- mean, dims)

scattermat(x::DenseMatrix, wv::AbstractWeights; mean=nothing, dims::Int=1) =
    _scattermatm(x, wv, mean, dims)
_scattermatm(x::DenseMatrix, wv::AbstractWeights, ::Nothing, dims::Int) =
    _unscaled_covzm(x .- mean(x, wv, dims=dims), wv, dims)
_scattermatm(x::DenseMatrix, wv::AbstractWeights, mean, dims::Int) =
    _unscaled_covzm(x .- mean, wv, dims)

## weighted cov
covm(x::DenseMatrix, mean, w::AbstractWeights, dims::Int=1;
     corrected::Union{Bool, Nothing}=nothing) =
    rmul!(scattermat(x, w, mean=mean, dims=dims), varcorrection(w, depcheck(:covm, :corrected, corrected)))


cov(x::DenseMatrix, w::AbstractWeights, dims::Int=1; corrected::Union{Bool, Nothing}=nothing) =
    covm(x, mean(x, w, dims=dims), w, dims; corrected=depcheck(:cov, :corrected, corrected))

function corm(x::DenseMatrix, mean, w::AbstractWeights, vardim::Int=1)
    c = covm(x, mean, w, vardim; corrected=false)
    s = std(x, w, vardim; mean=mean, corrected=false)
    cov2cor!(c, s)
end

"""
    cor(X, w::AbstractWeights, dims=1)

Compute the Pearson correlation matrix of `X` along the dimension
`dims` with a weighting `w` .
"""
cor(x::DenseMatrix, w::AbstractWeights, dims::Int=1) =
    corm(x, mean(x, w, dims=dims), w, dims)

function mean_and_cov(x::DenseMatrix, dims::Int=1; corrected::Bool=true)
    m = mean(x, dims=dims)
    return m, covm(x, m, dims, corrected=corrected)
end
function mean_and_cov(x::DenseMatrix, wv::AbstractWeights, dims::Int=1;
                      corrected::Union{Bool, Nothing}=nothing)
    m = mean(x, wv, dims=dims)
    return m, cov(x, wv, dims; corrected=depcheck(:mean_and_cov, :corrected, corrected))
end


"""
    cov2cor(C::AbstractMatrix, [s::AbstractArray])

Compute the correlation matrix from the covariance matrix `C` and, optionally, a vector
of standard deviations `s`. Use [`StatsBase.cov2cor!`](@ref) for an in-place version.
"""
function cov2cor(C::AbstractMatrix, s::AbstractArray = map(sqrt, view(C, diagind(C))))
    zs = zero(eltype(s))
    T = typeof(zero(eltype(C)) / (zs * zs))
    return cov2cor!(copyto!(similar(C, T), C), s)
end

# Original implementation: https://github.com/JuliaStats/Statistics.jl/blob/22dee82f9824d6045e87aa4b97e1d64fe6f01d8d/src/Statistics.jl#L633-L657
"""
    cov2cor!(C::AbstractMatrix, [s::AbstractArray])

Convert the covariance matrix `C` to a correlation matrix in-place, optionally using a vector of
standard deviations `s`.
"""
function cov2cor!(C::AbstractMatrix, s::AbstractArray = map(sqrt, view(C, diagind(C))))
    Base.require_one_based_indexing(C, s)
    n = length(s)
    size(C) == (n, n) || throw(DimensionMismatch("inconsistent dimensions"))
    for j = 1:n
        sj = s[j]
        for i = 1:(j-1)
            C[i,j] = adjoint(C[j,i])
        end
        C[j,j] = oneunit(C[j,j])
        for i = (j+1):n
            C[i,j] = _clampcor(C[i,j] / (s[i] * sj))
        end
    end
    return C
end
_clampcor(x::Real) = clamp(x, -1, 1)
_clampcor(x) = x

# Preserve structure of Symmetric and Hermitian covariance matrices
function cov2cor!(C::Union{Symmetric{<:Real},Hermitian}, s::AbstractArray)
    n = length(s)
    size(C) == (n, n) || throw(DimensionMismatch("inconsistent dimensions"))
    A = parent(C)
    if C.uplo === 'U'
        for j = 1:n
            sj = s[j]
            for i = 1:(j-1)
                A[i,j] = _clampcor(A[i,j] / (s[i] * sj))
            end
            A[j,j] = oneunit(A[j,j])
        end
    else
        for j = 1:n
            sj = s[j]
            A[j,j] = oneunit(A[j,j])
            for i = (j+1):n
                A[i,j] = _clampcor(A[i,j] / (s[i] * sj))
            end
        end
    end
    return C
end

"""
    cor2cov(C, s)

Compute the covariance matrix from the correlation matrix `C` and a vector of standard
deviations `s`. Use [`StatsBase.cor2cov!`](@ref) for an in-place version.
"""
function cor2cov(C::AbstractMatrix, s::AbstractArray)
    zs = zero(eltype(s))
    T = typeof(zero(eltype(C)) * (zs * zs))
    return cor2cov!(copyto!(similar(C, T), C), s)
end

"""
    cor2cov!(C, s)

Convert the correlation matrix `C` to a covariance matrix in-place using a vector of
standard deviations `s`.
"""
function cor2cov!(C::AbstractMatrix, s::AbstractArray)
    Base.require_one_based_indexing(C, s)
    n = length(s)
    size(C) == (n, n) || throw(DimensionMismatch("inconsistent dimensions"))
    for j in 1:n
        sj = s[j]
        for i in 1:(j-1)
            C[i,j] = adjoint(C[j,i])
        end
        C[j,j] = sj^2
        for i in (j+1):n
            C[i,j] *= s[i] * sj
        end
    end
    return C
end

# Preserve structure of Symmetric and Hermitian correlation matrices
function cor2cov!(C::Union{Symmetric{<:Real},Hermitian}, s::AbstractArray)
    n = length(s)
    size(C) == (n, n) || throw(DimensionMismatch("inconsistent dimensions"))
    A = parent(C)
    if C.uplo === 'U'
        for j in 1:n
            sj = s[j]
            for i in 1:(j-1)
                A[i,j] *= s[i] * sj
            end
            A[j,j] = sj^2
        end
    else
        for j in 1:n
            sj = s[j]
            A[j,j] = sj^2
            for i in (j+1):n
                A[i,j] *= s[i] * sj
            end
        end
    end
    return C
end

"""
    CovarianceEstimator

Abstract type for covariance estimators.
"""
abstract type CovarianceEstimator end

"""
    cov(ce::CovarianceEstimator, x::AbstractVector; mean=nothing)

Compute a variance estimate from the observation vector `x` using the  estimator `ce`.
"""
cov(ce::CovarianceEstimator, x::AbstractVector; mean=nothing) =
    error("cov is not defined for $(typeof(ce)) and $(typeof(x))")

"""
    cov(ce::CovarianceEstimator, x::AbstractVector, y::AbstractVector)

Compute the covariance of the vectors `x` and `y` using estimator `ce`.
"""
cov(ce::CovarianceEstimator, x::AbstractVector, y::AbstractVector) =
    error("cov is not defined for $(typeof(ce)), $(typeof(x)) and $(typeof(y))")

"""
    cov(ce::CovarianceEstimator, X::AbstractMatrix, [w::AbstractWeights];
        mean=nothing, dims::Int=1)

Compute the covariance matrix of the matrix `X` along dimension `dims`
using estimator `ce`. A weighting vector `w` can be specified.
The keyword argument `mean` can be:

* `nothing` (default) in which case the mean is estimated and subtracted
  from the data `X`,
* a precalculated mean in which case it is subtracted from the data `X`.
  Assuming `size(X)` is `(N,M)`, `mean` can either be:
  * when `dims=1`, an `AbstractMatrix` of size `(1,M)`,
  * when `dims=2`, an `AbstractVector` of length `N` or an `AbstractMatrix`
    of size `(N,1)`.
"""
cov(ce::CovarianceEstimator, X::AbstractMatrix; mean=nothing, dims::Int=1) =
    error("cov is not defined for $(typeof(ce)) and $(typeof(X))")

cov(ce::CovarianceEstimator, X::AbstractMatrix, w::AbstractWeights; mean=nothing, dims::Int=1) =
    error("cov is not defined for $(typeof(ce)), $(typeof(X)) and $(typeof(w))")

"""
    var(ce::CovarianceEstimator, x::AbstractVector; mean=nothing)

Compute the variance of the vector `x` using the estimator `ce`.
"""
var(ce::CovarianceEstimator, x::AbstractVector; kwargs...) = cov(ce, x; kwargs...)

"""
    std(ce::CovarianceEstimator, x::AbstractVector; mean=nothing)

Compute the standard deviation of the vector `x` using the estimator `ce`.
"""
std(ce::CovarianceEstimator, x::AbstractVector; kwargs...) = sqrt(var(ce, x; kwargs...))

"""
    cor(ce::CovarianceEstimator, x::AbstractVector, y::AbstractVector)

Compute the correlation of the vectors `x` and `y` using estimator `ce`.
"""
function cor(ce::CovarianceEstimator, x::AbstractVector, y::AbstractVector)
    # Here we allow `ce` to see both `x` and `y` simultaneously, and allow it to compute
    # a full covariance matrix, from which we will extract the correlation.
    #
    # Whilst in some cases it might be equivalent (and possibly more efficient) to use:
    #   cov(ce, x, y) / (std(ce, x) * std(ce, y)),
    # this need not apply in general.
    return cor(ce, hcat(x, y))[1, 2]
end

"""
    cor(ce::CovarianceEstimator, X::AbstractMatrix, [w::AbstractWeights];
        mean=nothing, dims::Int=1)

Compute the correlation matrix of the matrix `X` along dimension `dims`
using estimator `ce`. A weighting vector `w` can be specified.
The keyword argument `mean` can be:

* `nothing` (default) in which case the mean is estimated and subtracted
  from the data `X`,
* a precalculated mean in which case it is subtracted from the data `X`.
  Assuming `size(X)` is `(N,M)`, `mean` can either be:
  * when `dims=1`, an `AbstractMatrix` of size `(1,M)`,
  * when `dims=2`, an `AbstractVector` of length `N` or an `AbstractMatrix`
    of size `(N,1)`.
"""
function cor(ce::CovarianceEstimator, X::AbstractMatrix; kwargs...)
    return cov2cor!(cov(ce, X; kwargs...))
end
function cor(ce::CovarianceEstimator, X::AbstractMatrix, w::AbstractWeights; kwargs...)
    return cov2cor!(cov(ce, X, w; kwargs...))
end

"""
    SimpleCovariance(;corrected::Bool=false)

Simple covariance estimator. Estimation calls `cov(x; corrected=corrected)`,
`cov(x, y; corrected=corrected)` or `cov(X, w, dims; corrected=corrected)`
where `x`, `y` are vectors, `X` is a matrix and `w` is a weighting vector.
"""
struct SimpleCovariance <: CovarianceEstimator
    corrected::Bool
    SimpleCovariance(;corrected::Bool=false) = new(corrected)
end

cov(sc::SimpleCovariance, x::AbstractVector) =
    cov(x; corrected=sc.corrected)

cov(sc::SimpleCovariance, x::AbstractVector, y::AbstractVector) =
    cov(x, y; corrected=sc.corrected)

function cov(sc::SimpleCovariance, X::AbstractMatrix; dims::Int=1, mean=nothing)
    dims ∈ (1, 2) || throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    if mean === nothing
        return cov(X; dims=dims, corrected=sc.corrected)
    else
        return covm(X, mean, dims, corrected=sc.corrected)
    end
end

function cov(sc::SimpleCovariance, X::AbstractMatrix, w::AbstractWeights; dims::Int=1, mean=nothing)
    dims ∈ (1, 2) || throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    if mean === nothing
        return cov(X, w, dims, corrected=sc.corrected)
    else
        return covm(X, mean, w, dims, corrected=sc.corrected)
    end
end
