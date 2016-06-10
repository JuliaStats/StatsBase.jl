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

scattermat_zm(x::DenseMatrix, wv::WeightVec, vardim::Int) =
    _symmetrize!(Base.unscaled_covzm(x, _scalevars(x, values(wv), vardim), vardim))

if VERSION < v"0.5.0-dev+679"
"""
    scattermat(x, [wv,] [mean=nothing,] [vardim=1]) -> Matrix

Compute the scatter matrix, optionally specifying a weighting vector.
The function accepts two optional keyword arguments:

  • `mean`: a known mean value. By default it is set to `nothing`,
    which indicates that the mean is unknown, and the function will
    compute the mean. Specifying `mean=0` indicates that the data
    are centered and hence there's no need to subtract the mean.
  • `vardim`: the dimension of the variables. When `vardim = 1`, the
    variables are considered columns with observations in rows; when
    `vardim = 2`, variables are in rows with observations in columns.
    By default it's set to 1.
"""
    function scattermat(x::DenseMatrix; mean=nothing, vardim::Int=1)
        mean == 0 ? scattermat_zm(x, vardim) :
        mean == nothing ? scattermat_zm(x .- Base.mean(x, vardim), vardim) :
        scattermat_zm(x .- mean, vardim)
    end

    function scattermat(x::DenseMatrix, wv::WeightVec; mean=nothing, vardim::Int=1)
        mean == 0 ? scattermat_zm(x, wv, vardim) :
        mean == nothing ? scattermat_zm(x .- Base.mean(x, wv, vardim), wv, vardim) :
        scattermat_zm(x .- mean, wv, vardim)
    end

    ## weighted cov
    Base.cov(x::DenseMatrix, wv::WeightVec; mean=nothing, vardim::Int=1) =
        scale!(scattermat(x, wv; mean=mean, vardim=vardim), inv(sum(wv)))


"""
    mean_and_cov(x, [wv,] [vardim=1]) -> (Vector, Matrix)

Return the means and covariance matrix as a tuple. A weighting
vector `wv` can be specified. The function accepts a single
keyword argument `vardim` that designates whether the variables
are columns in the matrix (`1`, the default) or rows (`2`).
"""
    function mean_and_cov(x::DenseMatrix; vardim::Int=1)
        m = mean(x, vardim)
        return m, Base.covm(x, m; vardim=vardim)
    end
    function mean_and_cov(x::DenseMatrix, wv::WeightVec; vardim::Int=1)
        m = mean(x, wv, vardim)
        return m, Base.cov(x, wv; mean=m, vardim=vardim)
    end
else
    scattermatm(x::DenseMatrix, mean, vardim::Int=1) =
        scattermat_zm(x .- mean, vardim)

    scattermatm(x::DenseMatrix, mean, wv::WeightVec, vardim::Int=1) =
        scattermat_zm(x .- mean, wv, vardim)

    scattermat(x::DenseMatrix, vardim::Int=1) =
        scattermatm(x, Base.mean(x, vardim), vardim)

    scattermat(x::DenseMatrix, wv::WeightVec, vardim::Int=1) =
        scattermatm(x, Base.mean(x, wv, vardim), wv, vardim)

    ## weighted cov
    Base.covm(x::DenseMatrix, mean, wv::WeightVec, vardim::Int=1) =
        scale!(scattermatm(x, mean, wv, vardim), inv(sum(wv)))

    Base.cov(x::DenseMatrix, wv::WeightVec, vardim::Int=1) =
        Base.covm(x, Base.mean(x, wv, vardim), wv, vardim)

    function mean_and_cov(x::DenseMatrix, vardim::Int=1)
        m = mean(x, vardim)
        return m, Base.covm(x, m, vardim)
    end
    function mean_and_cov(x::DenseMatrix, wv::WeightVec, vardim::Int=1)
        m = mean(x, wv, vardim)
        return m, Base.cov(x, wv, vardim)
    end
end
