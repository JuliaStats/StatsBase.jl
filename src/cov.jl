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
    vardim == 1 ? scale(s, x) :
    vardim == 2 ? scale(x, s) :
    error("vardim should be either 1 or 2.")
end

## scatter matrix

scattermat_zm(x::DenseMatrix, vardim::Int) = Base.unscaled_covzm(x, vardim)

scattermat_zm(x::DenseMatrix, wv::WeightVec, vardim::Int) =
    _symmetrize!(Base.unscaled_covzm(x, _scalevars(x, values(wv), vardim), vardim))

if VERSION < v"0.5.0-dev+679"
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
    """
        mean_and_cov(x::DenseMatrix, wv::WeightVec, vardim::Int)

    ### Args:
    * `x`: A 2-d array of type `DenseMatrix` to compute the mean and covariance of.
    * `wv`: An optional weight vector of type `WeightVec`.
    * `vardim`: An optional parameter specifying the dimension along which to compute the mean.

    Jointly compute the mean and covariance of input matrix `x`.
    """
    function mean_and_cov end

    """
        scattermat(x::DenseMatrix, mean, vardim::Int=1)

    ### Args:
    * `x`: A `DenseMatrix`
    * `mean`: Pre-computed mean vector. Default value of mean is set to nothing, which        indicates that the function would compute the mean internally. One can also set mean t    o 0, which indicates that the input `X` has already been centralized. Otherwise, the s    upplied mean will be subtracted from `X`.
    
    Compute scatter matrix for the variables contained in X.

    A scatter matrix can be considered as a unnormalized version of the covariance matrix.

    By default, it considers each column as a variable (i.e each row as an observation), a    nd subtract the mean from each vector. One may change this default behavior by setting    the keyword arguments.
    """
    function scattermat end
    
end
