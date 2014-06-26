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

mean_and_cov(x::DenseMatrix; vardim::Int=1) = (m = mean(x, vardim); (m, Base.covm(x, m; vardim=vardim)))
mean_and_cov(x::DenseMatrix, wv::WeightVec; vardim::Int=1) = 
    (m = mean(x, wv, vardim); (m, Base.cov(x, wv; mean=m, vardim=vardim)))
