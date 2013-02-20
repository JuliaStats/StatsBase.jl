module Stats
    import Base.mean, Base.quantile

    export autocor,
           cor_spearman,
           cov_spearman,
           decile,
           distances,
           distances,
           ecdf,
           inverse_rle,
           iqr,
           kurtosis,
           mad,
           percentile,
           quantile,
           quartile,
           quintile,
           rle,
           skewness,
           tiedrank,
           weighted_mean

    mean(v::AbstractArray, dim::Int) = squeeze(sum(v, dim), dim) / size(v, dim)

    # tensor mean along a vector of dimensions
    # For example, if A is an array with dimension M x N x P, then mean(A, [2 3])
    # returns a length-M vector after taking the mean along the 2nd and 3rd dimension.
    function mean{N}(v::AbstractArray, dims::Array{Int,N})
       r = reducedim(+, v, dims, 0)
       s = size(v)
       r /= prod(s[dims])
       squeeze(r, dims)
   end

    weighted_mean(v::AbstractArray, w::AbstractArray) = sum(v .* w) / sum(w)

    # median absolute deviation with known center with consistency adjustment
    mad(v::AbstractArray, center::Number) = 1.4826 * median(abs(v - center))

    # median absolute deviation
    mad(v::AbstractArray) = mad(v, median(v))

    # maximum likelihood estimate of skewness with known mean m
    function skewness(v::AbstractVector, m::Number)
        n = length(v)
        empirical_third_centered_moment = 0.0
        empirical_variance = 0.0
        for x_i in v
            empirical_third_centered_moment += (x_i - m)^3
            empirical_variance += (x_i - m)^2
        end
        empirical_third_centered_moment /= n
        empirical_variance /= n
        return empirical_third_centered_moment / (empirical_variance^1.5)
    end

    # maximum likelihood estimate of skewness
    skewness(v::AbstractVector) = skewness(v, mean(v))

    # maximum likelihood estimate of kurtosis with known mean m
    function kurtosis(v::AbstractVector, m::Number)
        n = length(v)
        empirical_fourth_centered_moment = 0.0
        empirical_variance = 0.0
        for x_i in v
            empirical_fourth_centered_moment += (x_i - m)^4
            empirical_variance += (x_i - m)^2
        end
        empirical_fourth_centered_moment /= n
        empirical_variance /= n
        return (empirical_fourth_centered_moment / (empirical_variance^2)) - 3.0
    end

    # maximum likelihood estimate of kurtosis
    kurtosis(v::AbstractVector) = kurtosis(v, mean(v))

    # distance matrix
    function distances(m::AbstractMatrix)
        n = size(m, 1)
        d = Array(Float64, n, n)
        for i in 1:n
            d[i, i] = 0.0
            for j in (i + 1):n
                x = norm(m[i, :] - m[j, :])
                d[i, j] = x
                d[j, i] = x
            end
        end
        return d
    end

    # order (aka, rank), resolving ties using the mean rank
    function tiedrank(v::AbstractArray)
        n     = length(v)
        place = sortperm(v)
        ord   = Array(Float64, n)

        i = 1
        while i <= n
            j = i
            while j + 1 <= n && v[place[i]] == v[place[j + 1]]
                j += 1
            end

            if j > i
                m = sum(i:j) / (j - i + 1)
                for k = i:j
                    ord[place[k]] = m
                end
            else
                ord[place[i]] = i
            end

            i = j + 1
        end

        return ord
    end
    tiedrank(X::AbstractMatrix) = tiedrank(reshape(X, length(X)))
    function tiedrank(X::AbstractMatrix, dim::Int)
        retmat = apply(hcat, amap(tiedrank, X, 3 - dim))
        return dim == 1 ? retmat : retmat'
    end

    #
    # spearman covariance functions
    #

    # spearman covariance between two vectors
    function cov_spearman(x::AbstractVector, y::AbstractVector, corrected::Bool)
        return cov(tiedrank(x), tiedrank(y), corrected)
    end

    # spearman covariance over all pairs of columns of two matrices
    function cov_spearman(X::AbstractMatrix, Y::AbstractMatrix, corrected::Bool)
        return [cov_spearman(X[:,i], Y[:,j], corrected) for i = 1:size(X, 2), j = 1:size(Y,2)]
    end
    function cov_spearman(x::AbstractVector, Y::AbstractMatrix, corrected::Bool)
        return [cov_spearman(x, Y[:,i], corrected) for i = 1:size(Y, 2)]
    end
    function cov_spearman(X::AbstractMatrix, y::AbstractVector, corrected::Bool)
        return [cov_spearman(X[:,i], y, corrected) for i = 1:size(X, 2)]
    end
    # spearman covariance over all pairs of columns of a matrix
    cov_spearman(X::AbstractMatrix, corrected::Bool) = cov(tiedrank(X, 1), corrected)

    cov_spearman(x) = cov_spearman(x, true)
    cov_spearman(x, y) = cov_spearman(x, y, true)

    #
    # spearman correlation functions
    #

    # spearman correlation between two vectors
    function cor_spearman(x::AbstractVector, y::AbstractVector, corrected::Bool)
        return cor(tiedrank(x), tiedrank(y), corrected)
    end

    # spearman correlation over all pairs of columns of two matrices
    function cor_spearman(X::AbstractMatrix, Y::AbstractMatrix, corrected::Bool)
        return cor(tiedrank(X, 1), tiedrank(Y, 1))
    end
    function cor_spearman(X::AbstractMatrix, y::AbstractVector, corrected::Bool)
        return cor(tiedrank(X, 1), tiedrank(y))
    end
    function cor_spearman(x::AbstractVector, Y::AbstractMatrix, corrected::Bool)
        return cor(tiedrank(x), tiedrank(Y, 1))
    end

    # spearman correlation over all pairs of columns of a matrix
    function cor_spearman(X::AbstractMatrix, corrected::Bool)
        return cor(tiedrank(X, 1), corrected)
    end

    cor_spearman(x) = cor_spearman(x, true)
    cor_spearman(x, y) = cor_spearman(x, y, true)

    # autocorrelation at a specific lag
    function autocor(v::AbstractVector, lag::Real)
        return cor(v[1:end-lag], v[1+lag:end])
    end

    # autocorrelation at a default lag of 1
    autocor(v::AbstractVector) = autocor(v, 1)

    quantile(v::AbstractVector) = quantile(v, [.0, .25, .5, .75, 1.0])
    percentile(v::AbstractVector) = quantile(v, [1:99] / 100)
    quartile(v::AbstractVector) = quantile(v, [.25, .5, .75])
    quintile(v::AbstractVector) = quantile(v, [.2, .4, .6, .8])
    decile(v::AbstractVector) = quantile(v, [.1, .2, .3, .4, .5, .6, .7, .8, .9])
    iqr(v::AbstractVector) = quantile(v, [.25, .75])

    # run-length encoding
    function rle{T}(v::Vector{T})
        n = length(v)
        current_value = v[1]
        current_length = 1
        values = Array(T, n)
        total_values = 1
        lengths = Array(Int, n)
        total_lengths = 1
        for i in 2:n
            if v[i] == current_value
                current_length += 1
            else
                values[total_values] = current_value
                total_values += 1
                lengths[total_lengths] = current_length
                total_lengths += 1
                current_value = v[i]
                current_length = 1
            end
        end
        values[total_values] = current_value
        lengths[total_lengths] = current_length
        return (values[1:total_values], lengths[1:total_lengths])
    end

    # inverse run-length encoding
    function inverse_rle{T}(values::Vector{T}, lengths::Vector{Int})
        total_n = sum(lengths)
        pos = 0
        res = Array(T, total_n)
        n = length(values)
        for i in 1:n
            v = values[i]
            l = lengths[i]
            for j in 1:l
                pos += 1
                res[pos] = v
            end
        end
        return res
    end

    ## Empirical cummulative density function
    function ecdf{T}(X::AbstractVector{T})
        Xs = sort(X)
        isnan(Xs[end]) && error("ecdf undefined in presence of NaNs")
        n = length(X)
        e(x::Real) = searchsortedlast(Xs, x) / n
        function e(v::Vector)
            ord = sortperm(v)
            m = length(v)
            r = Array(T, m)
            r0 = 0
            i = 1
            for x in Xs
                if x > v[ord[i]]
                    r[ord[i]] = r0
                    i += 1
                end
                r0 += 1
                if i > m break end
            end
            if i == m r[ord[i]] = n end
            return r / n
        end
        return e
    end
end
