module Stats
    import Base.quantile

    export autocor,
           cor_spearman,
           cov_spearman,
           decile,
           distances,
           ecdf,
           inverse_rle,
           iqr,
           invlogit,
           kurtosis,
           logit,
           mad,
           mad!,
           mad!!,
           percentile,
           quantile,
           quartile,
           quintile,
           rle,
           skewness,
           tiedrank,
           weighted_mean,
           gmean,
           hmean,
           findat,
           table,
           range,
           variation,
           describe,
           sem

    # Weighted mean
    # NB: Weights should include 1/n factor
    function weighted_mean(v::AbstractArray, w::AbstractArray)
        sv, sw = 0.0, 0.0
        for i in 1:length(v)
            sv += v[i] * w[i]
            sw += w[i]
        end
        return sv / sw
    end

    # TODO: Support slicing along any dimensions
    function gmean(a::AbstractArray)
        s = 0.0
        n = length(a)
        for i in 1:n
            tmp = a[i]
            if tmp < 0.0
                throw(DomainError())
            elseif tmp == 0.0
                return 0.0
            else
                s += log(tmp)
            end
        end
        return exp(s / n)
    end

    # TODO: Support slicing along any dimensions
    function hmean(a::AbstractArray)
        s = 0.0
        n = length(a)
        for i in 1:n
            s += 1 / a[i]
        end
        return n / s
    end

    # median absolute deviation with consistency adjustment
    mad(v::AbstractArray, center::Number) = 1.4826 * median!(abs(v-center))
    mad(v::AbstractArray) = mad!!(copy(v))

    function mad!!(v::AbstractVector)
        center = median!(v)
        for i in 1:length(v)
            v[i] = abs(v[i]-center)
        end
        1.4826 * median!(v, checknan=false)
    end

    function mad!(v::AbstractVector)
        center = median!(v)
        by = x->abs(x-center)
        o = Sort.By(by)
        n = length(v)
        1.4826 * (isodd(n) ?
            by(select!(v,div(n+1,2),o)) :
           (by(select!(v,div(n,2),o))+by(select!(v,div(n,2)+1,o)))/2
        )
    end

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
    # Assumes vectors are column vectors
    function distances(m::AbstractMatrix)
        p, n = size(m)
        D = Array(Float64, n, n)
        for j in 1:n
            D[j, j] = 0.0
            for i in (j + 1):n
                x = 0.0
                for k in 1:p
                    x += (m[k, i] - m[k, j])^2
                end
                x = sqrt(x)
                D[i, j] = x
                D[j, i] = x
            end
        end
        return D
    end

    # distance matrix between entries of two distinct matrices
    # Assumes vectors are column vectors
    function distances(a::AbstractMatrix, b::AbstractMatrix)
        p_a, n_a = size(a)
        p_b, n_b = size(b)
        if p_a != p_b
            throw(DomainError("A and B must have the same dimension"))
        end
        D = Array(Float64, n_a, n_b)
        for j in 1:n_b
            for i in 1:n_a
                x = 0.0
                for k in 1:p_a
                    x += (a[k, i] - b[k, j])^2
                end
                D[i, j] = sqrt(x)
            end
        end
        return D
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
        retmat = apply(hcat, mapslices(tiedrank, X, 3 - dim))
        return dim == 1 ? retmat : retmat'
    end

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

    function logit(p::Real)
        if p < 0.0 || p > 1.0
            DomainError("logit(p) only defined when p lies in [0, 1]")
        else
            return log(p / (1.0 - p))
        end
    end

    invlogit(z::Real) = 1.0 / (1.0 + exp(-clamp(z, -709.0, 745.0)))


    # TODO: Support slicing along any dimensions
    function modes{T}(a::AbstractArray{T})
        res = Array(T, 0)
        tab = table(a)
        m = max(values(tab))
        for (k, v) in tab
            if v == m
                push!(res, k)
            end
        end
        return res
    end

    # TODO: Support slicing along any dimensions
    function findat!{T}(indices::Vector{Int},
                        a::AbstractArray,
                        b::AbstractArray{T})
        inds = Dict{T, Int}()
        for i in 1:length(b)
            tmp = b[i]
            if !haskey(inds, tmp)
                inds[tmp] = i
            end
        end
        for i = 1:length(a)
            indices[i] = get(inds, a[i], 0)
        end
        return
    end

    # TODO: Support slicing along any dimensions
    function findat(a::AbstractArray, b::AbstractArray)
        indices = Array(Int, length(a))
        findat!(indices, a, b)
        return indices
    end

    # TODO: Support slicing along any dimensions
    function table{T}(a::AbstractArray{T})
        counts = Dict{T, Int}()
        for i = 1:length(a)
            tmp = a[i]
            counts[tmp] = get(counts, tmp, 0) + 1
        end
        return counts
    end

    # TODO: Support slicing along any dimensions
    function range{T <: Real}(a::AbstractArray{T})
        minval, maxval = typemax(T), typemin(T)
        for i in 1:length(a)
            tmp = a[i]
            if tmp < minval
                minval = tmp
            end
            if tmp > maxval
                maxval = tmp
            end
        end
        return minval, maxval
    end

    # TODO: Support slicing along any dimensions
    variation(a::AbstractArray) = std(a) / mean(a)

    # TODO: Support slicing along any dimensions
    function describe(a::AbstractArray)
        q00, q25, q50, q75, q10 = quantile(a, [0.00, 0.25, 0.50, 0.75, 1.00])
        @printf "Min:          %.6f\n" q00
        @printf "1st Quartile: %.6f\n" q25
        @printf "Median:       %.6f\n" q50
        @printf "Mean:         %.6f\n" mean(a)
        @printf "3rd Quartile: %.6f\n" q75
        @printf "Max:          %.6f\n" q10
        return
    end

    # TODO: Support slicing along any dimensions
    sem(a::AbstractArray) = sqrt(var(a) / length(a))

    # TODO: Support slicing along any dimensions
    function midrange(a::AbstractArray)
        minval, maxval = range(a)
        return (maxvavl - minval) / 2
    end

end # module
