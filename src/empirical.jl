# Empirical estimation of CDF and PDF


## Empirical CDF
"""
    ecdf(X)

Return an empirical cumulative distribution function (ECDF) based on a vector of samples
given in `X`.

Note: this is a higher-level function that returns a function, which can then be applied
to evaluate CDF values on other samples.
"""
function ecdf(X::RealVector{T}) where T<:Real
    Xs = sort(X)
    n = length(X)

    ef(x::Real) = searchsortedlast(Xs, x) / n

    function ef(v::RealVector)
        ord = sortperm(v)
        m = length(v)
        r = Vector{T}(undef, m)
        r0 = 0

        i = 1
        for x in Xs
            while i <= m && x > v[ord[i]]
                r[ord[i]] = r0
                i += 1
            end
            r0 += 1
            if i > m
                break
            end
        end
        while i <= m
            r[ord[i]] = n
            i += 1
        end
        return r / n
    end

    return ef
end


