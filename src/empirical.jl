# Empirical estimation of CDF and PDF


## Empirical CDF

"""
    ecdf(X::RealVector)

Return an empirical cumulative distribution function based on a vector of samples given in `X`.
"""
function ecdf{T<:Real}(X::RealVector{T})
    Xs = sort(X)
    n = length(X)

    ef(x::Real) = searchsortedlast(Xs, x) / n

    function ef(v::RealVector)
        ord = sortperm(v)
        m = length(v)
        r = Array(T, m)
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


