# Empirical estimation of CDF and PDF


## Empirical CDF
"""
    ecdf(X)

Compute the empirical cumulative distribution function (ECDF)
of a real-valued vector.
"""
function ecdf{T<:Real}(X::RealVector{T})
    Xs = sort(X)
    n = length(X)

    ef(x::Real) = searchsortedlast(Xs, x) / n

    function ef(v::RealVector)
        ord = sortperm(v)
        m = length(v)
        r = Vector{T}(m)
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


## Empirical CCDF

function eccdf{T<:Real}(X::RealVector{T})
    Xs = sort(X)
    n = length(X)

    ef(x::Real) = (n - searchsortedfirst(Xs, x) + 1) / n

    function ef(v::Vector)
        ord = sortperm(v)
        m = length(v)
        r = Array(T, m)
        r0 = 0

        i = 1
        for x in Xs
            while i <= m && x >= v[ord[i]]
                r[ord[i]] = n - r0
                i += 1
            end
            r0 += 1
            if i > m 
            	break 
            end
        end
        while i <= m
            r[ord[i]] = 1
            i += 1
        end
        return r / n
    end

    return ef
end
