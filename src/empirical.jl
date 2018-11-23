# Empirical estimation of CDF and PDF


## Empirical CDF

struct ECDF{T <: AbstractVector{<:Real}}
    sorted_values::T
end

function (ecdf::ECDF)(x::Real)
    searchsortedlast(ecdf.sorted_values, x) / length(ecdf.sorted_values)
end

function (ecdf::ECDF)(v::RealVector)
    ord = sortperm(v)
    m = length(v)
    r = similar(ecdf.sorted_values, m)
    r0 = 0
    i = 1
    n = length(ecdf.sorted_values)
    for x in ecdf.sorted_values
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

"""
    ecdf(X)

Return an empirical cumulative distribution function (ECDF) based on a vector of samples
given in `X`.

Note: this function that returns a callable composite type, which can then be applied to
evaluate CDF values on other samples.

`extrema`, `minimum`, and `maximum` are supported to for obtaining the range over which
function is inside the interval ``(0,1)``; the function is defined for the whole real line.
"""
ecdf(X::RealVector{T}) where T<:Real = ECDF(sort(X))

minimum(ecdf::ECDF) = first(ecdf.sorted_values)

maximum(ecdf::ECDF) = last(ecdf.sorted_values)

extrema(ecdf::ECDF) = (minimum(ecdf), maximum(ecdf))
