# Empirical estimation of CDF and PDF


## Empirical CDF

struct ECDF{T <: AbstractVector{<:Real}, W <: AbstractWeights{<:Real}}
    sorted_values::T
    weights::W
end

function (ecdf::ECDF)(x::Real)
    weightsum = sum(ecdf.weights)
    n = searchsortedlast(ecdf.sorted_values, x)
    sum(view(ecdf.weights, 1:n)) / weightsum
end

function (ecdf::ECDF)(v::RealVector)
    ord = sortperm(v)
    m = length(v)
    r = similar(ecdf.sorted_values, m)
    r0 = 0
    i = 1
    weightsum = sum(ecdf.weights)
    for (j, x) in enumerate(ecdf.sorted_values)
        while i <= m && x > v[ord[i]]
            r[ord[i]] = r0
            i += 1
        end
        r0 += ecdf.weights[j]
        if i > m
            break
        end
    end
    while i <= m
        r[ord[i]] = weightsum
        i += 1
    end
    return r / weightsum
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
function ecdf(X::RealVector{T}, w::AbstractWeights{W} = weights(ones(length(X)))) where {T, W <: Real}
    length(X) == length(w) || throw(ArgumentError("data and weight vectors must be the same size," *
        "got $(length(X)) and $(length(w))"))
    ord = sortperm(X)
    ECDF(X[ord], weights(w.values[ord]))
end

minimum(ecdf::ECDF) = first(ecdf.sorted_values)

maximum(ecdf::ECDF) = last(ecdf.sorted_values)

extrema(ecdf::ECDF) = (minimum(ecdf), maximum(ecdf))
