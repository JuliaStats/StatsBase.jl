# Empirical estimation of CDF and PDF

## Empirical CDF

struct ECDF{T <: AbstractVector{<:Real}, W <: AbstractWeights{<:Real}}
    sorted_values::T
    weights::W
end

function (ecdf::ECDF)(x::Real)
    isnan(x) && return NaN
    n = searchsortedlast(ecdf.sorted_values, x)
    evenweights = isempty(ecdf.weights)
    weightsum = evenweights ? length(ecdf.sorted_values) : sum(ecdf.weights)
    partialsum = evenweights ? n : sum(view(ecdf.weights, 1:n))
    partialsum / weightsum
end

function (ecdf::ECDF)(v::RealVector)
    evenweights = isempty(ecdf.weights)
    weightsum = evenweights ? length(ecdf.sorted_values) : sum(ecdf.weights)
    ord = sortperm(v)
    m = length(v)
    r = similar(ecdf.sorted_values, m)
    r0 = zero(weightsum)
    i = 1
    for (j, x) in enumerate(ecdf.sorted_values)
        while i <= m && x > v[ord[i]]
            r[ord[i]] = r0
            i += 1
        end
        r0 += evenweights ? 1 : ecdf.weights[j]
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
    ecdf(X; weights::AbstractWeights)

Return an empirical cumulative distribution function (ECDF) based on a vector of samples
given in `X`. Optionally providing `weights` returns a weighted ECDF.

Note: this function that returns a callable composite type, which can then be applied to
evaluate CDF values on other samples.

`extrema`, `minimum`, and `maximum` are supported to for obtaining the range over which
function is inside the interval ``(0,1)``; the function is defined for the whole real line.
"""
function ecdf(X::RealVector; weights::AbstractVector{<:Real}=Weights(Float64[]))
    any(isnan, X) && throw(ArgumentError("ecdf can not include NaN values"))
    isempty(weights) || length(X) == length(weights) || throw(ArgumentError("data and weight vectors must be the same size," *
        "got $(length(X)) and $(length(weights))"))
    ord = sortperm(X)
    ECDF(X[ord], isempty(weights) ? weights : Weights(weights[ord]))
end

minimum(ecdf::ECDF) = first(ecdf.sorted_values)

maximum(ecdf::ECDF) = last(ecdf.sorted_values)

extrema(ecdf::ECDF) = (minimum(ecdf), maximum(ecdf))
