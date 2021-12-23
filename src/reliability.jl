struct CronbachAlpha{T <: Real}
    alpha::T
    dropped::Vector{T}
end

function Base.show(io::IO, x::CronbachAlpha)
    @printf(io, "Cronbach's alpha for all items: %.4f\n", x.alpha)
    isempty(x.dropped) && return
    println(io, "\nCronbach's alpha if an item is dropped:")
    for (idx, val) in enumerate(x.dropped)
        @printf(io, "item %i: %.4f\n", idx, val)
    end
end

"""
    cronbachalpha(covmatrix::AbstractMatrix{<:Real})

Calculate Cronbach's alpha (1951) from a covariance matrix `covmatrix` according to
the [formula](https://en.wikipedia.org/wiki/Cronbach%27s_alpha):

```math
\\rho = \\frac{k}{k-1} (1 - \\frac{\\sum^k_{i=1} \\sigma^2_i}{\\sum_{i=1}^k \\sum_{j=1}^k \\sigma_{ij}})
```

where ``k`` is the number of items, i.e. columns, ``\\sigma_i^2`` the item variance,
and ``\\sigma_{ij}`` the inter-item covariance.

Returns a `CronbachAlpha` object that holds:

* `alpha`: the Cronbach's alpha score for all items, i.e. columns, in `covmatrix`; and
* `dropped`: a vector giving Cronbach's alpha scores if a specific item,
  i.e. column, is dropped from `covmatrix`.

# Example
```jldoctest
julia> using StatsBase

julia> cov_X = [10 6 6 6;
                6 11 6 6;
                6 6 12 6;
                6 6 6 13];

julia> cronbachalpha(cov_X)
Cronbach's alpha for all items: 0.8136

Cronbach's alpha if an item is dropped:
item 1: 0.7500
item 2: 0.7606
item 3: 0.7714
item 4: 0.7826
```
"""
function cronbachalpha(covmatrix::AbstractMatrix{<:Real})
    if !isposdef(covmatrix)
        throw(ArgumentError("Covariance matrix must be positive definite. " *
                            "Maybe you passed the data matrix instead of its covariance matrix? " *
                            "If so, call `cronbachalpha(cov(...))` instead."))
    end
    k = size(covmatrix, 2)
    k > 1  || throw(ArgumentError("Covariance matrix must have more than one column."))
    v = vec(sum(covmatrix, dims=1))
    σ = sum(v)
    for i in axes(v, 1)
        v[i] -= covmatrix[i, i]
    end
    σ_diag = sum(i -> covmatrix[i, i], 1:k)

    alpha = k * (1 - σ_diag / σ) / (k - 1)
    if k > 2
        dropped = typeof(alpha)[(k - 1) * (1 - (σ_diag - covmatrix[i, i]) / (σ - 2*v[i] - covmatrix[i, i])) / (k - 2)
                                for i in 1:k]
    else
        # if k = 2 do not produce dropped; this has to be also
        # correctly handled in show
        dropped = Vector{typeof(alpha)}()
    end
    return CronbachAlpha(alpha, dropped)
end
