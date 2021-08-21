struct Reliability{T <: Real}
    alpha::T
    dropped::Vector{T}
end

function Base.show(io::IO, x::Reliability)
    @printf(io, "Reliability for all items: %.4f", x.alpha)
    isempty(x.dropped) && return
    println(io, "")
    println(io, "Reliability if an item is dropped:")
    for (idx, val) in enumerate(x.dropped)
        @printf(io, "item %i: %.4f\n", idx, val)
    end
end

"""
    crombachalpha(covmatrix::AbstractMatrix{<:Real})

Calculate Crombach's alpha (1951) from a covariance matrix `covmatrix` according to
the [formula](https://en.wikipedia.org/wiki/Cronbach%27s_alpha):

```math
\\rho = \\frac{k}{k-1} (1 - \\frac{\\sum^k_{i=1} \\sigma^2_i}{\\sigma^2_X})
```

where ``k`` is the number of items, i.e. columns; ``\\sigma_i`` denotes item variance;
and ``\\sigma^2_i`` consists of item variances and inter-item covariances.

Returns a `Reliability` object that holds:

* `alpha`: the reliability score for all items, i.e. columns, in `covmatrix`; and
* `dropped`: A `Vector{Pair{item, score}}` giving reliability scores if a specific item,
  i.e. column, is dropped from `covmatrix`.

# Example
```jldoctest
julia> using StatsBase

julia> cov_X = [10 6 6 6;
                6 11 6 6;
                6 6 12 6;
                6 6 6 13];

julia> crombachalpha(cov_X)
Reliability for all items: 0.8136

Reliability if an item is dropped:
item 1: 0.7500
item 2: 0.7606
item 3: 0.7714
item 4: 0.7836
```
"""
function crombachalpha(covmatrix::AbstractMatrix{<:Real})
    isposdef(covmatrix) || throw(ArgumentError("Covariance matrix must be positive definite."))
    k = size(covmatrix, 2)
    k >= 1  || throw(ArgumentError("Covariance matrix must have more than one column."))
    v = vec(sum(covmatrix, dims=1))
    σ = sum(v)
    for i in axes(v, 1)
        v[i] -= covmatrix[i, i]
    end
    σ_diag = sum(i -> covmatrix[i, i], 1:k)

    alpha = k * (1 - σ_diag / σ) / (k - 1)
    if k > 2
        dropped = Vector{typeof(alpha)}(undef, k)
        for i in eachindex(dropped)
            dropped[i] = (k - 1) * (1 - (σ_diag - covmatrix[i,i]) / (σ - 2*v[i] - covmatrix[i,i])) / (k - 2)
        end
    else
        # if k = 2 do not produce dropped; this has to be also
        # correctly handled in show
        dropped = Vector{typeof(alpha)}()
    end
    return Reliability(alpha, dropped)
end
