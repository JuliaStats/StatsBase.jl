using LinearAlgebra: isposdef

struct Reliability{T <: Real}
    alpha::T
    dropped::Vector{Pair{Int, T}}
end

function Base.show(io::IO, x::Reliability)
    @printf(io, "Reliability for all items: %.4f", x.alpha)
    isempty(x.dropped) && return
    println(io, "")
    println(io, "Reliability if an item is dropped:")
    for i ∈ x.dropped
        @printf(io, "item %i: %.4f\n, i.first, i.second")
    end
end

"""
    crombach_alpha(covmatrix::AbstractMatrix{<:Real})

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

julia> crombach_alpha(cov_X)
Reliability for all items: 0.8136

Reliability if an item is dropped:
item 1: 0.75
item 2: 0.7606
item 3: 0.7714
item 4: 0.7836
```
"""
function crombach_alpha(covmatrix::AbstractMatrix{<:Real})
    isposdef(covmatrix) || throw(ArgumentError("Covariance matrix must be positive definite."))
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
        dropped = Vector{Pair{Int, typeof(alpha)}}(undef, k)
        for i ∈ 1:k
            dropped[i] = i => (k - 1) * (1 - (σ_diag - covmatrix[i,i]) / (σ - 2*v[i] - covmatrix[i,i])) / (k - 2)
        end
    else
        # if k = 2 do not produce dropped; this has to be also
        # correctly handled in show
        dropped = Pair{Int, AbstractFloat}[]
    end
    return Reliability(alpha, dropped)
end
