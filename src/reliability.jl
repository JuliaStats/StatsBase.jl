using LinearAlgebra: Diagonal, isposdef

struct Reliability{T <: Real}
    alpha::T
    dropped::Vector{Pair{Int, T}}
end

function Base.show(io::IO, x::Reliability)
    println(io, "Reliability for all items: $(round(x.alpha; digits=4))")
    isempty(x.dropped) && return
    println(io, "")
    println(io, "Reliability if an item is dropped:")
    for i ∈ x.dropped
        println(io, "item $(i.first): $(round(i.second; digits=4))")
    end
end

"""
    crombach_alpha(covmatrix::AbstractMatrix{<:Real})

Calculate Crombach's alpha (1951) from a covariance matrix `covmatrix` according to
the Wikipedia formula (https://en.wikipedia.org/wiki/Cronbach%27s_alpha):

```math
\\rho = \\frac{k}{k-1} (1 - \\frac{\\sum^k_{i=1} \\sigma^2_i}{\\sigma^2_X})
```

where k is the number of items, i.e. columns; σᵢ denotes item variance;
and σ²ₓ consists of item variances and inter-item covariances.

Returns a `Reliability` object that holds:

* `alpha`: the reliability score for all items, i.e. columns, in `comatrix`; and
* `dropped`: A `Vector{Pair{item, score}}` reliability scores if a specific item, i.e. column/row, is dropped from `comatrix`.

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
    isposdef(covmatrix) || throw(ArgumentError("Covariance matrix is not positive definite!"))
    k = size(covmatrix, 2)
    k > 1  || throw(ArgumentError("Covariance matrix has only one columnn!"))
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
