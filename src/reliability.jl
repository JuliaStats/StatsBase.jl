using LinearAlgebra: Diagonal, isposdef

struct Reliability
    alpha::AbstractFloat
    dropped::Vector{Pair{Int, AbstractFloat}}
end

function Base.show(io::IO, x::Reliability)
    println(io, "Reliability for all items: $(round(x.alpha; digits=4))")
    println(io, "")
    println(io, "Reliability if an item is dropped:")
    for i ∈ x.dropped
        println(io, "item $(i.first): $(round(i.second; digits=4))")
    end
end

function _crombach_alpha(covmatrix::AbstractMatrix{T}) where T <: Real
    @assert isposdef(covmatrix) "Covariance matrix is not positive definite!"
    k = size(covmatrix, 2)
    @assert k > 1 "Covariance matrix has only one columnn!"
    σ = sum(covmatrix)
    ρ = k / (k - 1) * (1 - sum(i -> covmatrix[i, i], 1:k)/ σ)
    return ρ
end

"""
    crombach_alpha(covmatrix::AbstractMatrix{T}) where T <: Real

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
Reliability for all items: 0.814

Reliability if an item is dropped:
item 1: 0.75
item 2: 0.761
item 3: 0.771
item 4: 0.783
```
"""
function crombach_alpha(covmatrix::AbstractMatrix{T}) where T <: Real
    alpha = _crombach_alpha(covmatrix)
    k = size(covmatrix, 2)
    dropped = Vector{Pair{Int, AbstractFloat}}(undef, k)
    for i ∈ 1:k
        reduced_covmatrix = covmatrix[1:end .!= i, 1:end .!= i]
        dropped[i] = Pair{Int, AbstractFloat}(i, _crombach_alpha(reduced_covmatrix))
    end
    return Reliability(alpha, dropped)
end
