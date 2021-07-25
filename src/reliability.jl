using LinearAlgebra: Diagonal

struct Reliability
    alpha::Float64
    dropped::Vector{Pair{Int64, Float64}}
end

function Base.show(io::IO, x::Reliability)
    println(io, "Reliability for all items: $(round(x.alpha; digits=3))")
    println(io, "")
    println(io, "Reliability if an item is dropped:")
    for i ∈ x.dropped
        println(io, "item $(i.first): $(round(i.second; digits=3))")
    end
end

function _alpha(covmatrix::AbstractMatrix{T}) where T <: Real
    k = size(covmatrix, 2)
    σ = sum(covmatrix)
    σ_ij = sum(covmatrix - Diagonal(covmatrix)) / (k * (k - 1))
    ρ = k^2 * σ_ij / σ
    return ρ
end

"""
    alpha(covmatrix::AbstractMatrix{T}) where T <: Real

Calculate Crombach's alpha (1951) from a covariance matrix `covmatrix` according to
the Wikipedia formula (https://en.wikipedia.org/wiki/Cronbach%27s_alpha):

```math
\\rho = \\frac{k^2 \\bar{sigma}_{ij}}{\\sigma^2_X}
```

where k is the number of items, i.e. columns; \\bar{sigma}_{ij} denote the average of
the inter-item covariances; and \\sigma^2_X consists of item variances and inter-item
covariances.

Returns a `Reliability` object that holds:

* `alpha`: the reliability score for all items, i.e. columns, in `comatrix`; and
* `dropped`: reliability score if a specific item, i.e. column, is dropped from `comatrix`.

# Example
```jldoctest
julia> using StatsBase

julia> cov_X = [10 6 6 6;
                6 11 6 6;
                6 6 12 6;
                6 6 6 13];

julia> alpha(cov_X)
Reliability for all items: 0.814

Reliability if an item is dropped:
item 1: 0.75
item 2: 0.761
item 3: 0.771
item 4: 0.783
```
"""
function alpha(covmatrix::AbstractMatrix{T}) where T <: Real
    alpha = _alpha(covmatrix)
    k = size(covmatrix, 2)
    dropped = Vector{Pair{Int64, Float64}}(undef, k)
    for i ∈ 1:k
        reduced_covmatrix = covmatrix[1:end .!= i, 1:end .!= i]
        @inbounds dropped[i] = Pair{Int64, Float64}(i, _alpha(reduced_covmatrix))
    end
    return Reliability(alpha, dropped)
end
