# Temporary location for method definitions to make transformations work for 
# CuArrays with `CUDA.allowscalar(false)`.
using CUDA

## Changes required for `ZScoreTransform`
# Proposal: Chance in CUDA.jl as described.

# TODO Implement instead of `Statistics._var` in Statistics._var in CUDA/src/statistics.jl
# required to make `mean_and_std` work for `CuArray`s.
Statistics.varm(A::CuArray{<:Real},m::AbstractArray{<:Real}; dims, corrected::Bool=true) =
    sum((A .- m).^2, dims=dims)/(prod(size(A)[[dims...]])::Int-corrected)

# TODO Implement instead of `Statistics._std` in CUDA/src/statistics.jl
# required to make `mean_and_std` work for `CuArray`s.
Statistics.stdm(A::CuArray{<:Real},m::AbstractArray{<:Real}, dim::Int; corrected::Bool=true) =
    sqrt.(varm(A,m;dims=dim,corrected=corrected))

## Required for `UnitRangeTransform`
# How to integrate `StatsBase._compute_extrema` with StatsBase.jl or CUDA.jl? Options:
# 1) Add to CUDA.jl. Requires adding StatsBase as a dependency to CUDA.jl
# 2) Make this the generic version in StatsBase._compute_extrema, specialize the current implementation
#    for `Array`s?
# 3) Add to StatsBase.jl. Requires adding CUDA.jl as a (significant!) dependency of StatsBase

function StatsBase._compute_extrema(X::Union{CuMatrix{<:Real},Adjoint{T,CuMatrix{T}} where T<:Real})
    l = size(X,2)
    tmin = minimum(X, dims=1)[:]
    tmax = maximum(X, dims=1)[:]
    return l, tmin, tmax
end
