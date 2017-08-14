
###### Weight vector #####

abstract type AbstractWeights{S<:Real, T<:Real, V<:AbstractVector{T}} <: AbstractVector{T} end

"""
    @weights name

Generates a new generic weight type with specified `name`, which subtypes `AbstractWeights`
and stores the `values` (`V<:RealVector`) and `sum` (`S<:Real`).
"""
macro weights(name)
    return quote
        struct $name{S<:Real, T<:Real, V<:AbstractVector{T}} <: AbstractWeights{S, T, V}
            values::V
            sum::S
        end
    end
end

eltype(wv::AbstractWeights) = eltype(wv.values)
length(wv::AbstractWeights) = length(wv.values)
values(wv::AbstractWeights) = wv.values
sum(wv::AbstractWeights) = wv.sum
isempty(wv::AbstractWeights) = isempty(wv.values)

Base.getindex(wv::AbstractWeights, i) = getindex(wv.values, i)
Base.size(wv::AbstractWeights) = size(wv.values)

"""
    varcorrection(n::Integer, corrected=false)

Compute a bias correction factor for calculating `var`, `std` and `cov` with
`n` observations. Returns ``\\frac{1}{n - 1}`` when `corrected=true`
(i.e. [Bessel's correction](https://en.wikipedia.org/wiki/Bessel's_correction)),
otherwise returns ``\\frac{1}{n}`` (i.e. no correction).
"""
@inline varcorrection(n::Integer, corrected::Bool=false) = 1 / (n - Int(corrected))

@weights Weights

"""
    Weights(vs, wsum=sum(vs))

Construct a `Weights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

The `Weights` type describes a generic weights vector which does not support
all operations possible for [`FrequencyWeights`](@ref), [`AnalyticWeights`](@ref)
and [`ProbabilityWeights`](@ref).
"""
Weights(vs::V, s::S=sum(vs)) where {S<:Real, V<:RealVector} = Weights{S, eltype(vs), V}(vs, s)

"""
    weights(vs)

Construct a `Weights` vector from array `vs`.
See the documentation for [`Weights`](@ref) for more details.
"""
weights(vs::RealVector) = Weights(vs)
weights(vs::RealArray) = Weights(vec(vs))

"""
    varcorrection(w::Weights, corrected=false)

Returns ``\\frac{1}{\\sum w}`` when `corrected=false` and throws an `ArgumentError`
if `corrected=true`.
"""
@inline function varcorrection(w::Weights, corrected::Bool=false)
    corrected && throw(ArgumentError("Weights type does not support bias correction: " *
                                     "use FrequencyWeights, AnalyticWeights or ProbabilityWeights if applicable."))
    1 / w.sum
end

@weights AnalyticWeights

"""
    AnalyticWeights(vs, wsum=sum(vs))

Construct an `AnalyticWeights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

Analytic weights describe a non-random relative importance (usually between 0 and 1)
for each observation. These weights may also be referred to as reliability weights,
precision weights or inverse variance weights. These are typically used when the observations
being weighted are aggregate values (e.g., averages) with differing variances.
"""
AnalyticWeights(vs::V, s::S=sum(vs)) where {S<:Real, V<:RealVector} =
    AnalyticWeights{S, eltype(vs), V}(vs, s)

"""
    aweights(vs)

Construct an `AnalyticWeights` vector from array `vs`.
See the documentation for [`AnalyticWeights`](@ref) for more details.
"""
aweights(vs::RealVector) = AnalyticWeights(vs)
aweights(vs::RealArray) = AnalyticWeights(vec(vs))

"""
    varcorrection(w::AnalyticWeights, corrected=false)

* `corrected=true`: ``\\frac{1}{\\sum w - \\sum {w^2} / \\sum w}``
* `corrected=false`: ``\\frac{1}{\\sum w}``
"""
@inline function varcorrection(w::AnalyticWeights, corrected::Bool=false)
    s = w.sum

    if corrected
        sum_sn = sum(x -> (x / s) ^ 2, w)
        1 / (s * (1 - sum_sn))
    else
        1 / s
    end
end

@weights FrequencyWeights

"""
    FrequencyWeights(vs, wsum=sum(vs))

Construct a `FrequencyWeights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

Frequency weights describe the number of times (or frequency) each observation
was observed. These weights may also be referred to as case weights or repeat weights.
"""
FrequencyWeights(vs::V, s::S=sum(vs)) where {S<:Real, V<:RealVector} =
    FrequencyWeights{S, eltype(vs), V}(vs, s)

"""
    fweights(vs)

Construct a `FrequencyWeights` vector from a given array.
See the documentation for [`FrequencyWeights`](@ref) for more details.
"""
fweights(vs::RealVector) = FrequencyWeights(vs)
fweights(vs::RealArray) = FrequencyWeights(vec(vs))

"""
    varcorrection(w::FrequencyWeights, corrected=false)

* `corrected=true`: ``\\frac{1}{\\sum{w} - 1}``
* `corrected=false`: ``\\frac{1}{\\sum w}``
"""
@inline function varcorrection(w::FrequencyWeights, corrected::Bool=false)
    s = w.sum

    if corrected
        1 / (s - 1)
    else
        1 / s
    end
end

@weights ProbabilityWeights

"""
    ProbabilityWeights(vs, wsum=sum(vs))

Construct a `ProbabilityWeights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

Probability weights represent the inverse of the sampling probability for each observation,
providing a correction mechanism for under- or over-sampling certain population groups.
These weights may also be referred to as sampling weights.
"""
ProbabilityWeights(vs::V, s::S=sum(vs)) where {S<:Real, V<:RealVector} =
    ProbabilityWeights{S, eltype(vs), V}(vs, s)

"""
    pweights(vs)

Construct a `ProbabilityWeights` vector from a given array.
See the documentation for [`ProbabilityWeights`](@ref) for more details.
"""
pweights(vs::RealVector) = ProbabilityWeights(vs)
pweights(vs::RealArray) = ProbabilityWeights(vec(vs))

"""
    varcorrection(w::ProbabilityWeights, corrected=false)

* `corrected=true`: ``\\frac{n}{(n - 1) \\sum w}`` where ``n`` equals `count(!iszero, w)`
* `corrected=false`: ``\\frac{1}{\\sum w}``
"""
@inline function varcorrection(w::ProbabilityWeights, corrected::Bool=false)
    s = w.sum

    if corrected
        n = count(!iszero, w)
        n / (s * (n - 1))
    else
        1 / s
    end
end

##### Weighted sum #####

## weighted sum over vectors

"""
    wsum(v, w::AbstractVector, [dim])

Compute the weighted sum of an array `v` with weights `w`, optionally over the dimension `dim`.
"""
wsum(v::AbstractVector, w::AbstractVector) = dot(v, w)
wsum(v::AbstractArray, w::AbstractVector) = dot(vec(v), w)

# Note: the methods for BitArray and SparseMatrixCSC are to avoid ambiguities
Base.sum(v::BitArray, w::AbstractWeights) = wsum(v, values(w))
Base.sum(v::SparseMatrixCSC, w::AbstractWeights) = wsum(v, values(w))
Base.sum(v::AbstractArray, w::AbstractWeights) = dot(v, values(w))

## wsum along dimension
#
#  Brief explanation of the algorithm:
#  ------------------------------------
#
#  1. _wsum! provides the core implementation, which assumes that
#     the dimensions of all input arguments are consistent, and no
#     dimension checking is performed therein.
#
#     wsum and wsum! perform argument checking and call _wsum!
#     internally.
#
#  2. _wsum! adopt a Cartesian based implementation for general
#     sub types of AbstractArray. Particularly, a faster routine
#     that keeps a local accumulator will be used when dim = 1.
#
#     The internal function that implements this is _wsum_general!
#
#  3. _wsum! is specialized for following cases:
#     (a) A is a vector: we invoke the vector version wsum above.
#         The internal function that implements this is _wsum1!
#
#     (b) A is a dense matrix with eltype <: BlasReal: we call gemv!
#         The internal function that implements this is _wsum2_blas!
#
#     (c) A is a contiguous array with eltype <: BlasReal:
#         dim == 1: treat A like a matrix of size (d1, d2 x ... x dN)
#         dim == N: treat A like a matrix of size (d1 x ... x d(N-1), dN)
#         otherwise: decompose A into multiple pages, and apply _wsum2!
#         for each
#
#     (d) A is a general dense array with eltype <: BlasReal:
#         dim <= 2: delegate to (a) and (b)
#         otherwise, decompose A into multiple pages
#

function _wsum1!(R::AbstractArray, A::AbstractVector, w::AbstractVector, init::Bool)
    r = wsum(A, w)
    if init
        R[1] = r
    else
        R[1] += r
    end
    return R
end

function _wsum2_blas!(R::StridedVector{T}, A::StridedMatrix{T}, w::StridedVector{T}, dim::Int, init::Bool) where T<:BlasReal
    beta = ifelse(init, zero(T), one(T))
    trans = dim == 1 ? 'T' : 'N'
    BLAS.gemv!(trans, one(T), A, w, beta, R)
    return R
end

function _wsumN!(R::StridedArray{T}, A::StridedArray{T,N}, w::StridedVector{T}, dim::Int, init::Bool) where {T<:BlasReal,N}
    if dim == 1
        m = size(A, 1)
        n = div(length(A), m)
        _wsum2_blas!(view(R,:), reshape(A, (m, n)), w, 1, init)
    elseif dim == N
        n = size(A, N)
        m = div(length(A), n)
        _wsum2_blas!(view(R,:), reshape(A, (m, n)), w, 2, init)
    else # 1 < dim < N
        m = 1
        for i = 1:dim-1; m *= size(A, i); end
        n = size(A, dim)
        k = 1
        for i = dim+1:N; k *= size(A, i); end
        Av = reshape(A, (m, n, k))
        Rv = reshape(R, (m, k))
        for i = 1:k
            _wsum2_blas!(view(Rv,:,i), view(Av,:,:,i), w, 2, init)
        end
    end
    return R
end

function _wsumN!(R::StridedArray{T}, A::DenseArray{T,N}, w::StridedVector{T}, dim::Int, init::Bool) where {T<:BlasReal,N}
    @assert N >= 3
    if dim <= 2
        m = size(A, 1)
        n = size(A, 2)
        npages = 1
        for i = 3:N
            npages *= size(A, i)
        end
        rlen = ifelse(dim == 1, n, m)
        Rv = reshape(R, (rlen, npages))
        for i = 1:npages
            _wsum2_blas!(view(Rv,:,i), view(A,:,:,i), w, dim, init)
        end
    else
        _wsum_general!(R, identity, A, w, dim, init)
    end
    return R
end

# General Cartesian-based weighted sum across dimensions
@generated function _wsum_general!(R::AbstractArray{RT}, f::supertype(typeof(abs)),
                                   A::AbstractArray{T,N}, w::AbstractVector{WT}, dim::Int, init::Bool) where {T,RT,WT,N}
    quote
        init && fill!(R, zero(RT))
        wi = zero(WT)
        if dim == 1
            @nextract $N sizeR d->size(R,d)
            sizA1 = size(A, 1)
            @nloops $N i d->(d>1? (1:size(A,d)) : (1:1)) d->(j_d = sizeR_d==1 ? 1 : i_d) begin
                @inbounds r = (@nref $N R j)
                for i_1 = 1:sizA1
                    @inbounds r += f(@nref $N A i) * w[i_1]
                end
                @inbounds (@nref $N R j) = r
            end
        else
            @nloops $N i A d->(if d == dim
                                   wi = w[i_d]
                                   j_d = 1
                               else
                                   j_d = i_d
                               end) @inbounds (@nref $N R j) += f(@nref $N A i) * wi
        end
        return R
    end
end

@generated function _wsum_centralize!(R::AbstractArray{RT}, f::supertype(typeof(abs)),
                                      A::AbstractArray{T,N}, w::AbstractVector{WT}, means,
                                      dim::Int, init::Bool) where {T,RT,WT,N}
    quote
        init && fill!(R, zero(RT))
        wi = zero(WT)
        if dim == 1
            @nextract $N sizeR d->size(R,d)
            sizA1 = size(A, 1)
            @nloops $N i d->(d>1? (1:size(A,d)) : (1:1)) d->(j_d = sizeR_d==1 ? 1 : i_d) begin
                @inbounds r = (@nref $N R j)
                @inbounds m = (@nref $N means j)
                for i_1 = 1:sizA1
                    @inbounds r += f((@nref $N A i) - m) * w[i_1]
                end
                @inbounds (@nref $N R j) = r
            end
        else
            @nloops $N i A d->(if d == dim
                                   wi = w[i_d]
                                   j_d = 1
                               else
                                   j_d = i_d
                               end) @inbounds (@nref $N R j) += f((@nref $N A i) - (@nref $N means j)) * wi
        end
        return R
    end
end


# N = 1
_wsum!(R::StridedArray{T}, A::DenseArray{T,1}, w::StridedVector{T}, dim::Int, init::Bool) where {T<:BlasReal} =
    _wsum1!(R, A, w, init)

# N = 2
_wsum!(R::StridedArray{T}, A::DenseArray{T,2}, w::StridedVector{T}, dim::Int, init::Bool) where {T<:BlasReal} =
    (_wsum2_blas!(view(R,:), A, w, dim, init); R)

# N >= 3
_wsum!(R::StridedArray{T}, A::DenseArray{T,N}, w::StridedVector{T}, dim::Int, init::Bool) where {T<:BlasReal,N} =
    _wsumN!(R, A, w, dim, init)

_wsum!(R::AbstractArray, A::AbstractArray, w::AbstractVector, dim::Int, init::Bool) =
    _wsum_general!(R, identity, A, w, dim, init)

## wsum! and wsum

wsumtype(::Type{T}, ::Type{W}) where {T,W} = typeof(zero(T) * zero(W) + zero(T) * zero(W))
wsumtype(::Type{T}, ::Type{T}) where {T<:BlasReal} = T


"""
    wsum!(R, A, w, dim; init=true)

Compute the weighted sum of `A` with weights `w` over the dimension `dim` and store
the result in `R`. If `init=false`, the sum is added to `R` rather than starting
from zero.
"""
function wsum!(R::AbstractArray, A::AbstractArray{T,N}, w::AbstractVector, dim::Int; init::Bool=true) where {T,N}
    1 <= dim <= N || error("dim should be within [1, $N]")
    ndims(R) <= N || error("ndims(R) should not exceed $N")
    length(w) == size(A,dim) || throw(DimensionMismatch("Inconsistent array dimension."))
    # TODO: more careful examination of R's size
    _wsum!(R, A, w, dim, init)
end

function wsum(A::AbstractArray{T}, w::AbstractVector{W}, dim::Int) where {T<:Number,W<:Real}
    length(w) == size(A,dim) || throw(DimensionMismatch("Inconsistent array dimension."))
    _wsum!(similar(A, wsumtype(T,W), Base.reduced_indices(indices(A), dim)), A, w, dim, true)
end

# extended sum! and wsum

Base.sum!(R::AbstractArray, A::AbstractArray, w::AbstractWeights{<:Real}, dim::Int; init::Bool=true) =
    wsum!(R, A, values(w), dim; init=init)

Base.sum(A::AbstractArray{<:Number}, w::AbstractWeights{<:Real}, dim::Int) = wsum(A, values(w), dim)


###### Weighted means #####

"""
    wmean(v, w::AbstractVector)

Compute the weighted mean of an array `v` with weights `w`.
"""
function wmean(v::AbstractArray{<:Number}, w::AbstractVector)
    Base.depwarn("wmean is deprecated, use mean(v, weights(w)) instead.", :wmean)
    mean(v, weights(w))
end

"""
    mean(A::AbstractArray, w::AbstractWeights[, dim::Int])

Compute the weighted mean of array `A` with weight vector `w`
(of type `AbstractWeights`). If `dim` is provided, compute the
weighted mean along dimension `dim`.

# Examples
```julia
w = rand(n)
mean(x, weights(w))
```
"""
Base.mean(A::AbstractArray, w::AbstractWeights) = sum(A, w) / sum(w)

"""
    mean(R::AbstractArray, , A::AbstractArray, w::AbstractWeights[, dim::Int])

Compute the weighted mean of array `A` with weight vector `w`
(of type `AbstractWeights`) along dimension `dim`, and write results to `R`.
"""
Base.mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights, dim::Int) =
    scale!(Base.sum!(R, A, w, dim), inv(sum(w)))

wmeantype(::Type{T}, ::Type{W}) where {T,W} = typeof((zero(T)*zero(W) + zero(T)*zero(W)) / one(W))
wmeantype(::Type{T}, ::Type{T}) where {T<:BlasReal} = T

Base.mean(A::AbstractArray{T}, w::AbstractWeights{W}, dim::Int) where {T<:Number,W<:Real} =
    mean!(similar(A, wmeantype(T, W), Base.reduced_indices(indices(A), dim)), A, w, dim)


###### Weighted median #####
function Base.median(v::AbstractArray, w::AbstractWeights)
    throw(MethodError(median, (v, w)))
end

"""
    median(v::RealVector, w::AbstractWeights)

Compute the weighted median of `x`, using weights given by a weight vector `w`
(of type `AbstractWeights`). The weight and data vectors must have the same length.

The weighted median ``x_k`` is the element of `x` that satisfies
``\\sum_{x_i < x_k} w_i \\le \\frac{1}{2} \\sum_{j} w_j`` and
``\\sum_{x_i > x_k} w_i \\le \\frac{1}{2} \\sum_{j} w_j``.

If a weight has value zero, then its associated data point is ignored.
If none of the weights are positive, an error is thrown.
`NaN` is returned if `x` contains any `NaN` values. 
An error is raised if `w` contains any `NaN` values.
"""
function Base.median(v::RealVector, w::AbstractWeights{<:Real})
    isempty(v) && error("median of an empty array is undefined")
    if length(v) != length(w)
        error("data and weight vectors must be the same size")
    end
    @inbounds for x in w.values
        isnan(x) && error("weight vector cannot contain NaN entries")
    end
    @inbounds for x in v
        isnan(x) && return x
    end
    mask = w.values .!= 0
    if !any(mask)
        error("all weights are zero")
    end
    if all(w.values .<= 0)
        error("no positive weights found")
    end
    v = v[mask]
    wt = w[mask]
    midpoint = w.sum / 2
    maxval, maxind = findmax(wt)
    if maxval > midpoint
        v[maxind]
    else
        permute = sortperm(v)
        cumulative_weight = zero(eltype(wt))
        i = 0
        for (i, p) in enumerate(permute)
            if cumulative_weight == midpoint
                i += 1
                break
            elseif cumulative_weight > midpoint
                cumulative_weight -= wt[p]
                break
            end
            cumulative_weight += wt[p]
        end
        if cumulative_weight == midpoint
            middle(v[permute[i-2]], v[permute[i-1]])
        else
            middle(v[permute[i-1]])
        end
    end
end


"""
    wmedian(v, w)

Compute the weighted median of an array `v` with weights `w`, given as either a
vector or an `AbstractWeights` vector.
"""
wmedian(v::RealVector, w::RealVector) = median(v, weights(w))
wmedian(v::RealVector, w::AbstractWeights{<:Real}) = median(v, w)

###### Weighted quantile #####

# http://stats.stackexchange.com/questions/13169/defining-quantiles-over-a-weighted-sample
# In the non weighted version, we compute a vector of index h(N, p)
# and take interpolation between floor and ceil of this index
# Here there is a supplementary function from index to weighted index k -> Sk

"""
    quantile(v, w::AbstractWeights, p)

Compute the weighted quantiles of a vector `x` at a specified set of probability
values `p`, using weights given by a weight vector `w` (of type `AbstractWeights`).
Weights must not be negative. The weights and data vectors must have the same length.

The quantile for `p` is defined as follows. Denoting
``S_k = (k-1)w_k + (n-1) \\sum_{i<k}w_i``, define ``x_{k+1}`` the smallest element of `x`
such that ``S_{k+1}/S_{n}`` is strictly superior to `p`. The function returns
``(1-\\gamma) x_k + \\gamma x_{k+1}`` with  ``\\gamma = (pS_n- S_k)/(S_{k+1}-S_k)``.

This corresponds to  R-7, Excel, SciPy-(1,1) and Maple-6 when `w` contains only ones
(see [Wikipedia](https://en.wikipedia.org/wiki/Quantile)).
"""
function quantile(v::RealVector{V}, w::AbstractWeights{W}, p::RealVector) where {V,W<:Real}
    # checks
    isempty(v) && error("quantile of an empty array is undefined")
    isempty(p) && throw(ArgumentError("empty quantile array"))

    w.sum == 0 && error("weight vector cannot sum to zero")
    length(v) == length(w) || error("data and weight vectors must be the same size, got $(length(v)) and $(length(w))")
    for x in w.values
        isnan(x) && error("weight vector cannot contain NaN entries")
        x < 0 && error("weight vector cannot contain negative entries")
    end

    # full sort
    vw = sort!(collect(zip(v, w.values)))

    wsum = w.sum

    # prepare percentiles
    ppermute = sortperm(p)
    p = p[ppermute]
    p = bound_quantiles(p)

    # prepare out vector
    N = length(vw)
    out = Vector{typeof(zero(V)/1)}(length(p))
    fill!(out, vw[end][1])

    # start looping on quantiles
    cumulative_weight, Sk, Skold =  zero(W), zero(W), zero(W)
    vk, vkold = zero(V), zero(V)
    k = 1
    for i in 1:length(p)
        h = p[i] * (N - 1) * wsum
        if h == 0
            # happens when N or p or wsum equal zero
            out[ppermute[i]] = vw[1][1]
        else
            while Sk <= h
                # happens in particular when k == 1
                vk, wk = vw[k]
                cumulative_weight += wk
                if k >= N
                    # out was initialized with maximum v
                    return out
                end
                k += 1
                Skold, vkold = Sk, vk
                vk, wk = vw[k]
                Sk = (k - 1) * wk + (N - 1) * cumulative_weight
            end
            # in particular, Sk is different from Skold
            g = (h - Skold) / (Sk - Skold)
            out[ppermute[i]] = vkold + g * (vk - vkold)
        end
    end
    return out
end

# similarly to statistics.jl in Base
function bound_quantiles(qs::AbstractVector{T}) where T<:Real
    epsilon = 100 * eps()
    if (any(qs .< -epsilon) || any(qs .> 1+epsilon))
        throw(ArgumentError("quantiles out of [0,1] range"))
    end
    T[min(one(T), max(zero(T), q)) for q = qs]
end

quantile(v::RealVector, w::AbstractWeights{<:Real}, p::Number) = quantile(v, w, [p])[1]


"""
    wquantile(v, w, p)

Compute the `p`th quantile(s) of `v` with weights `w`, given as either a vector
or an `AbstractWeights` vector.
"""
wquantile(v::RealVector, w::AbstractWeights{<:Real}, p::RealVector) = quantile(v, w, p)
wquantile(v::RealVector, w::AbstractWeights{<:Real}, p::Number) = quantile(v, w, [p])[1]
wquantile(v::RealVector, w::RealVector, p::RealVector) = quantile(v, weights(w), p)
wquantile(v::RealVector, w::RealVector, p::Number) = quantile(v, weights(w), [p])[1]
