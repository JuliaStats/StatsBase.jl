###### Weight vector #####
abstract type AbstractWeights{S<:Real, T<:Real, V<:AbstractVector{T}} <: AbstractVector{T} end

"""
    @weights name

Generates a new generic weight type with specified `name`, which subtypes `AbstractWeights`
and stores the `values` (`V<:RealVector`) and `sum` (`S<:Real`).
"""
macro weights(name)
    return quote
        mutable struct $name{S<:Real, T<:Real, V<:AbstractVector{T}} <: AbstractWeights{S, T, V}
            values::V
            sum::S
        end
        $(esc(name))(vs) = $(esc(name))(vs, sum(vs))
    end
end

eltype(wv::AbstractWeights) = eltype(wv.values)
length(wv::AbstractWeights) = length(wv.values)
values(wv::AbstractWeights) = wv.values
sum(wv::AbstractWeights) = wv.sum
isempty(wv::AbstractWeights) = isempty(wv.values)

Base.getindex(wv::AbstractWeights, i) = getindex(wv.values, i)
Base.size(wv::AbstractWeights) = size(wv.values)

@propagate_inbounds function Base.setindex!(wv::AbstractWeights, v::Real, i::Int)
    s = v - wv[i]
    wv.values[i] = v
    wv.sum += s
    v
end

"""
    varcorrection(n::Integer, corrected=false)

Compute a bias correction factor for calculating `var`, `std` and `cov` with
`n` observations. Returns ``\\frac{1}{n - 1}`` when `corrected=true`
(i.e. [Bessel's correction](https://en.wikipedia.org/wiki/Bessel's_correction)),
otherwise returns ``\\frac{1}{n}`` (i.e. no correction).
"""
@inline varcorrection(n::Integer, corrected::Bool=false) = 1 / (n - Int(corrected))

@weights Weights

@doc """
    Weights(vs, wsum=sum(vs))

Construct a `Weights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

The `Weights` type describes a generic weights vector which does not support
all operations possible for [`FrequencyWeights`](@ref), [`AnalyticWeights`](@ref)
and [`ProbabilityWeights`](@ref).
""" Weights

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

@doc """
    AnalyticWeights(vs, wsum=sum(vs))

Construct an `AnalyticWeights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

Analytic weights describe a non-random relative importance (usually between 0 and 1)
for each observation. These weights may also be referred to as reliability weights,
precision weights or inverse variance weights. These are typically used when the observations
being weighted are aggregate values (e.g., averages) with differing variances.
""" AnalyticWeights

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

@doc """
    FrequencyWeights(vs, wsum=sum(vs))

Construct a `FrequencyWeights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

Frequency weights describe the number of times (or frequency) each observation
was observed. These weights may also be referred to as case weights or repeat weights.
""" FrequencyWeights

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

@doc """
    ProbabilityWeights(vs, wsum=sum(vs))

Construct a `ProbabilityWeights` vector with weight values `vs`.
A precomputed sum may be provided as `wsum`.

Probability weights represent the inverse of the sampling probability for each observation,
providing a correction mechanism for under- or over-sampling certain population groups.
These weights may also be referred to as sampling weights.
""" ProbabilityWeights

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

##### Equality tests #####

for w in (AnalyticWeights, FrequencyWeights, ProbabilityWeights, Weights)
    @eval begin
        Base.isequal(x::$w, y::$w) = isequal(x.sum, y.sum) && isequal(x.values, y.values)
        Base.:(==)(x::$w, y::$w)   = (x.sum == y.sum) && (x.values == y.values)
    end
end

Base.isequal(x::AbstractWeights, y::AbstractWeights) = false
Base.:(==)(x::AbstractWeights, y::AbstractWeights)   = false

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
Base.sum(v::SparseArrays.SparseMatrixCSC, w::AbstractWeights) = wsum(v, values(w))
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
            @nloops $N i d->(d>1 ? (1:size(A,d)) : (1:1)) d->(j_d = sizeR_d==1 ? 1 : i_d) begin
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
            @nloops $N i d->(d>1 ? (1:size(A,d)) : (1:1)) d->(j_d = sizeR_d==1 ? 1 : i_d) begin
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
    _wsum!(similar(A, wsumtype(T,W), Base.reduced_indices(axes(A), dim)), A, w, dim, true)
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
    mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights[; dims=nothing])

Compute the weighted mean of array `A` with weight vector `w`
(of type `AbstractWeights`) along dimension `dims`, and write results to `R`.
"""
mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights;
      dims::Union{Nothing,Int}=nothing) = _mean!(R, A, w, dims)
_mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights, dims::Nothing) = throw(ArgumentError("dims argument must be provided"))
_mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights, dims::Int) =
    rmul!(Base.sum!(R, A, w, dims), inv(sum(w)))

wmeantype(::Type{T}, ::Type{W}) where {T,W} = typeof((zero(T)*zero(W) + zero(T)*zero(W)) / one(W))
wmeantype(::Type{T}, ::Type{T}) where {T<:BlasReal} = T

"""
    mean(A::AbstractArray, w::AbstractWeights[, dims::Int])

Compute the weighted mean of array `A` with weight vector `w`
(of type `AbstractWeights`). If `dim` is provided, compute the
weighted mean along dimension `dims`.

# Examples
```julia
w = rand(n)
mean(x, weights(w))
```
"""
mean(A::AbstractArray, w::AbstractWeights; dims::Union{Nothing,Int}=nothing) =
    _mean(A, w, dims)
_mean(A::AbstractArray, w::AbstractWeights, dims::Nothing) =
    sum(A, w) / sum(w)
_mean(A::AbstractArray{T}, w::AbstractWeights{W}, dims::Int) where {T,W} =
    _mean!(similar(A, wmeantype(T, W), Base.reduced_indices(axes(A), dims)), A, w, dims)


###### Weighted quantile #####
"""
    quantile(v, w::AbstractWeights, p)

Compute the weighted quantiles of a vector `v` at a specified set of probability
values `p`, using weights given by a weight vector `w` (of type `AbstractWeights`).
Weights must not be negative. The weights and data vectors must have the same length.
`NaN` is returned if `x` contains any `NaN` values. An error is raised if `w` contains 
any `NaN` values.

With [`FrequencyWeights`](@ref), the function returns the same result as
`quantile` for a vector with repeated values. Weights must be integers.

With non `FrequencyWeights`,  denote ``N`` the length of the vector, ``w`` the vector of weights,
``h = p (\\sum_{i<= N} w_i - w_1) + w_1`` the cumulative weight corresponding to the
probability ``p`` and ``S_k = \\sum_{i<=k} w_i`` the cumulative weight for each
observation, define ``v_{k+1}`` the smallest element of `v` such that ``S_{k+1}``
is strictly superior to ``h``. The weighted ``p`` quantile is given by ``v_k + \\gamma (v_{k+1} - v_k)``
with  ``\\gamma = (h - S_k)/(S_{k+1} - S_k)``. In particular, when all weights are equal,
the function returns the same result as the unweighted `quantile`.
"""
function quantile(v::RealVector{V}, w::AbstractWeights{W}, p::RealVector) where {V,W<:Real}
    # checks
    isempty(v) && throw(ArgumentError("quantile of an empty array is undefined"))
    isempty(p) && throw(ArgumentError("empty quantile array"))
    all(x -> 0 <= x <= 1, p) || throw(ArgumentError("input probability out of [0,1] range"))

    w.sum == 0 && throw(ArgumentError("weight vector cannot sum to zero"))
    length(v) == length(w) || throw(ArgumentError("data and weight vectors must be the same size," * 
        "got $(length(v)) and $(length(w))"))
    for x in w.values
        isnan(x) && throw(ArgumentError("weight vector cannot contain NaN entries"))
        x < 0 && throw(ArgumentError("weight vector cannot contain negative entries"))
    end

    isa(w, FrequencyWeights) && !(eltype(w) <: Integer) && any(!isinteger, w) && 
        throw(ArgumentError("The values of the vector of `FrequencyWeights` must be numerically" * 
                            "equal to integers. Use `ProbabilityWeights` or `AnalyticWeights` instead."))

    # remove zeros weights and sort
    wsum = sum(w)
    nz = .!iszero.(w)
    vw = sort!(collect(zip(view(v, nz), view(w, nz))))
    N = length(vw)

    # prepare percentiles
    ppermute = sortperm(p)
    p = p[ppermute]

    # prepare out vector
    out = Vector{typeof(zero(V)/1)}(undef, length(p))
    fill!(out, vw[end][1])

    @inbounds for x in v
        isnan(x) && return fill!(out, x)
    end

    # loop on quantiles
    Sk, Skold = zero(W), zero(W)
    vk, vkold = zero(V), zero(V)
    k = 0

    w1 = vw[1][2]
    for i in 1:length(p)
        if isa(w, FrequencyWeights)
            h = p[i] * (wsum - 1) + 1
        else
            h = p[i] * (wsum - w1) + w1
        end
        while Sk <= h
            k += 1
            if k > N
               # out was initialized with maximum v
               return out
            end
            Skold, vkold = Sk, vk
            vk, wk = vw[k]
            Sk += wk
        end
        if isa(w, FrequencyWeights)
            out[ppermute[i]] = vkold + min(h - Skold, 1) * (vk - vkold)
        else
            out[ppermute[i]] = vkold + (h - Skold) / (Sk - Skold) * (vk - vkold)
        end
    end
    return out
end

quantile(v::RealVector, w::AbstractWeights{<:Real}, p::Number) = quantile(v, w, [p])[1]



###### Weighted median #####
"""
    median(v::RealVector, w::AbstractWeights)

Compute the weighted median of `v` with weights `w`
(of type `AbstractWeights`). See the documentation for [`quantile`](@ref) for more details.
"""
median(v::RealVector, w::AbstractWeights{<:Real}) = quantile(v, w, 0.5)
