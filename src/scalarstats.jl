# Descriptive Statistics


#############################
#
#   Location
#
#############################

# Geometric mean
"""
    geomean(a)

Return the geometric mean of a real-valued array.
"""
function geomean(a::RealArray)
    s = 0.0
    n = length(a)
    for i = 1 : n
        @inbounds s += log(a[i])
    end
    return exp(s / n)
end

# Harmonic mean
"""
    harmmean(a)

Return the harmonic mean of a real-valued array.
"""
function harmmean(a::RealArray)
    s = 0.0
    n = length(a)
    for i in 1 : n
        @inbounds s += inv(a[i])
    end
    return n / s
end

# Generalized mean
"""
    genmean(a, p)

Return the generalized/power mean with exponent `p` of a real-valued array,
i.e. ``\\left( \\frac{1}{n} \\sum_{i=1}^n a_i^p \\right)^{\\frac{1}{p}}``, where `n = length(a)`.
It is taken to be the geometric mean when `p == 0`.
"""
function genmean(a::RealArray, p::Real)
    if p == 0
        return geomean(a)
    end
    s = 0.0
    n = length(a)
    for x in a
        #= At least one of `x` or `p` must not be an int to avoid domain errors when `p` is a negative int.
        We choose `x` in order to exploit exponentiation by squaring when `p` is an int. =#
        @inbounds s += convert(Float64, x)^p
    end
    return (s/n)^(1/p)
end

# compute mode, given the range of integer values
"""
    mode(a, [r])

Return the mode (most common number) of an array, optionally
over a specified range `r`. If several modes exist, the first
one (in order of appearance) is returned.
"""
function mode{T<:Integer}(a::AbstractArray{T}, r::UnitRange{T})
    isempty(a) && error("mode: input array cannot be empty.")
    len = length(a)
    r0 = r[1]
    r1 = r[end]
    cnts = zeros(Int, length(r))
    mc = 0    # maximum count
    mv = r0   # a value corresponding to maximum count
    for i = 1:len
        @inbounds x = a[i]
        if r0 <= x <= r1
            @inbounds c = (cnts[x - r0 + 1] += 1)
            if c > mc
                mc = c
                mv = x
            end
        end
    end
    return mv
end

"""
    modes(a, [r])::Vector

Return all modes (most common numbers) of an array, optionally over a
specified range `r`.
"""
function modes{T<:Integer}(a::AbstractArray{T}, r::UnitRange{T})
    r0 = r[1]
    r1 = r[end]
    n = length(r)
    cnts = zeros(Int, n)
    # find the maximum count
    mc = 0
    for i = 1:length(a)
        @inbounds x = a[i]
        if r0 <= x <= r1
            @inbounds c = (cnts[x - r0 + 1] += 1)
            if c > mc
                mc = c
            end
        end
    end
    # find all values corresponding to maximum count
    ms = T[]
    for i = 1:n
        @inbounds if cnts[i] == mc
            push!(ms, r[i])
        end
    end
    return ms
end

# compute mode over arbitrary array
function mode{T}(a::AbstractArray{T})
    isempty(a) && error("mode: input array cannot be empty.")
    cnts = Dict{T,Int}()
    # first element
    mc = 1
    mv = a[1]
    cnts[mv] = 1
    # find the mode along with table construction
    for i = 2 : length(a)
        @inbounds x = a[i]
        if haskey(cnts, x)
            c = (cnts[x] += 1)
            if c > mc
                mc = c
                mv = x
            end
        else
            cnts[x] = 1
            # in this case: c = 1, and thus c > mc won't happen
        end
    end
    return mv
end

function modes{T}(a::AbstractArray{T})
    isempty(a) && error("modes: input array cannot be empty.")
    cnts = Dict{T,Int}()
    # first element
    mc = 1
    cnts[a[1]] = 1
    # find the mode along with table construction
    for i = 2 : length(a)
        @inbounds x = a[i]
        if haskey(cnts, x)
            c = (cnts[x] += 1)
            if c > mc
                mc = c
            end
        else
            cnts[x] = 1
            # in this case: c = 1, and thus c > mc won't happen
        end
    end
    # find values corresponding to maximum counts
    ms = T[]
    for (x, c) in cnts
        if c == mc
            push!(ms, x)
        end
    end
    return ms
end


#############################
#
#   quantile and friends
#
#############################

"""
    percentile(v, p)

Return the `p`th percentile of a real-valued array `v`, i.e. `quantile(x, p / 100)`.
"""
percentile{T<:Real}(v::AbstractArray{T}, p) = quantile(v, p * 0.01)

quantile{T<:Real}(v::AbstractArray{T}) = quantile(v, [.0, .25, .5, .75, 1.0])

"""
    nquantile(v, n)

Return the n-quantiles of a real-valued array, i.e. the values which
partition `v` into `n` subsets of nearly equal size.

Equivalent to `quantile(v, [0:n]/n)`. For example, `nquantiles(x, 5)`
returns a vector of quantiles, respectively at `[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]`.
"""
nquantile{T<:Real}(v::AbstractArray{T}, n::Integer) = quantile(v, (0:n)/n)


#############################
#
#   Dispersion
#
#############################

# span, i.e. the range minimum(x):maximum(x)
"""
    span(x)

Return the span of an integer array, i.e. the range `minimum(x):maximum(x)`.
The minimum and maximum of `x` are computed in one-pass using `extrema`.
"""
span{T<:Integer}(x::AbstractArray{T}) = ((a, b) = extrema(x); a:b)

# Variation coefficient: std / mean
"""
    variation(x, m=mean(x))

Return the coefficient of variation of an array `x`, optionally specifying
a precomputed mean `m`. The coefficient of variation is the ratio of the
standard deviation to the mean.
"""
variation{T<:Real}(x::AbstractArray{T}, m::Real) = stdm(x, m) / m
variation{T<:Real}(x::AbstractArray{T}) = variation(x, mean(x))

# Standard error of the mean: std(a) / sqrt(len)
"""
    sem(a)

Return the standard error of the mean of `a`, i.e. `sqrt(var(a) / length(a))`.
"""
sem{T<:Real}(a::AbstractArray{T}) = sqrt(var(a) / length(a))

# Median absolute deviation
"""
    mad(v)

Compute the median absolute deviation of `v`.
"""
function mad{T<:Real}(v::AbstractArray{T})
    isempty(v) && throw(ArgumentError("mad is not defined for empty arrays"))

    S = promote_type(T, typeof(middle(first(v))))

    mad!(LinAlg.copy_oftype(v, S))
end


"""
    StatsBase.mad!(v, center=median!(v); constant=k)

Compute the maximum absolute deviation (MAD) of `v` about a precomputed center
`center`, overwriting `v` in the process. Using the MAD as a consistent estimator
of the standard deviation requires a scaling factor that depends on the underlying
distribution. For normally distributed data, `k` is chosen as
`1 / quantile(Normal(), 3/4) ≈ 1.4826`, which is used as the default here.
"""
function mad!{T<:Real}(v::AbstractArray{T}, center::Real=median!(v);
                       constant::Real = 1 / (-sqrt(2 * one(T)) * erfcinv(3 * one(T) / 2)))
    for i in 1:length(v)
        @inbounds v[i] = abs(v[i]-center)
    end
    constant * median!(v)
end

# Interquartile range
"""
    iqr(v)

Compute the interquartile range (IQR) of an array, i.e. the 75th percentile
minus the 25th percentile.
"""
iqr{T<:Real}(v::AbstractArray{T}) = (q = quantile(v, [.25, .75]); q[2] - q[1])


#############################
#
#   Z-scores
#
#############################

function _zscore!(Z::AbstractArray, X::AbstractArray, μ::Real, σ::Real)
    # Z and X are assumed to have the same size
    iσ = inv(σ)
    if μ == zero(μ)
        for i = 1 : length(X)
            @inbounds Z[i] = X[i] * iσ
        end
    else
        for i = 1 : length(X)
            @inbounds Z[i] = (X[i] - μ) * iσ
        end
    end
    return Z
end

@ngenerate N typeof(Z) function _zscore!{S,T,N}(Z::AbstractArray{S,N}, X::AbstractArray{T,N}, μ::AbstractArray, σ::AbstractArray)
    # Z and X are assumed to have the same size
    # μ and σ are assumed to have the same size, that is compatible with size(X)
    siz1 = size(X, 1)
    @nextract N ud d->size(μ, d)
    if size(μ, 1) == 1 && siz1 > 1
        @nloops N i d->(d>1 ? (1:size(X,d)) : (1:1)) d->(j_d = ud_d ==1 ? 1 : i_d) begin
            v = (@nref N μ j)
            c = inv(@nref N σ j)
            for i_1 = 1:siz1
                (@nref N Z i) = ((@nref N X i) - v) * c
            end
        end
    else
        @nloops N i X d->(j_d = ud_d ==1 ? 1 : i_d) begin
            (@nref N Z i) = ((@nref N X i) - (@nref N μ j)) / (@nref N σ j)
        end
    end
    return Z
end

function _zscore_chksize(X::AbstractArray, μ::AbstractArray, σ::AbstractArray)
    size(μ) == size(σ) || throw(DimensionMismatch("μ and σ should have the same size."))
    for i=1:ndims(X)
        dμ_i = size(μ,i)
        (dμ_i == 1 || dμ_i == size(X,i)) || throw(DimensionMismatch("X and μ have incompatible sizes."))
    end
end


"""
    zscore!([Z], X, μ, σ)

Compute the z-scores of an array `X` with mean `μ` and standard deviation `σ`.
z-scores are the signed number of standard deviations above the mean that an
observation lies, i.e. ``(x - μ) / σ``.

If a destination array `Z` is provided, the scores are stored
in `Z` and it must have the same shape as `X`. Otherwise `X` is overwritten.
"""
function zscore!{ZT<:AbstractFloat,T<:Real}(Z::AbstractArray{ZT}, X::AbstractArray{T}, μ::Real, σ::Real)
    size(Z) == size(X) || throw(DimensionMismatch("Z and X must have the same size."))
    _zscore!(Z, X, μ, σ)
end

function zscore!{ZT<:AbstractFloat,T<:Real,U<:Real,S<:Real}(Z::AbstractArray{ZT}, X::AbstractArray{T},
                                                            μ::AbstractArray{U}, σ::AbstractArray{S})
    size(Z) == size(X) || throw(DimensionMismatch("Z and X must have the same size."))
    _zscore_chksize(X, μ, σ)
    _zscore!(Z, X, μ, σ)
end

zscore!{T<:AbstractFloat}(X::AbstractArray{T}, μ::Real, σ::Real) = _zscore!(X, X, μ, σ)

zscore!{T<:AbstractFloat,U<:Real,S<:Real}(X::AbstractArray{T}, μ::AbstractArray{U}, σ::AbstractArray{S}) =
    (_zscore_chksize(X, μ, σ); _zscore!(X, X, μ, σ))


"""
    zscore(X, [μ, σ])

Compute the z-scores of `X`, optionally specifying a precomputed mean `μ` and
standard deviation `σ`. z-scores are the signed number of standard deviations
above the mean that an observation lies, i.e. ``(x - μ) / σ``.

`μ` and `σ` should be both scalars or both arrays. The computation is broadcasting.
In particular, when `μ` and `σ` are arrays, they should have the same size, and
`size(μ, i) == 1  || size(μ, i) == size(X, i)` for each dimension.
"""
function zscore{T<:Real}(X::AbstractArray{T}, μ::Real, σ::Real)
    ZT = typeof((zero(T) - zero(μ)) / one(σ))
    _zscore!(Array{ZT}(size(X)), X, μ, σ)
end

function zscore{T<:Real,U<:Real,S<:Real}(X::AbstractArray{T}, μ::AbstractArray{U}, σ::AbstractArray{S})
    _zscore_chksize(X, μ, σ)
    ZT = typeof((zero(T) - zero(U)) / one(S))
    _zscore!(Array{ZT}(size(X)), X, μ, σ)
end

zscore{T<:Real}(X::AbstractArray{T}) = ((μ, σ) = mean_and_std(X); zscore(X, μ, σ))
zscore{T<:Real}(X::AbstractArray{T}, dim::Int) = ((μ, σ) = mean_and_std(X, dim); zscore(X, μ, σ))



#############################
#
#   entropy and friends
#
#############################

"""
    entropy(p, [b])

Compute the entropy of an array `p`, optionally specifying a real number
`b` such that the entropy is scaled by `1/log(b)`.
"""
function entropy{T<:Real}(p::AbstractArray{T})
    s = zero(T)
    z = zero(T)
    for i = 1:length(p)
        @inbounds pi = p[i]
        if pi > z
            s += pi * log(pi)
        end
    end
    return -s
end

function entropy{T<:Real}(p::AbstractArray{T}, b::Real)
    return entropy(p) / log(b)
end

"""
    renyientropy(p, α)

Compute the Rényi (generalized) entropy of order `α` of an array `p`.
"""
function renyientropy{T<:Real, U<:Real}(p::AbstractArray{T}, α::U)
    α < 0 && throw(ArgumentError("Order of Rényi entropy not legal, $(α) < 0."))

    s = zero(T)
    z = zero(T)
    scale = sum(p)

    if α ≈ 0
        for i = 1:length(p)
            @inbounds pi = p[i]
            if pi > z
                s += 1
            end
        end
        s = log(s / scale)
    elseif α ≈ 1
        for i = 1:length(p)
            @inbounds pi = p[i]
            if pi > z
                s -= pi * log(pi)
            end
        end
        s = s / scale
    elseif (isinf(α))
        s = -log(maximum(p))
    else # a normal Rényi entropy
        for i = 1:length(p)
            @inbounds pi = p[i]
            if pi > z
                s += pi ^ α
            end
        end
        s = log(s / scale) / (1 - α)
    end
    return s
end

"""
    crossentropy(p, q, [b])

Compute the cross entropy between `p` and `q`, optionally specifying a real
number `b` such that the result is scaled by `1/log(b)`.
"""
function crossentropy{T<:Real}(p::AbstractArray{T}, q::AbstractArray{T})
    length(p) == length(q) || throw(DimensionMismatch("Inconsistent array length."))
    s = 0.
    z = zero(T)
    for i = 1:length(p)
        @inbounds pi = p[i]
        @inbounds qi = q[i]
        if pi > z
            s += pi * log(qi)
        end
    end
    return -s
end

function crossentropy{T<:Real}(p::AbstractArray{T}, q::AbstractArray{T}, b::Real)
    return crossentropy(p,q) / log(b)
end


"""
    kldivergence(p, q, [b])

Compute the Kullback-Leibler divergence of `q` from `p`, optionally specifying
a real number `b` such that the divergence is scaled by `1/log(b)`.
"""
function kldivergence{T<:Real}(p::AbstractArray{T}, q::AbstractArray{T})
    length(p) == length(q) || throw(DimensionMismatch("Inconsistent array length."))
    s = 0.
    z = zero(T)
    for i = 1:length(p)
        @inbounds pi = p[i]
        @inbounds qi = q[i]
        if pi > z
            s += pi * log(pi / qi)
        end
    end
    return s
end

function kldivergence{T<:Real}(p::AbstractArray{T}, q::AbstractArray{T}, b::Real)
    return kldivergence(p,q) / log(b)
end

#############################
#
#   summary
#
#############################

immutable SummaryStats{T<:AbstractFloat}
    mean::T
    min::T
    q25::T
    median::T
    q75::T
    max::T
end


"""
    summarystats(a)

Compute summary statistics for a real-valued array `a`. Returns a
`SummaryStats` object containing the mean, minimum, 25th percentile,
median, 75th percentile, and maxmimum.
"""
function summarystats{T<:Real}(a::AbstractArray{T})
    m = mean(a)
    qs = quantile(a, [0.00, 0.25, 0.50, 0.75, 1.00])
    R = typeof(convert(AbstractFloat, zero(T)))
    SummaryStats{R}(
        convert(R, m),
        convert(R, qs[1]),
        convert(R, qs[2]),
        convert(R, qs[3]),
        convert(R, qs[4]),
        convert(R, qs[5]))
end

function Base.show(io::IO, ss::SummaryStats)
    println(io, "Summary Stats:")
    @printf(io, "Mean:           %.6f\n", ss.mean)
    @printf(io, "Minimum:        %.6f\n", ss.min)
    @printf(io, "1st Quartile:   %.6f\n", ss.q25)
    @printf(io, "Median:         %.6f\n", ss.median)
    @printf(io, "3rd Quartile:   %.6f\n", ss.q75)
    @printf(io, "Maximum:        %.6f\n", ss.max)
end


"""
    describe(a)

Pretty-print the summary statistics provided by [`summarystats`](@ref):
the mean, minimum, 25th percentile, median, 75th percentile, and
maximum.
"""
describe(a::AbstractArray) = describe(STDOUT, a)
function describe{T<:Real}(io::IO, a::AbstractArray{T})
    show(io, summarystats(a))
    println(io, "Length:         $(length(a))")
    println(io, "Type:           $(string(eltype(a)))")
end
function describe(io::IO, a::AbstractArray)
    println(io, "Summary Stats:")
    println(io, "Length:         $(length(a))")
    println(io, "Type:           $(string(eltype(a)))")
    println(io, "Number Unique:  $(length(unique(a)))")
    return
end
