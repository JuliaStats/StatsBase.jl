# Descriptive Statistics


#############################
#
#   Location
#
#############################

# Geometric mean
"""
    geomean(a)

Return the geometric mean of a collection.
"""
geomean(a) = exp(mean(log, a))

# Harmonic mean
"""
    harmmean(a)

Return the harmonic mean of a collection.
"""
harmmean(a) = inv(mean(inv, a))

# Generalized mean
"""
    genmean(a, p)

Return the generalized/power mean with exponent `p` of a real-valued array,
i.e. ``\\left( \\frac{1}{n} \\sum_{i=1}^n a_i^p \\right)^{\\frac{1}{p}}``, where `n = length(a)`.
It is taken to be the geometric mean when `p == 0`.
"""
function genmean(a, p::Real)
    if p == 0
        return geomean(a)
    end

    # At least one of `x` or `p` must not be an int to avoid domain errors when `p` is a negative int.
    # We choose `x` in order to exploit exponentiation by squaring when `p` is an int.
    r = mean(a) do x
        float(x)^p
    end
    return r^inv(p)
end

# compute mode, given the range of integer values
"""
    mode(a, [r])

Return the mode (most common number) of an array, optionally
over a specified range `r`. If several modes exist, the first
one (in order of appearance) is returned.
"""
function mode(a::AbstractArray{T}, r::UnitRange{T}) where T<:Integer
    isempty(a) && throw(ArgumentError("mode is not defined for empty collections"))
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
function modes(a::AbstractArray{T}, r::UnitRange{T}) where T<:Integer
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

# compute mode over arbitrary iterable
function mode(a)
    isempty(a) && throw(ArgumentError("mode is not defined for empty collections"))
    cnts = Dict{eltype(a),Int}()
    # first element
    mc = 1
    mv, st = iterate(a)
    cnts[mv] = 1
    # find the mode along with table construction
    y = iterate(a, st)
    while y !== nothing
        x, st = y
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
        y = iterate(a, st)
    end
    return mv
end

function modes(a)
    isempty(a) && throw(ArgumentError("mode is not defined for empty collections"))
    cnts = Dict{eltype(a),Int}()
    # first element
    mc = 1
    x, st = iterate(a)
    cnts[x] = 1
    # find the mode along with table construction
    y = iterate(a, st)
    while y !== nothing
        x, st = y
        if haskey(cnts, x)
            c = (cnts[x] += 1)
            if c > mc
                mc = c
            end
        else
            cnts[x] = 1
            # in this case: c = 1, and thus c > mc won't happen
        end
        y = iterate(a, st)
    end
    # find values corresponding to maximum counts
    return [x for (x, c) in cnts if c == mc]
end


#############################
#
#   quantile and friends
#
#############################

"""
    percentile(x, p)

Return the `p`th percentile of a collection `x`, i.e. `quantile(x, p / 100)`.
"""
percentile(x, p) = quantile(x, p * 0.01)

"""
    nquantile(x, n::Integer)

Return the n-quantiles of collection `x`, i.e. the values which
partition `v` into `n` subsets of nearly equal size.

Equivalent to `quantile(x, [0:n]/n)`. For example, `nquantiles(x, 5)`
returns a vector of quantiles, respectively at `[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]`.
"""
nquantile(x, n::Integer) = quantile(x, (0:n)/n)


#############################
#
#   Dispersion
#
#############################

# span, i.e. the range minimum(x):maximum(x)
"""
    span(x)

Return the span of a collection, i.e. the range `minimum(x):maximum(x)`.
The minimum and maximum of `x` are computed in one pass using `extrema`.
"""
span(x) = ((a, b) = extrema(x); a:b)

# Variation coefficient: std / mean
"""
    variation(x, m=mean(x))

Return the coefficient of variation of collection `x`, optionally specifying
a precomputed mean `m`. The coefficient of variation is the ratio of the
standard deviation to the mean.
"""
variation(x, m) = stdm(x, m) / m
variation(x) = ((m, s) = mean_and_std(x); s/m)

# Standard error of the mean: std / sqrt(len)
# Code taken from var in the Statistics stdlib module

# faster computation of real(conj(x)*y)
realXcY(x::Real, y::Real) = x*y
realXcY(x::Complex, y::Complex) = real(x)*real(y) + imag(x)*imag(y)

"""
    sem(x)

Return the standard error of the mean of collection `x`,
i.e. `sqrt(var(x, corrected=true) / length(x))`.
"""
function sem(x)
    y = iterate(x)
    if y === nothing
        T = eltype(x)
        # Return the NaN of the type that we would get, had this collection
        # contained any elements (this is consistent with std)
        return oftype(sqrt((abs2(zero(T)) + abs2(zero(T)))/2), NaN)
    end
    count = 1
    value, state = y
    y = iterate(x, state)
    # Use Welford algorithm as seen in (among other places)
    # Knuth's TAOCP, Vol 2, page 232, 3rd edition.
    M = value / 1
    S = real(zero(M))
    while y !== nothing
        value, state = y
        y = iterate(x, state)
        count += 1
        new_M = M + (value - M) / count
        S = S + realXcY(value - M, value - new_M)
        M = new_M
    end
    var = S / (count - 1)
    return sqrt(var/count)
end

# Median absolute deviation
@irrational mad_constant 1.4826022185056018 BigFloat("1.482602218505601860547076529360423431326703202590312896536266275245674447622701")

"""
    mad(x; center=median(x), normalize=true)

Compute the median absolute deviation (MAD) of collection `x` around `center`
(by default, around the median).

If `normalize` is set to `true`, the MAD is multiplied by
`1 / quantile(Normal(), 3/4) ≈ 1.4826`, in order to obtain a consistent estimator
of the standard deviation under the assumption that the data is normally distributed.
"""
function mad(x; center=nothing, normalize::Union{Bool, Nothing}=nothing, constant=nothing)
    if normalize === nothing
        Base.depwarn("the `normalize` keyword argument will be false by default in future releases: set it explicitly to silence this deprecation", :mad)
        normalize = true
    end

    isempty(x) && throw(ArgumentError("mad is not defined for empty arrays"))
    T = eltype(x)
    # Knowing the eltype allows allocating a single array able to hold both original values
    # and differences from the center, instead of two arrays
    S = isconcretetype(T) ? promote_type(T, typeof(middle(zero(T)))) : T
    x2 = x isa AbstractArray ? LinearAlgebra.copy_oftype(x, S) : collect(S, x)
    c = center === nothing ? median!(x2) : center
    if isconcretetype(T)
        x2 .= abs.(x2 .- c)
    else
        x2 = abs.(x2 .- c)
    end
    m = median!(x2)
    if normalize isa Nothing
        Base.depwarn("the `normalize` keyword argument will be false by default in future releases: set it explicitly to silence this deprecation", :mad)
        normalize = true
    end
    if !isa(constant, Nothing)
        Base.depwarn("keyword argument `constant` is deprecated, use `normalize` instead or apply the multiplication directly", :mad)
        m * constant
    elseif normalize
        m * mad_constant
    else
        m
    end
end

"""
    StatsBase.mad!(x; center=median!(x), normalize=true)

Compute the median absolute deviation (MAD) of array `x` around `center`
(by default, around the median), overwriting `x` in the process.
`x` must be able to hold values of generated by calling `middle` on its elements
(for example an integer vector is not appropriate since `middle` can produce
non-integer values).

If `normalize` is set to `true`, the MAD is multiplied by
`1 / quantile(Normal(), 3/4) ≈ 1.4826`, in order to obtain a consistent estimator
of the standard deviation under the assumption that the data is normally distributed.
"""
function mad!(x::AbstractArray;
              center=median!(x),
              normalize::Union{Bool,Nothing}=true,
              constant=nothing)
    isempty(x) && throw(ArgumentError("mad is not defined for empty arrays"))
    x .= abs.(x .- center)
    m = median!(x)
    if normalize isa Nothing
        Base.depwarn("the `normalize` keyword argument will be false by default in future releases: set it explicitly to silence this deprecation", :mad)
        normalize = true
    end
    if !isa(constant, Nothing)
        Base.depwarn("keyword argument `constant` is deprecated, use `normalize` instead or apply the multiplication directly", :mad)
        m * constant
    elseif normalize
        m * mad_constant
    else
        m
    end
end

# Interquartile range
"""
    iqr(x)

Compute the interquartile range (IQR) of collection `x`, i.e. the 75th percentile
minus the 25th percentile.
"""
iqr(x) = (q = quantile(x, [.25, .75]); q[2] - q[1])

# Generalized variance
"""
    genvar(X)

Compute the generalized sample variance of `X`. If `X` is a vector, one-column matrix,
or other iterable, this is equivalent to the sample variance.
Otherwise if `X` is a matrix, this is equivalent to the determinant of the covariance
matrix of `X`.

!!! note
    The generalized sample variance will be 0 if the columns of the matrix of deviations
    are linearly dependent.
"""
genvar(X::AbstractMatrix) = size(X, 2) == 1 ? var(vec(X)) : det(cov(X))
genvar(itr) = var(itr)

# Total variation
"""
    totalvar(X)

Compute the total sample variance of `X`. If `X` is a vector, one-column matrix,
or other iterable, this is equivalent to the sample variance.
Otherwise if `X` is a matrix, this is equivalent to the sum of the diagonal elements
of the covariance matrix of `X`.
"""
totalvar(X::AbstractMatrix) = sum(var(X, dims=1))
totalvar(itr) = var(itr)

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

@generated function _zscore!(Z::AbstractArray{S,N}, X::AbstractArray{T,N},
                             μ::AbstractArray, σ::AbstractArray) where {S,T,N}
    quote
        # Z and X are assumed to have the same size
        # μ and σ are assumed to have the same size, that is compatible with size(X)
        siz1 = size(X, 1)
        @nextract $N ud d->size(μ, d)
        if size(μ, 1) == 1 && siz1 > 1
            @nloops $N i d->(d>1 ? (1:size(X,d)) : (1:1)) d->(j_d = ud_d ==1 ? 1 : i_d) begin
                v = (@nref $N μ j)
                c = inv(@nref $N σ j)
                for i_1 = 1:siz1
                    (@nref $N Z i) = ((@nref $N X i) - v) * c
                end
            end
        else
            @nloops $N i X d->(j_d = ud_d ==1 ? 1 : i_d) begin
                (@nref $N Z i) = ((@nref $N X i) - (@nref $N μ j)) / (@nref $N σ j)
            end
        end
        return Z
    end
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
function zscore!(Z::AbstractArray{ZT}, X::AbstractArray{T}, μ::Real, σ::Real) where {ZT<:AbstractFloat,T<:Real}
    size(Z) == size(X) || throw(DimensionMismatch("Z and X must have the same size."))
    _zscore!(Z, X, μ, σ)
end

function zscore!(Z::AbstractArray{<:AbstractFloat}, X::AbstractArray{<:Real},
                 μ::AbstractArray{<:Real}, σ::AbstractArray{<:Real})
    size(Z) == size(X) || throw(DimensionMismatch("Z and X must have the same size."))
    _zscore_chksize(X, μ, σ)
    _zscore!(Z, X, μ, σ)
end

zscore!(X::AbstractArray{<:AbstractFloat}, μ::Real, σ::Real) = _zscore!(X, X, μ, σ)

zscore!(X::AbstractArray{<:AbstractFloat}, μ::AbstractArray{<:Real}, σ::AbstractArray{<:Real}) =
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
function zscore(X::AbstractArray{T}, μ::Real, σ::Real) where T<:Real
    ZT = typeof((zero(T) - zero(μ)) / one(σ))
    _zscore!(Array{ZT}(undef, size(X)), X, μ, σ)
end

function zscore(X::AbstractArray{T}, μ::AbstractArray{U}, σ::AbstractArray{S}) where {T<:Real,U<:Real,S<:Real}
    _zscore_chksize(X, μ, σ)
    ZT = typeof((zero(T) - zero(U)) / one(S))
    _zscore!(Array{ZT}(undef, size(X)), X, μ, σ)
end

zscore(X::AbstractArray{<:Real}) = ((μ, σ) = mean_and_std(X); zscore(X, μ, σ))
zscore(X::AbstractArray{<:Real}, dim::Int) = ((μ, σ) = mean_and_std(X, dim); zscore(X, μ, σ))



#############################
#
#   entropy and friends
#
#############################

"""
    entropy(p, [b])

Compute the entropy of a collection of probabilities `p`,
optionally specifying a real number `b` such that the entropy is scaled by `1/log(b)`.
Elements with probability 0 or 1 add 0 to the entropy.
"""
entropy(p) = -sum(pᵢ -> iszero(pᵢ) ? zero(pᵢ) : pᵢ * log(pᵢ), p)

entropy(p, b::Real) = entropy(p) / log(b)

"""
    renyientropy(p, α)

Compute the Rényi (generalized) entropy of order `α` of an array `p`.
"""
function renyientropy(p::AbstractArray{T}, α::Real) where T<:Real
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
function crossentropy(p::AbstractArray{T}, q::AbstractArray{T}) where T<:Real
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

crossentropy(p::AbstractArray{T}, q::AbstractArray{T}, b::Real) where {T<:Real} =
    crossentropy(p,q) / log(b)


"""
    kldivergence(p, q, [b])

Compute the Kullback-Leibler divergence from `q` to `p`,
also called the relative entropy of `p` with respect to `q`,
that is the sum `pᵢ * log(pᵢ / qᵢ)`. Optionally a real number `b`
can be specified such that the divergence is scaled by `1/log(b)`.
"""
function kldivergence(p::AbstractArray{T}, q::AbstractArray{T}) where T<:Real
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

kldivergence(p::AbstractArray{T}, q::AbstractArray{T}, b::Real) where {T<:Real} =
    kldivergence(p,q) / log(b)

#############################
#
#   summary
#
#############################

struct SummaryStats{T<:Union{AbstractFloat,Missing}}
    mean::T
    min::T
    q25::T
    median::T
    q75::T
    max::T
    nobs::Int
    nmiss::Int
end


"""
    summarystats(a)

Compute summary statistics for a real-valued array `a`. Returns a
`SummaryStats` object containing the mean, minimum, 25th percentile,
median, 75th percentile, and maxmimum.
"""
function summarystats(a::AbstractArray{T}) where T<:Union{Real,Missing}
    # `mean` doesn't fail on empty input but rather returns `NaN`, so we can use the
    # return type to populate the `SummaryStats` structure.
    s = T >: Missing ? collect(skipmissing(a)) : a
    m = mean(s)
    R = typeof(m)
    n = length(a)
    ns = length(s)
    qs = if m == 0 || n == 0
        R[NaN, NaN, NaN, NaN, NaN]
    elseif T >: Missing
        quantile!(s, [0.00, 0.25, 0.50, 0.75, 1.00])
    else
        quantile(s, [0.00, 0.25, 0.50, 0.75, 1.00])
    end
    SummaryStats{R}(m, qs..., n, n - ns)
end

function Base.show(io::IO, ss::SummaryStats)
    println(io, "Summary Stats:")
    @printf(io, "Length:         %i\n", ss.nobs)
    ss.nobs > 0 || return
    @printf(io, "Missing Count:  %i\n", ss.nmiss)
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
describe(a::AbstractArray) = describe(stdout, a)
function describe(io::IO, a::AbstractArray{T}) where T<:Union{Real,Missing}
    show(io, summarystats(a))
    println(io, "Type:           $(string(eltype(a)))")
end
function describe(io::IO, a::AbstractArray)
    println(io, "Summary Stats:")
    println(io, "Length:         $(length(a))")
    println(io, "Type:           $(string(eltype(a)))")
    println(io, "Number Unique:  $(length(unique(a)))")
    return
end
