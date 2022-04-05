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
    mode(a::AbstractArray, wv::AbstractWeights)

Return the mode (most common number) of an array, optionally
over a specified range `r` or weighted via a vector `wv`.
If several modes exist, the first one (in order of appearance) is returned.
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
    mode(a::AbstractArray, wv::AbstractWeights)::Vector

Return all modes (most common numbers) of an array, optionally over a
specified range `r` or weighted via vector `wv`.
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

# Weighted mode of arbitrary vectors of values
function mode(a::AbstractVector, wv::AbstractWeights{T}) where T <: Real
    isempty(a) && throw(ArgumentError("mode is not defined for empty collections"))
    length(a) == length(wv) ||
        throw(ArgumentError("data and weight vectors must be the same size, got $(length(a)) and $(length(wv))"))

    # Iterate through the data
    mv = first(a)
    mw = first(wv)
    weights = Dict{eltype(a), T}()
    for (x, w) in zip(a, wv)
        _w = get!(weights, x, zero(T)) + w
        if _w > mw
            mv = x
            mw = _w
        end
        weights[x] = _w
    end

    return mv
end

function modes(a::AbstractVector, wv::AbstractWeights{T}) where T <: Real
    isempty(a) && throw(ArgumentError("mode is not defined for empty collections"))
    length(a) == length(wv) ||
        throw(ArgumentError("data and weight vectors must be the same size, got $(length(a)) and $(length(wv))"))

    # Iterate through the data
    mw = first(wv)
    weights = Dict{eltype(a), T}()
    for (x, w) in zip(a, wv)
        _w = get!(weights, x, zero(T)) + w
        if _w > mw
            mw = _w
        end
        weights[x] = _w
    end

    # find values corresponding to maximum counts
    return [x for (x, w) in weights if w == mw]
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

"""
    quantilerank(itr, value; method=:inc)

Compute the quantile position in the [0, 1] interval of `value` relative to collection `itr`.

Different definitions can be chosen via the `method` keyword argument.
Let `count_less` be the number of elements of `itr` that are less than `value`, 
`count_equal` the number of elements of `itr` that are equal to `value`, `n` the length of `itr`, 
`greatest_smaller` the highest value below `value` and `smallest_greater` the lowest value above `value`. 
Then `method` supports the following definitions:

- `:inc` (default): Return a value in the range 0 to 1 inclusive. 
Return `count_less / (n - 1)` if `value ∈ itr`, otherwise apply interpolation based on 
definition 7 of quantile in Hyndman and Fan (1996)
(equivalent to Excel `PERCENTRANK` and `PERCENTRANK.INC`).
This definition corresponds to the lower semi-continuous inverse of
[`quantile`](@ref) with its default parameters.

- `:exc`: Return a value in the range 0 to 1 exclusive.
Return `(count_less + 1) / (n + 1)` if `value ∈ itr` otherwise apply interpolation 
based on definition 6 of quantile in Hyndman and Fan (1996)
(equivalent to Excel `PERCENTRANK.EXC`).

- `:compete`: Return `count_less / (n - 1)` if `value ∈ itr`, otherwise 
return `(count_less - 1) / (n - 1)`, without interpolation
(equivalent to MariaDB `PERCENT_RANK`, dplyr `percent_rank`).

- `:tied`: Return `(count_less + count_equal/2) / n`, without interpolation.
Based on the definition in Roscoe, J. T. (1975)
(equivalent to `"mean"` kind of SciPy `percentileofscore`).

- `:strict`: Return `count_less / n`, without interpolation
(equivalent to `"strict"` kind of SciPy `percentileofscore`).

- `:weak`: Return `(count_less + count_equal) / n`, without interpolation
(equivalent to `"weak"` kind of SciPy `percentileofscore`).

!!! note
    An `ArgumentError` is thrown if `itr` contains `NaN` or `missing` values
    or if `itr` contains fewer than two elements.

# References
Roscoe, J. T. (1975). [Fundamental Research Statistics for the Behavioral Sciences]
(http://www.bryanburnham.net/wp-content/uploads/2014/07/Fundamental-Statistics-for-the-Behavioral-Sciences-v2.0.pdf#page=57)",
2nd ed., New York : Holt, Rinehart and Winston.

Hyndman, R.J and Fan, Y. (1996) "[Sample Quantiles in Statistical Packages]
(https://www.amherst.edu/media/view/129116/original/Sample+Quantiles.pdf)",
*The American Statistician*, Vol. 50, No. 4, pp. 361-365.

# Examples
```julia
julia> using StatsBase

julia> v1 = [1, 1, 1, 2, 3, 4, 8, 11, 12, 13];

julia> v2 = [1, 2, 3, 5, 6, missing, 8];

julia> v3 = [1, 2, 3, 4, 4, 5, 6, 7, 8, 9];

julia> quantilerank(v1, 2)
0.3333333333333333

julia> quantilerank(v1, 2, method=:exc), quantilerank(v1, 2, method=:tied)
(0.36363636363636365, 0.35)

# use `skipmissing` for vectors with missing entries.
julia> quantilerank(skipmissing(v2), 4)
0.5

# use broadcasting with `Ref` to compute quantile rank for multiple values
julia> quantilerank.(Ref(v3), [4, 8])
2-element Vector{Float64}:
 0.3333333333333333
 0.8888888888888888
```
"""
function quantilerank(itr, value; method::Symbol=:inc)
    ((value isa Number && isnan(value)) || ismissing(value)) &&
        throw(ArgumentError("`value` cannot be NaN or missing"))
    any(x -> ismissing(x) || (x isa Number && isnan(x)), itr) &&
        throw(ArgumentError("`itr` cannot contain missing or NaN entries"))

    count_less = count_equal = n = 0
    greatest_smaller = smallest_greater = value
    for x in itr
        if x == value
            count_equal += 1
        elseif x < value
            count_less += 1
            if greatest_smaller == value || greatest_smaller < x
                greatest_smaller = x
            end
        else
            if smallest_greater == value || smallest_greater > x
                smallest_greater = x
            end
        end
        n += 1
    end

    n == 0 && throw(ArgumentError("`itr` is empty. Pass a collection with at least two elements"))
    n == 1 && throw(ArgumentError("`itr` has only 1 value. Pass a collection with at least two elements"))

    if method == :inc
        if greatest_smaller == value
            return 0.0
        elseif count_equal > 0
            return count_less / (n - 1)
        elseif smallest_greater == value
            return 1.0
        else
            lower = (count_less - 1) / (n - 1)
            upper = count_less / (n - 1)
            ratio = (value - greatest_smaller) / (smallest_greater - greatest_smaller)
            return lower + ratio * (upper - lower)
        end
    elseif method == :exc
        if count_less == 0 && count_equal == 0
            return 0.0
        elseif count_less == 0
            return 1.0 / (n + 1)
        elseif count_equal > 0
            return (count_less + 1) / (n + 1)
        elseif smallest_greater == value
            return 1.0
        else
            lower = count_less / (n + 1)
            upper = (count_less + 1) / (n + 1)
            ratio = (value - greatest_smaller) / (smallest_greater - greatest_smaller)
            return lower + ratio * (upper - lower)
        end
    elseif method == :compete
        if value > maximum(itr)
            return 1.0
        elseif value ≤ minimum(itr) 
            return 0.0
        else
            value ∈ itr && (count_less += 1)
            return (count_less - 1) / (n - 1)
        end 
    elseif method == :tied
        return (count_less + count_equal/2) / n
    elseif method == :strict
        return count_less / n
    elseif method == :weak
        return (count_less + count_equal) / n
    else
        throw(ArgumentError("method=:$method is not valid. Pass :inc, :exc, :compete, :tied, :strict or :weak."))
    end
end

"""
    percentilerank(itr, value; method=:inc)

Return the `q`th percentile of `value` in collection `itr`, i.e. [`quantilerank(itr, value)`](@ref) * 100.

See the [`quantilerank`](@ref) docstring for more details.
"""
percentilerank(itr, value; method::Symbol=:inc) = quantilerank(itr, value, method=method) * 100

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
    sem(x; mean=nothing)
    sem(x::AbstractArray[, weights::AbstractWeights]; mean=nothing)

Return the standard error of the mean for a collection `x`.
A pre-computed `mean` may be provided.

When not using weights, this is the (sample) standard deviation
divided by the sample size. If weights are used, the 
variance of the sample mean is calculated as follows:

* `AnalyticWeights`: Not implemented.
* `FrequencyWeights`: ``\\frac{\\sum_{i=1}^n w_i (x_i - \\bar{x_i})^2}{(\\sum w_i) (\\sum w_i - 1)}``
* `ProbabilityWeights`: ``\\frac{n}{n-1} \\frac{\\sum_{i=1}^n w_i^2 (x_i - \\bar{x_i})^2}{\\left( \\sum w_i \\right)^2}``

The standard error is then the square root of the above quantities.

# References

Carl-Erik Särndal, Bengt Swensson, Jan Wretman (1992). Model Assisted Survey Sampling.
New York: Springer. pp. 51-53.
"""
function sem(x; mean=nothing)
    if isempty(x)
        # Return the NaN of the type that we would get for a nonempty x
        T = eltype(x)
        _mean = mean === nothing ? zero(T) / 1 : mean
        z = abs2(zero(T) - _mean)
        return oftype((z + z) / 2, NaN)
    elseif mean === nothing
        n = 0
        y = iterate(x)
        value, state = y
        # Use Welford algorithm as seen in (among other places)
        # Knuth's TAOCP, Vol 2, page 232, 3rd edition.
        _mean = value / 1
        sse = real(zero(_mean))
        while y !== nothing
            value, state = y
            y = iterate(x, state)
            n += 1
            new_mean = _mean + (value - _mean) / n
            sse += realXcY(value - _mean, value - new_mean)
            _mean = new_mean
        end
    else
        n = 1
        y = iterate(x)
        value, state = y
        sse = abs2(value - mean)
        while (y = iterate(x, state)) !== nothing
            value, state = y
            n += 1
            sse += abs2(value - mean)
        end
    end
    variance = sse / (n - 1)
    return sqrt(variance / n)
end

function sem(x::AbstractArray; mean=nothing) 
    if isempty(x)
        # Return the NaN of the type that we would get for a nonempty x
        T = eltype(x)
        _mean = mean === nothing ? zero(T) / 1 : mean
        z = abs2(zero(T) - _mean)
        return oftype((z + z) / 2, NaN)
    end
    return sqrt(var(x; mean=mean, corrected=true) / length(x))
end

function sem(x::AbstractArray, weights::UnitWeights; mean=nothing)
    if length(x) ≠ length(weights)
        throw(DimensionMismatch("array and weights do not have the same length"))
    end
    return sem(x; mean=mean)
end


# Weighted methods for the above
sem(x::AbstractArray, weights::FrequencyWeights; mean=nothing) =
    sqrt(var(x, weights; mean=mean, corrected=true) / sum(weights))

function sem(x::AbstractArray, weights::ProbabilityWeights; mean=nothing)
    if isempty(x)
        # Return the NaN of the type that we would get for a nonempty x
        return var(x, weights; mean=mean, corrected=true) / 0
    else
        _mean = mean === nothing ? Statistics.mean(x, weights) : mean
        # sum of squared errors = sse
        sse = sum(Broadcast.instantiate(Broadcast.broadcasted(x, weights) do x_i, w
            return abs2(w * (x_i - _mean))
        end))
        n = count(!iszero, weights)
        return sqrt(sse * n / (n - 1)) / sum(weights)
    end
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
    isempty(x) && throw(ArgumentError("mad is not defined for empty arrays"))
    T = eltype(x)
    # Knowing the eltype allows allocating a single array able to hold both original values
    # and differences from the center, instead of two arrays
    S = isconcretetype(T) ? promote_type(T, typeof(middle(zero(T)))) : T
    x2 = x isa AbstractArray ? copyto!(similar(x, S), x) : collect(S, x)
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
function entropy(p)
    if isempty(p)
        throw(ArgumentError("empty collections are not supported since they do not " *
                            "represent proper probability distributions"))
    end
    return -sum(xlogx, p)
end

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
function crossentropy(p::AbstractArray{<:Real}, q::AbstractArray{<:Real})
    length(p) == length(q) || throw(DimensionMismatch("Inconsistent array length."))

    # handle empty collections
    if isempty(p)
        Base.depwarn(
            "support for empty collections will be removed since they do not " *
            "represent proper probability distributions",
            :crossentropy,
        )
        # return zero for empty arrays
        return xlogy(zero(eltype(p)), zero(eltype(q)))
    end

    # use pairwise summation (https://github.com/JuliaLang/julia/pull/31020)
    broadcasted = Broadcast.broadcasted(xlogy, vec(p), vec(q))
    return - sum(Broadcast.instantiate(broadcasted))
end

crossentropy(p::AbstractArray{<:Real}, q::AbstractArray{<:Real}, b::Real) =
    crossentropy(p,q) / log(b)


"""
    kldivergence(p, q, [b])

Compute the Kullback-Leibler divergence from `q` to `p`,
also called the relative entropy of `p` with respect to `q`,
that is the sum `pᵢ * log(pᵢ / qᵢ)`. Optionally a real number `b`
can be specified such that the divergence is scaled by `1/log(b)`.
"""
function kldivergence(p::AbstractArray{<:Real}, q::AbstractArray{<:Real})
    length(p) == length(q) || throw(DimensionMismatch("Inconsistent array length."))

    # handle empty collections
    if isempty(p)
        Base.depwarn(
            "support for empty collections will be removed since they do not "*
            "represent proper probability distributions",
            :kldivergence,
        )
        # return zero for empty arrays
        pzero = zero(eltype(p))
        qzero = zero(eltype(q))
        return xlogy(pzero, zero(pzero / qzero))
    end

    # use pairwise summation (https://github.com/JuliaLang/julia/pull/31020)
    broadcasted = Broadcast.broadcasted(vec(p), vec(q)) do pi, qi
        # handle pi = qi = 0, otherwise `NaN` is returned
        piqi = iszero(pi) && iszero(qi) ? zero(pi / qi) : pi / qi
        return xlogy(pi, piqi)
    end
    return sum(Broadcast.instantiate(broadcasted))
end

kldivergence(p::AbstractArray{<:Real}, q::AbstractArray{<:Real}, b::Real) =
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
    qs = if ns == 0
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
DataAPI.describe(x) = describe(stdout, x)
function DataAPI.describe(io::IO, a::AbstractArray{T}) where T<:Union{Real,Missing}
    show(io, summarystats(a))
    println(io, "Type:           $(string(eltype(a)))")
end
function DataAPI.describe(io::IO, a::AbstractArray)
    println(io, "Summary Stats:")
    println(io, "Length:         $(length(a))")
    println(io, "Type:           $(string(eltype(a)))")
    println(io, "Number Unique:  $(length(unique(a)))")
    return
end
