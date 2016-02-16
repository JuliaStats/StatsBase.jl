# Descriptive Statistics


#############################
#
#   Location
#
#############################

# Geometric mean
"""
    geomean(a::RealArray)

Compute the geometric mean of `a`.
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
    harmmean(a::RealArray)

Compute the harmonic mean of `a`.
"""
function harmmean(a::RealArray)
    s = 0.0
    n = length(a)
    for i in 1 : n
        @inbounds s += inv(a[i])
    end
    return n / s
end


# Trimmed mean
"""
    trimmean(a::RealArray)

Compute the trimmed mean of `a`.
"""
function trimmean(x::RealArray, p::Real)
    n = length(x)
    n > 0 || error("x can not be empty.")
    0 <= p < 1 || error("p must be non-negative and less than 1.")
    rn = min(round(Int, n * p), n-1)

    sx = sort(x)
    nl = rn >> 1
    nh = (rn - nl)
    s = 0.0
    for i = (1+nl) : (n-nh)
        @inbounds s += sx[i]
    end
    return s / (n - rn)
end

# compute mode, given the range of integer values
function mode{T<:Integer}(a::AbstractArray{T}, rgn::UnitRange{T})
    isempty(a) && error("mode: input array cannot be empty.")
    len = length(a)
    r0 = rgn[1]
    r1 = rgn[end]
    cnts = zeros(Int, length(rgn))
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

function modes{T<:Integer}(a::AbstractArray{T}, rgn::UnitRange{T})
    r0 = rgn[1]
    r1 = rgn[end]
    n = length(rgn)
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
            push!(ms, rgn[i])
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

"""
    mode{T}(a::AbstractArray{T}, rgn::UnitRange{T})

### Args:
* `a`: An arbitrary array of type `AbstractArray`.
* `rgn`: A range of integer values, which is optional.

Compute the first mode of `a`, over the range `rng` if specified.
"""
function mode end

"""
    modes{T}(a::AbstractArray{T}, rgn::UnitRange{T})

### Args:
* `a`: An arbitrary array of type `AbstractArray`.
* `rgn`: A range of integer values, which is optional.

Compute all the modes of `a`, over the range `rng` if specified.
"""
function mode end


#############################
#
#   quantile and friends
#
#############################

"""
    percentile(v::AbstractArray{T}, p)

### Args:
* `v`: An `AbstracArray`.
* `p`: Percentage value.

Computes the quantile using percentage (instead of fraction) as argument.
"""
percentile{T<:Real}(v::AbstractArray{T}, p) = quantile(v, p * 0.01)

quantile{T<:Real}(v::AbstractArray{T}) = quantile(v, [.0, .25, .5, .75, 1.0])
"""
    nquantile(v::AbstractArray{T}, n::Integer)

### Args:
* `v`: An array of type `AbstractArray`.
* `n`: Integer value to compute the quantile at.

Compute quantiles at [0:n]/n. For example, nquantiles(x, 5) returns a vector of quantiles, respectively at 0.0, 0.2, 0.4, 0.6, 0.8, 1.0.
"""
nquantile{T<:Real}(v::AbstractArray{T}, n::Integer) = quantile(v, (0:n)/n)


#############################
#
#   Dispersion
#
#############################

# span, i.e. the range minimum(x):maximum(x)
"""
    span(x::AbstractArray{T})

Returns the range of `x`, i.e., minimum(x):maximum(x).
"""
span{T<:Integer}(x::AbstractArray{T}) = ((a, b) = extrema(x); a:b)

# Variation coefficient: std / mean
variation{T<:Real}(x::AbstractArray{T}, m::Real) = stdm(x, m) / m
variation{T<:Real}(x::AbstractArray{T}) = variation(x, mean(x))
"""
    variation(x::AbstractArray{T}, m::Real)

### Args:
* `x`: An array of type `AbstractArray`.
* `m`: An optional paramater, which specifies the mean.

Computes the ratio of standard deviation to mean.
"""
function variation end

# Standard error of the mean: std(a) / sqrt(len)
"""
    sem(a::AbstractArray{T})

### Args:
* `a`: An array of type `AbstractArray`.

Computes the standard error of the mean, i.e. sqrt(var / n).
"""
sem{T<:Real}(a::AbstractArray{T}) = sqrt(var(a) / length(a))

# Median absolute deviation
mad{T<:Real}(v::AbstractArray{T}, args...;arg...) = mad!(copy(v), args...;arg...)
mad{T<:Real}(v::Range{T}, args...;arg...) = mad!([v;], args...;arg...)
"""
    mad(x, center;constant=1.4826)
### Args:
* `x`: An array of type `AbstractArray`.
* `center`: Optional parameter which specifies the center of `x`.
* `constant`: Optional paramter, Constant scale factor, default is 1.4826.

Compute the median absolute deviation of x. One can optionally supply the center. By default, constant=1.4826 for consistent estimation of the standard deviation of a normal distribution.
"""
function mad end    
function mad!{T<:Real}(v::AbstractArray{T}, center::Real=median!(v); constant::Real=1.4826)
    for i in 1:length(v)
        @inbounds v[i] = abs(v[i]-center)
    end
    constant * median!(v)
end

# Interquartile range
"""
    iqr(v::AbstractArray{T})

### Args:
* `v`: An array of type `AbstractArray`.

Compute the interquartile range of x, i.e. quantile(x, 0.75) - quantile(x, 0.25).
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

function zscore{T<:Real}(X::AbstractArray{T}, μ::Real, σ::Real)
    ZT = typeof((zero(T) - zero(μ)) / one(σ))
    _zscore!(Array(ZT, size(X)), X, μ, σ)
end

function zscore{T<:Real,U<:Real,S<:Real}(X::AbstractArray{T}, μ::AbstractArray{U}, σ::AbstractArray{S})
    _zscore_chksize(X, μ, σ)
    ZT = typeof((zero(T) - zero(U)) / one(S))
    _zscore!(Array(ZT, size(X)), X, μ, σ)
end

zscore{T<:Real}(X::AbstractArray{T}) = ((μ, σ) = mean_and_std(X); zscore(X, μ, σ))
zscore{T<:Real}(X::AbstractArray{T}, dim::Int) = ((μ, σ) = mean_and_std(X, dim); zscore(X, μ, σ))

"""
    zscore(X::AbstractArray, μ::AbstractArray, σ::AbstractArray)

### Args:
* `X`: The input matrix to compute the Z-score.
* `μ`: The mean, this could be a scalar or a an array.
* `σ`: The standard deviation, this could be a scalar or a an array.

Compute the Z-scores, given the mean μ and standard deviation σ, which is defined as (x - μ) / σ. This function returns an array Z of the same size as X. Here, μ and σ should be both scalars or both arrays.
"""
function zscore end

"""
    zscore(X::AbstractArray, μ::AbstractArray, σ::AbstractArray)

### Args:
* `X`: The input matrix to compute the Z-score.
* `μ`: The mean, this could be a scalar or a an array.
* `σ`: The standard deviation, this could be a scalar or a an array.

Compute the Z-scores inline, given the mean μ and standard deviation σ, which is defined as (x - μ) / σ. Here, μ and σ should be both scalars or both arrays.
"""
function zscore! end


#############################
#
#   entropy and friends
#
#############################

function entropy{T<:Real}(p::AbstractArray{T})
    s = 0.
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
    entropy(p::AbstractArray{T}, b::Real)

### Args:
* `p`: Probability vector of type `AbstractArray{T}`.
* `b`: Optional paramater, which specifies the base of logarithm, default is natural logarithm.

Compute the entropy of the probability vector p using logarithms of base b (e.g. entropy(p,2) returns the entropy in bits).
"""
function entropy end   

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
    crossentropy(p::AbstractArray{T}, q::AbstractArray{T}, b::Real)

### Args:
* `p`: First probability vector.
* `q`: Second probability vector.
* `b`: Base of logarithm, an optional parameter.

Compute the cross entropy between p and q using logarithms of base b.
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
"""
    kldivergence(p::AbstractArray{T}, q::AbstractArray{T}, b::Real)

### Args:
* `p`: First probability vector.
* `q`: Second probability vector.
* `b`: Base of logarithm, an optional parameter.

Computes the Kullback-Leiver divergence between the two probability vectors. 
"""

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
    summarystats(a::AbstractArray)

Compute a set of statistics over x and return a struct of type SummaryStats defined as below:
```
immutable SummaryStats{T<:AbstractFloat}
    mean::T
    min::T
    q25::T
    median::T
    q75::T
    max::T
end
```
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
    @printf(io, "Mean:         %.6f\n", ss.mean)
    @printf(io, "Minimum:      %.6f\n", ss.min)
    @printf(io, "1st Quartile: %.6f\n", ss.q25)
    @printf(io, "Median:       %.6f\n", ss.median)
    @printf(io, "3rd Quartile: %.6f\n", ss.q75)
    @printf(io, "Maximum:      %.6f\n", ss.max)
end

"""
    describe(a::AbstractArray)

Print the summary statistics of `a`.
"""
describe{T<:Real}(a::AbstractArray{T}) = show(summarystats(a))
