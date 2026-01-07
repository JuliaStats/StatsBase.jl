##### Weighted var & std

## var
"""
    var(x::AbstractArray, w::AbstractWeights, [dim]; mean=nothing, corrected=false)

Compute the variance of a real-valued array `x`, optionally over a dimension `dim`.
Observations in `x` are weighted using weight vector `w`.
The uncorrected (when `corrected=false`) sample variance is defined as:
```math
\\frac{1}{\\sum{w}} \\sum_{i=1}^n {w_i\\left({x_i - μ}\\right)^2 }
```
where ``n`` is the length of the input and ``μ`` is the mean.
The unbiased estimate (when `corrected=true`) of the population variance is computed by
replacing ``\\frac{1}{\\sum{w}}`` with a factor dependent on the type of weights used:
* `AnalyticWeights`: ``\\frac{1}{\\sum w - \\sum {w^2} / \\sum w}``
* `FrequencyWeights`: ``\\frac{1}{\\sum{w} - 1}``
* `ProbabilityWeights`: ``\\frac{n}{(n - 1) \\sum w}`` where ``n`` equals `count(!iszero, w)`
* `Weights`: `ArgumentError` (bias correction not supported)
"""
function var(v::AbstractArray{<:Real}, w::AbstractWeights; mean=nothing,
             corrected::Union{Bool,Nothing}=nothing)
    length(w) == length(v) || throw(DimensionMismatch("Inconsistent array lengths."))
    corrected = depcheck(:var, :corrected, corrected)
    if mean === nothing
        mean = Statistics.mean(v, w)
    end
    return _moment2(v, w, mean; corrected)
end
function var(v::AbstractArray{<:Real}, w::UnitWeights; mean=nothing,
             corrected::Union{Bool,Nothing}=nothing)
    length(w) == length(v) || throw(DimensionMismatch("Inconsistent array lengths."))
    corrected = depcheck(:var, :corrected, corrected)
    return var(v; mean, corrected)
end

## var along dim

function var!(R::AbstractArray, A::AbstractArray{<:Real}, w::AbstractWeights, dims::Int;
              mean=nothing, corrected::Union{Bool,Nothing}=nothing)
    corrected = depcheck(:var!, :corrected, corrected)

    if mean == 0
        mean = Base.reducedim_initarray(A, dims, 0, eltype(R))
    elseif mean === nothing
        mean = Statistics.mean(A, w; dims=dims)
    else
        # check size of mean
        for i in 1:ndims(A)
            dA = size(A, i)
            dM = size(mean, i)
            if i == dims
                dM == 1 || throw(DimensionMismatch("Incorrect size of mean."))
            else
                dM == dA || throw(DimensionMismatch("Incorrect size of mean."))
            end
        end
    end
    return rmul!(_wsum_centralize!(R, abs2, A, convert(Vector, w), mean, dims, true),
                 varcorrection(w, corrected))
end

function var(A::AbstractArray{<:Real}, w::AbstractWeights, dim::Int; mean=nothing,
             corrected::Union{Bool,Nothing}=nothing)
    corrected = depcheck(:var, :corrected, corrected)
    if mean === nothing
        z = (zero(eltype(w)) * zero(eltype(A))^2) / zero(eltype(w))
    else
        z = (zero(eltype(w)) * zero(zero(eltype(A)) - zero(eltype(mean)))^2) /
            zero(eltype(w))
    end
    return var!(similar(A, typeof(z), Base.reduced_indices(axes(A), dim)), A, w, dim;
                mean=mean, corrected=corrected)
end
function var(v::AbstractArray{<:Real}, w::UnitWeights, dim::Int; mean=nothing,
             corrected::Union{Bool,Nothing}=nothing)
    length(w) == length(v) || throw(DimensionMismatch("Inconsistent array lengths."))
    corrected = depcheck(:var, :corrected, corrected)
    return var(v; mean, corrected, dims=dim)
end

## std
"""
    std(x::AbstractArray, w::AbstractWeights, [dim]; mean=nothing, corrected=false)

Compute the standard deviation of a real-valued array `x`,
optionally over a dimension `dim`. Observations in `x` are weighted using weight vector `w`.
The uncorrected (when `corrected=false`) sample standard deviation is defined as:
```math
\\sqrt{\\frac{1}{\\sum{w}} \\sum_{i=1}^n {w_i\\left({x_i - μ}\\right)^2 }}
```
where ``n`` is the length of the input and ``μ`` is the mean.
The unbiased estimate (when `corrected=true`) of the population standard deviation is
computed by replacing ``\\frac{1}{\\sum{w}}`` with a factor dependent on the type of
weights used:
* `AnalyticWeights`: ``\\frac{1}{\\sum w - \\sum {w^2} / \\sum w}``
* `FrequencyWeights`: ``\\frac{1}{\\sum{w} - 1}``
* `ProbabilityWeights`: ``\\frac{n}{(n - 1) \\sum w}`` where ``n`` equals `count(!iszero, w)`
* `Weights`: `ArgumentError` (bias correction not supported)
"""
function std(v::AbstractArray{<:Real}, w::AbstractWeights; mean=nothing,
             corrected::Union{Bool,Nothing}=nothing)
    return sqrt.(var(v, w; mean=mean, corrected=depcheck(:std, :corrected, corrected)))
end

function std(v::AbstractArray{<:Real}, w::AbstractWeights, dim::Int;
             mean=nothing, corrected::Union{Bool,Nothing}=nothing)
    return sqrt.(var(v, w, dim; mean=mean, corrected=depcheck(:std, :corrected, corrected)))
end

##### Fused statistics
"""
    mean_and_var(x, [w::AbstractWeights], [dim]; corrected=true) -> (mean, var)

Return the mean and variance of collection `x`. If `x` is an `AbstractArray`,
`dim` can be specified as a tuple to compute statistics over these dimensions.
A weighting vector `w` can be specified to weight the estimates.
Finally, bias correction is be applied to the variance calculation if `corrected=true`.
See [`var`](@ref) documentation for more details.
"""
function mean_and_var(x; corrected::Bool=true)
    m = mean(x)
    v = var(x; mean=m, corrected=corrected)
    return m, v
end

"""
    mean_and_std(x, [w::AbstractWeights], [dim]; corrected=true) -> (mean, std)

Return the mean and standard deviation of collection `x`. If `x` is an `AbstractArray`,
`dim` can be specified as a tuple to compute statistics over these dimensions.
A weighting vector `w` can be specified to weight the estimates.
Finally, bias correction is applied to the
standard deviation calculation if `corrected=true`.
See [`std`](@ref) documentation for more details.
"""
function mean_and_std(x; corrected::Bool=true)
    m = mean(x)
    s = std(x; mean=m, corrected=corrected)
    return m, s
end

function mean_and_var(x::AbstractArray{<:Real}, w::AbstractWeights;
                      corrected::Union{Bool,Nothing}=nothing)
    m = mean(x, w)
    v = var(x, w; mean=m, corrected=depcheck(:mean_and_var, :corrected, corrected))
    return m, v
end
function mean_and_std(x::AbstractArray{<:Real}, w::AbstractWeights;
                      corrected::Union{Bool,Nothing}=nothing)
    m = mean(x, w)
    s = std(x, w; mean=m, corrected=depcheck(:mean_and_std, :corrected, corrected))
    return m, s
end

function mean_and_var(x::AbstractArray{<:Real}, dim::Int; corrected::Bool=true)
    m = mean(x; dims=dim)
    v = var(x; dims=dim, mean=m, corrected=corrected)
    return m, v
end
function mean_and_std(x::AbstractArray{<:Real}, dim::Int; corrected::Bool=true)
    m = mean(x; dims=dim)
    s = std(x; dims=dim, mean=m, corrected=corrected)
    return m, s
end

function mean_and_var(x::AbstractArray{<:Real}, w::AbstractWeights, dims::Int;
                      corrected::Union{Bool,Nothing}=nothing)
    m = mean(x, w; dims=dims)
    v = var(x, w, dims; mean=m, corrected=depcheck(:mean_and_var, :corrected, corrected))
    return m, v
end
function mean_and_std(x::AbstractArray{<:Real}, w::AbstractWeights, dims::Int;
                      corrected::Union{Bool,Nothing}=nothing)
    m = mean(x, w; dims=dims)
    s = std(x, w, dims; mean=m, corrected=depcheck(:mean_and_std, :corrected, corrected))
    return m, s
end

##### General central moment
function _moment2(v::AbstractArray{<:Real}, m::Real; corrected::Bool)
    n = length(v)
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        s = z^2
    else
        s = let m = m
            sum(v) do vi
                zi = vi - m
                return zi^2
            end
        end
    end
    return oftype(float(s), varcorrection(n, corrected)) * s
end

function _moment2(v::AbstractArray{<:Real}, wv::AbstractWeights, m::Real; corrected::Bool)
    if isempty(v)
        z = zero(zero(eltype(v)) - m)
        s = z^2 * zero(eltype(wv))
    else
        broadcasted = let m = m
            Broadcast.broadcasted(vec(v), wv) do vi, wvi
                zi = vi - m
                return zi^2 * wvi
            end
        end
        s = sum(Broadcast.instantiate(broadcasted))
    end
    return oftype(float(s), varcorrection(wv, corrected)) * s
end

function _moment3(v::AbstractArray{<:Real}, m::Real)
    n = length(v)
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        s = z^3
    else
        s = let m = m
            sum(v) do vi
                zi = vi - m
                return zi^3
            end
        end
    end
    return s / n
end

function _moment3(v::AbstractArray{<:Real}, wv::AbstractWeights, m::Real)
    if isempty(v)
        z = zero(zero(eltype(v)) - m)
        s = zero(z^3 * zero(eltype(wv)))
    else
        broadcasted = let m = m
            Broadcast.broadcasted(vec(v), wv) do vi, wvi
                zi = vi - m
                return zi^3 * wvi
            end
        end
        s = sum(Broadcast.instantiate(broadcasted))
    end
    return s / sum(wv)
end

function _moment4(v::AbstractArray{<:Real}, m::Real)
    n = length(v)
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        s = zero(z^4)
    else
        s = let m = m
            sum(v) do vi
                zi = vi - m
                return zi^4
            end
        end
    end
    return s / n
end

function _moment4(v::AbstractArray{<:Real}, wv::AbstractWeights, m::Real)
    if isempty(v)
        z = zero(zero(eltype(v)) - m)
        s = zero(z^4 * zero(eltype(wv)))
    else
        broadcasted = let m = m
            Broadcast.broadcasted(vec(v), wv) do vi, wvi
                zi = vi - m
                return zi^4 * wvi
            end
        end
        s = sum(Broadcast.instantiate(broadcasted))
    end
    return s / sum(wv)
end

function _momentk(v::AbstractArray{<:Real}, k::Int, m::Real)
    n = length(v)
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        s = zero(z^k)
    else
        s = let m = m, k = k
            sum(v) do vi
                zi = vi - m
                return zi^k
            end
        end
    end
    return s / n
end

function _momentk(v::AbstractArray{<:Real}, k::Int, wv::AbstractWeights, m::Real)
    if isempty(v)
        z = zero(zero(eltype(v)) - m)
        s = zero(z^k * zero(eltype(wv)))
    else
        broadcasted = let m = m, k = k
            Broadcast.broadcasted(vec(v), wv) do vi, wvi
                zi = vi - m
                return zi^k * wvi
            end
        end
        s = sum(Broadcast.instantiate(broadcasted))
    end
    return s / sum(wv)
end

"""
    moment(v, k, [wv::AbstractWeights], m=mean(v))

Return the `k`th order central moment of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a center `m`.
"""
function moment(v::AbstractArray{<:Real}, k::Int, m::Real=mean(v))
    return k == 2 ? _moment2(v, m; corrected=false) :
           k == 3 ? _moment3(v, m) :
           k == 4 ? _moment4(v, m) :
           _momentk(v, k, m)
end

function moment(v::AbstractArray{<:Real}, k::Int, wv::AbstractWeights, m::Real=mean(v, wv))
    length(wv) == length(v) || throw(DimensionMismatch("Inconsistent array lengths."))
    return k == 2 ? _moment2(v, wv, m; corrected=false) :
           k == 3 ? _moment3(v, wv, m) :
           k == 4 ? _moment4(v, wv, m) :
           _momentk(v, k, wv, m)
end
function moment(v::AbstractArray{<:Real}, k::Int, wv::UnitWeights, m::Real)
    length(wv) == length(v) || throw(DimensionMismatch("Inconsistent array lengths."))
    return moment(v, k, m)
end

##### Skewness and Kurtosis

# Skewness
# This is Type 1 definition according to Joanes and Gill (1998)
"""
    skewness(v, [wv::AbstractWeights], m=mean(v))

Compute the standardized skewness of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a center `m`.
"""
function skewness(v::AbstractArray{<:Real}, m::Real=mean(v))
    n = length(v)
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        cm2 = z^2 # empirical 2nd centered moment (variance)
        cm3 = cm2 * z # empirical 3rd centered moment
    else
        cm2, cm3 = let m = m
            mapreduce(_add, v) do vi
                zi = vi - m
                z2i = zi^2
                z3i = z2i * zi
                return z2i, z3i
            end
        end
    end
    return (cm3/n) / sqrt((cm2/n)^3)  # this is much faster than cm2^1.5
end

function skewness(v::AbstractArray{<:Real}, wv::AbstractWeights, m::Real=mean(v, wv))
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        cm2 = z^2 * zero(eltype(wv)) # empirical 2nd centered moment (variance)
        cm3 = cm2 * z # empirical 3rd centered moment
    else
        broadcasted = let m = m
            Broadcast.broadcasted(vec(v), wv) do vi, wvi
                zi = vi - m
                z2wi = zi^2 * wvi
                z3wi = z2wi * zi
                return z2wi, z3wi
            end
        end
        cm2, cm3 = reduce(_add, Broadcast.instantiate(broadcasted))
    end
    sw = sum(wv)
    return (cm3/sw) / sqrt((cm2/sw)^3)  # this is much faster than cm2^1.5
end
function skewness(v::AbstractArray{<:Real}, wv::UnitWeights, m::Real)
    length(wv) == length(v) || throw(DimensionMismatch("Inconsistent array lengths."))
    return skewness(v, m)
end

# (excessive) Kurtosis
# This is Type 1 definition according to Joanes and Gill (1998)
"""
    kurtosis(v, [wv::AbstractWeights], m=mean(v))

Compute the excess kurtosis of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a center `m`.
"""
function kurtosis(v::AbstractArray{<:Real}, m::Real=mean(v))
    n = length(v)
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        cm2 = z^2 # empirical 2nd centered moment (variance)
        cm4 = cm2^2 # empirical 4th centered moment
    else
        cm2, cm4 = let m = m
            mapreduce(_add, v) do vi
                zi = vi - m
                z2i = zi^2
                z4i = z2i^2
                return z2i, z4i
            end
        end
    end
    return (cm4/n) / (cm2/n)^2 - 3
end

function kurtosis(v::AbstractArray{<:Real}, wv::AbstractWeights, m::Real=mean(v, wv))
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    if iszero(n)
        z = zero(zero(eltype(v)) - m)
        z2 = z^2
        cm2 = z2 * zero(eltype(wv)) # empirical 2nd centered moment (variance)
        cm4 = cm2 * z2 # empirical 4th centered moment
    else
        broadcasted = let m = m
            Broadcast.broadcasted(vec(v), wv) do vi, wvi
                zi = vi - m
                z2i = zi^2
                z2wi = z2i * wvi
                z4wi = z2wi * z2i
                return z2wi, z4wi
            end
        end
        cm2, cm4 = reduce(_add, Broadcast.instantiate(broadcasted))
    end
    sw = sum(wv)
    return (cm4/sw) / (cm2/sw)^2 - 3
end
function kurtosis(v::AbstractArray{<:Real}, wv::UnitWeights, m::Real)
    length(wv) == length(v) || throw(DimensionMismatch("Inconsistent array lengths."))
    return kurtosis(v, m)
end

"""
    cumulant(v, k, [wv::AbstractWeights], m=mean(v))

Return the `k`th order cumulant of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a pre-computed mean `m`.

If `k` is a range of `Integer`s, then return all the cumulants of orders in this range as a vector.

This quantity is calculated using a recursive definition on lower-order cumulants and central moments.

Reference: Smith, P. J. 1995. A Recursive Formulation of the Old Problem of Obtaining
Moments from Cumulants and Vice Versa. The American Statistician, 49(2), 217–218.
https://doi.org/10.2307/2684642
"""
function cumulant(v::AbstractArray{<:Real}, krange::Union{Integer,AbstractRange{<:Integer}},
                  wv::AbstractWeights,
                  m::Real=mean(v, wv))
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    kmin, kmax = extrema(krange)
    if kmin <= 0
        throw(ArgumentError("Cumulant orders must be strictly positive."))
    end
    cmoms = [moment(v, i, wv, m) for i in 2:kmax]
    cumls = Vector{eltype(cmoms)}(undef, kmax)
    cumls[1] = m
    for i in 2:kmax
        kn = cmoms[i-1]
        for j in 2:(i - 2)
            kn -= binomial(i-1, j)*cmoms[j-1]*cumls[i-j]
        end
        cumls[i] = kn
    end
    return cumls[krange]
end

function cumulant(v::AbstractArray{<:Real}, krange::Union{Integer,AbstractRange{<:Integer}},
                  m::Real=mean(v))
    return cumulant(v, krange, uweights(length(v)), m)
end
