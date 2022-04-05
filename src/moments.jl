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
function var(v::RealArray, w::AbstractWeights; mean=nothing,
                  corrected::DepBool=nothing)
    corrected = depcheck(:var, :corrected, corrected)

    if mean == nothing
        _moment2(v, w, Statistics.mean(v, w); corrected=corrected)
    else
        _moment2(v, w, mean; corrected=corrected)
    end
end

## var along dim

function var!(R::AbstractArray, A::RealArray, w::AbstractWeights, dims::Int;
              mean=nothing, corrected::DepBool=nothing)
    corrected = depcheck(:var!, :corrected, corrected)

    if mean == 0
        mean = Base.reducedim_initarray(A, dims, 0, eltype(R))
    elseif mean === nothing
        mean = Statistics.mean(A, w, dims=dims)
    else
        # check size of mean
        for i = 1:ndims(A)
            dA = size(A,i)
            dM = size(mean,i)
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

function var(A::RealArray, w::AbstractWeights, dim::Int; mean=nothing,
                  corrected::DepBool=nothing)
    corrected = depcheck(:var, :corrected, corrected)
    var!(similar(A, Float64, Base.reduced_indices(axes(A), dim)), A, w, dim;
         mean=mean, corrected=corrected)
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
std(v::RealArray, w::AbstractWeights; mean=nothing, corrected::DepBool=nothing) =
    sqrt.(var(v, w; mean=mean, corrected=depcheck(:std, :corrected, corrected)))

std(v::RealArray, w::AbstractWeights, dim::Int;
    mean=nothing, corrected::DepBool=nothing) =
    sqrt.(var(v, w, dim; mean=mean, corrected=depcheck(:std, :corrected, corrected)))

##### Fused statistics
"""
    mean_and_var(x, [w::AbstractWeights], [dim]; corrected=false) -> (mean, var)

Return the mean and variance of collection `x`. If `x` is an `AbstractArray`,
`dim` can be specified as a tuple to compute statistics over these dimensions.
A weighting vector `w` can be specified to weight the estimates.
Finally, bias correction is be applied to the variance calculation if `corrected=true`.
See [`var`](@ref) documentation for more details.
"""
function mean_and_var(x; corrected::Bool=true)
    m = mean(x)
    v = var(x, mean=m, corrected=corrected)
    m, v
end

"""
    mean_and_std(x, [w::AbstractWeights], [dim]; corrected=false) -> (mean, std)

Return the mean and standard deviation of collection `x`. If `x` is an `AbstractArray`,
`dim` can be specified as a tuple to compute statistics over these dimensions.
A weighting vector `w` can be specified to weight the estimates.
Finally, bias correction is applied to the
standard deviation calculation if `corrected=true`.
See [`std`](@ref) documentation for more details.
"""
function mean_and_std(x; corrected::Bool=true)
    m = mean(x)
    s = std(x, mean=m, corrected=corrected)
    m, s
end

function mean_and_var(x::RealArray, w::AbstractWeights; corrected::DepBool=nothing)
    m = mean(x, w)
    v = var(x, w, mean=m, corrected=depcheck(:mean_and_var, :corrected, corrected))
    m, v
end
function mean_and_std(x::RealArray, w::AbstractWeights; corrected::DepBool=nothing)
    m = mean(x, w)
    s = std(x, w, mean=m, corrected=depcheck(:mean_and_std, :corrected, corrected))
    m, s
end


function mean_and_var(x::RealArray, dim::Int; corrected::Bool=true)
    m = mean(x, dims=dim)
    v = var(x, dims=dim, mean=m, corrected=corrected)
    m, v
end
function mean_and_std(x::RealArray, dim::Int; corrected::Bool=true)
    m = mean(x, dims=dim)
    s = std(x, dims=dim, mean=m, corrected=corrected)
    m, s
end


function mean_and_var(x::RealArray, w::AbstractWeights, dims::Int;
                      corrected::DepBool=nothing)
    m = mean(x, w, dims=dims)
    v = var(x, w, dims, mean=m, corrected=depcheck(:mean_and_var, :corrected, corrected))
    m, v
end
function mean_and_std(x::RealArray, w::AbstractWeights, dims::Int;
                      corrected::DepBool=nothing)
    m = mean(x, w, dims=dims)
    s = std(x, w, dims, mean=m, corrected=depcheck(:mean_and_std, :corrected, corrected))
    m, s
end



##### General central moment
function _moment2(v::RealArray, m::Real; corrected=false)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        s += z * z
    end
    varcorrection(n, corrected) * s
end

function _moment2(v::RealArray, wv::AbstractWeights, m::Real; corrected=false)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z * z) * wv[i]
    end

    varcorrection(wv, corrected) * s
end

function _moment3(v::RealArray, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        s += z * z * z
    end
    s / n
end

function _moment3(v::RealArray, wv::AbstractWeights, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z * z * z) * wv[i]
    end
    s / sum(wv)
end

function _moment4(v::RealArray, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        s += abs2(z * z)
    end
    s / n
end

function _moment4(v::RealArray, wv::AbstractWeights, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += abs2(z * z) * wv[i]
    end
    s / sum(wv)
end

function _momentk(v::RealArray, k::Int, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        s += (z ^ k)
    end
    s / n
end

function _momentk(v::RealArray, k::Int, wv::AbstractWeights, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z ^ k) * wv[i]
    end
    s / sum(wv)
end


"""
    moment(v, k, [wv::AbstractWeights], m=mean(v))

Return the `k`th order central moment of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a center `m`.
"""
function moment(v::RealArray, k::Int, m::Real)
    k == 2 ? _moment2(v, m) :
    k == 3 ? _moment3(v, m) :
    k == 4 ? _moment4(v, m) :
    _momentk(v, k, m)
end

function moment(v::RealArray, k::Int, wv::AbstractWeights, m::Real)
    k == 2 ? _moment2(v, wv, m) :
    k == 3 ? _moment3(v, wv, m) :
    k == 4 ? _moment4(v, wv, m) :
    _momentk(v, k, wv, m)
end

moment(v::RealArray, k::Int) = moment(v, k, mean(v))
function moment(v::RealArray, k::Int, wv::AbstractWeights)
    moment(v, k, wv, mean(v, wv))
end


##### Skewness and Kurtosis

# Skewness
# This is Type 1 definition according to Joanes and Gill (1998)
"""
    skewness(v, [wv::AbstractWeights], m=mean(v))

Compute the standardized skewness of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a center `m`.
"""
function skewness(v::RealArray, m::Real)
    n = length(v)
    cm2 = 0.0   # empirical 2nd centered moment (variance)
    cm3 = 0.0   # empirical 3rd centered moment
    for i = 1:n
        @inbounds z = v[i] - m
        z2 = z * z

        cm2 += z2
        cm3 += z2 * z
    end
    cm3 /= n
    cm2 /= n
    return cm3 / sqrt(cm2 * cm2 * cm2)  # this is much faster than cm2^1.5
end

function skewness(v::RealArray, wv::AbstractWeights, m::Real)
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    cm2 = 0.0   # empirical 2nd centered moment (variance)
    cm3 = 0.0   # empirical 3rd centered moment

    @inbounds for i = 1:n
        x_i = v[i]
        w_i = wv[i]
        z = x_i - m
        z2w = z * z * w_i
        cm2 += z2w
        cm3 += z2w * z
    end
    sw = sum(wv)
    cm3 /= sw
    cm2 /= sw
    return cm3 / sqrt(cm2 * cm2 * cm2)  # this is much faster than cm2^1.5
end

skewness(v::RealArray) = skewness(v, mean(v))
skewness(v::RealArray, wv::AbstractWeights) = skewness(v, wv, mean(v, wv))

# (excessive) Kurtosis
# This is Type 1 definition according to Joanes and Gill (1998)
"""
    kurtosis(v, [wv::AbstractWeights], m=mean(v))

Compute the excess kurtosis of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a center `m`.
"""
function kurtosis(v::RealArray, m::Real)
    n = length(v)
    cm2 = 0.0  # empirical 2nd centered moment (variance)
    cm4 = 0.0  # empirical 4th centered moment
    for i = 1:n
        @inbounds z = v[i] - m
        z2 = z * z
        cm2 += z2
        cm4 += z2 * z2
    end
    cm4 /= n
    cm2 /= n
    return (cm4 / (cm2 * cm2)) - 3.0
end

function kurtosis(v::RealArray, wv::AbstractWeights, m::Real)
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    cm2 = 0.0  # empirical 2nd centered moment (variance)
    cm4 = 0.0  # empirical 4th centered moment

    @inbounds for i = 1 : n
        x_i = v[i]
        w_i = wv[i]
        z = x_i - m
        z2 = z * z
        z2w = z2 * w_i
        cm2 += z2w
        cm4 += z2w * z2
    end
    sw = sum(wv)
    cm4 /= sw
    cm2 /= sw
    return (cm4 / (cm2 * cm2)) - 3.0
end

kurtosis(v::RealArray) = kurtosis(v, mean(v))
kurtosis(v::RealArray, wv::AbstractWeights) = kurtosis(v, wv, mean(v, wv))
