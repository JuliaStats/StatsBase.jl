##### Weighted var & std

## var

"""
    varm(x, wv::WeightVec, m, [dim])

Return the variance of a real-valued array `x` with a known mean `m`, optionally
over a dimension `dim`. The weighting vector `wv` specifies frequency weights
(also called case weights) for the result.

This function differs from its counterpart in Base in that Bessel's correction
is not used. That is, here the denominator for the variance is `sum(wv)`,
whereas it's `length(x)-1` in `Base.varm`. The impact is that this is not a
weighted estimate of the population variance based on the sample; it's the weighted
variance of the sample.
"""
Base.varm(v::RealArray, wv::WeightVec, m::Real) = _moment2(v, wv, m)

"""
    var(x, wv::WeightVec, [dim]; mean=nothing)

Return the variance of a real-valued array `x`, optionally over a dimension `dim`.
The weighting vector `wv` specifies frequency weights (also called case weights)
for the estimate.

This function differs from its counterpart in Base in that Bessel's correction
is not used. That is, here the denominator for the variance is `sum(wv)`,
whereas it's `length(x)-1` in `Base.var`. The impact is that this is not a
weighted estimate of the population variance based on the sample; it's the weighted
variance of the sample.
"""
function Base.var(v::RealArray, wv::WeightVec; mean=nothing)
    mean == 0 ? Base.varm(v, wv, 0) :
    mean == nothing ? varm(v, wv, Base.mean(v, wv)) :
    varm(v, wv, mean)
end

## var along dim

Base.varm!(R::AbstractArray, A::RealArray, wv::WeightVec, M::RealArray, dim::Int) =
    scale!(_wsum_centralize!(R, @functorize(abs2), A, values(wv), M, dim, true), inv(sum(wv)))

function var!(R::AbstractArray, A::RealArray, wv::WeightVec, dim::Int; mean=nothing)
    if mean == 0
        Base.varm!(R, A, wv,
            Base.reducedim_initarray(A, dim, 0, eltype(R)), dim)
    elseif mean == nothing
        Base.varm!(R, A, wv, Base.mean(A, wv, dim), dim)
    else
        # check size of mean
        for i = 1:ndims(A)
            dA = size(A,i)
            dM = size(mean,i)
            if i == dim
                dM == 1 || throw(DimensionMismatch("Incorrect size of mean."))
            else
                dM == dA || throw(DimensionMismatch("Incorrect size of mean."))
            end
        end
        Base.varm!(R, A, wv, mean, dim)
    end
end

Base.varm(A::RealArray, wv::WeightVec, M::RealArray, dim::Int) =
    @static if VERSION < v"0.6.0-dev.1121"
        Base.varm!(similar(A, Float64, Base.reduced_dims(size(A), dim)), A, wv, M, dim)
    else
        Base.varm!(similar(A, Float64, Base.reduced_indices(indices(A), dim)), A, wv, M, dim)
    end

Base.var(A::RealArray, wv::WeightVec, dim::Int; mean=nothing) =
    @static if VERSION < v"0.6.0-dev.1121"
        var!(similar(A, Float64, Base.reduced_dims(size(A), dim)), A, wv, dim; mean=mean)
    else
        var!(similar(A, Float64, Base.reduced_indices(indices(A), dim)), A, wv, dim; mean=mean)
    end

## std
"""
    stdm(v, wv::WeightVec, m, [dim])

Return the standard deviation of a real-valued array `v` with a known mean `m`,
optionally over a dimension `dim`. The weighting vector `wv` specifies frequency
weights (also called case weights) for the estimate.
"""
Base.stdm(v::RealArray, wv::WeightVec, m::Real) = sqrt(varm(v, wv, m))

"""
    std(v, wv::WeightVec, [dim]; mean=nothing)

Return the standard deviation of a real-valued array `v`, optionally over a
dimension `dim`. The weighting vector `wv` specifies frequency weights (also
called case weights) for the estimate.
"""
Base.std(v::RealArray, wv::WeightVec; mean=nothing) = sqrt.(var(v, wv; mean=mean))

Base.stdm(v::RealArray, m::RealArray, dim::Int) = Base.sqrt!(varm(v, m, dim))
Base.stdm(v::RealArray, wv::WeightVec, m::RealArray, dim::Int) = sqrt.(varm(v, wv, m, dim))
Base.std(v::RealArray, wv::WeightVec, dim::Int; mean=nothing) = sqrt.(var(v, wv, dim; mean=mean))


##### Fused statistics
"""
    mean_and_var(x, [wv::WeightVec], [dim]) -> (mean, var)

Return the mean and variance of a real-valued array `x`, optionally over a dimension
`dim`, as a tuple. A weighting vector `wv` can be specified to weight the estimates.
The weights are assumed to be frequency weights, also called case weights.
"""
mean_and_var(A::RealArray) = (m = mean(A); (m, varm(A, m)))

"""
    mean_and_std(x, [wv::WeightVec], [dim]) -> (mean, std)

Return the mean and standard deviation of a real-valued array `x`, optionally
over a dimension `dim`, as a tuple. A weighting vector `wv` can be specified
to weight the estimates. The weights are assumed to be frequency weights, also
called case weights.
"""
mean_and_std(A::RealArray) = (m = mean(A); (m, stdm(A, m)))

mean_and_var(A::RealArray, wv::WeightVec) = (m = mean(A, wv); (m, varm(A, wv, m)))
mean_and_std(A::RealArray, wv::WeightVec) = (m = mean(A, wv); (m, stdm(A, wv, m)))

mean_and_var(A::RealArray, dim::Int) = (m = mean(A, dim); (m, varm(A, m, dim)))
mean_and_std(A::RealArray, dim::Int) = (m = mean(A, dim); (m, stdm(A, m, dim)))

mean_and_var(A::RealArray, wv::WeightVec, dim::Int) = (m = mean(A, wv, dim); (m, varm(A, wv, m, dim)))
mean_and_std(A::RealArray, wv::WeightVec, dim::Int) = (m = mean(A, wv, dim); (m, stdm(A, wv, m, dim)))


##### General central moment

function _moment2(v::RealArray, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        s += z * z
    end
    s / n
end

function _moment2(v::RealArray, wv::WeightVec, m::Real)
    n = length(v)
    s = 0.0
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z * z) * w[i]
    end
    s / sum(wv)
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

function _moment3(v::RealArray, wv::WeightVec, m::Real)
    n = length(v)
    s = 0.0
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z * z * z) * w[i]
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

function _moment4(v::RealArray, wv::WeightVec, m::Real)
    n = length(v)
    s = 0.0
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += abs2(z * z) * w[i]
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

function _momentk(v::RealArray, k::Int, wv::WeightVec, m::Real)
    n = length(v)
    s = 0.0
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z ^ k) * w[i]
    end
    s / sum(wv)
end


"""
    moment(v, k, [wv::WeightVec], m=mean(v))

Return the `k`th order central moment of a real-valued array `v`, optionally
specifying a weighting vector `wv` and a center `m`.
"""
function moment(v::RealArray, k::Int, m::Real)
    k == 2 ? _moment2(v, m) :
    k == 3 ? _moment3(v, m) :
    k == 4 ? _moment4(v, m) :
    _momentk(v, k, m)
end

function moment(v::RealArray, k::Int, wv::WeightVec, m::Real)
    k == 2 ? _moment2(v, wv, m) :
    k == 3 ? _moment3(v, wv, m) :
    k == 4 ? _moment4(v, wv, m) :
    _momentk(v, k, wv, m)
end

moment(v::RealArray, k::Int) = moment(v, k, mean(v))
moment(v::RealArray, k::Int, wv::WeightVec) = moment(v, k, wv, mean(v, wv))


##### Skewness and Kurtosis

# Skewness
# This is Type 1 definition according to Joanes and Gill (1998)
"""
    skewness(v, [wv::WeightVec], m=mean(v))

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

function skewness(v::RealArray, wv::WeightVec, m::Real)
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    cm2 = 0.0   # empirical 2nd centered moment (variance)
    cm3 = 0.0   # empirical 3rd centered moment
    w = values(wv)

    @inbounds for i = 1:n
        x_i = v[i]
        w_i = w[i]
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
skewness(v::RealArray, wv::WeightVec) = skewness(v, wv, mean(v, wv))

# (excessive) Kurtosis
# This is Type 1 definition according to Joanes and Gill (1998)
"""
    kurtosis(v, [wv::WeightVec], m=mean(v))

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

function kurtosis(v::RealArray, wv::WeightVec, m::Real)
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    cm2 = 0.0  # empirical 2nd centered moment (variance)
    cm4 = 0.0  # empirical 4th centered moment
    w = values(wv)

    @inbounds for i = 1 : n
        x_i = v[i]
        w_i = w[i]
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
kurtosis(v::RealArray, wv::WeightVec) = kurtosis(v, wv, mean(v, wv))

