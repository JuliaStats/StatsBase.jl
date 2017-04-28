##### Weighted var & std

## var

"""
    varm(x, wv::AbstractWeights, m, [dim])

Return the variance of a real-valued array `x` with a known mean `m`, optionally
over a dimension `dim`. The weighting vector `wv` specifies frequency weights
(also called case weights) for the result.

This function differs from its counterpart in Base in that Bessel's correction
is not used. That is, here the denominator for the variance is `sum(wv)`,
whereas it's `length(x)-1` in `Base.varm`. The impact is that this is not a
weighted estimate of the population variance based on the sample; it's the weighted
variance of the sample.
"""
Base.varm(v::RealArray, wv::AbstractWeights, m::Real; corrected=false) =
    _moment2(v, wv, m, corrected=corrected)

"""
    var(x, wv::AbstractWeights, [dim]; mean=nothing)

Return the variance of a real-valued array `x`, optionally over a dimension `dim`.
The weighting vector `wv` specifies frequency weights (also called case weights)
for the estimate.

This function differs from its counterpart in Base in that Bessel's correction
is not used. That is, here the denominator for the variance is `sum(wv)`,
whereas it's `length(x)-1` in `Base.var`. The impact is that this is not a
weighted estimate of the population variance based on the sample; it's the weighted
variance of the sample.
"""
function Base.var(v::RealArray, wv::AbstractWeights; mean=nothing, corrected=false)
    if mean == 0
        return varm(v, wv, 0; corrected=corrected)
    elseif mean == nothing
        return varm(v, wv, Base.mean(v, wv); corrected=corrected)
    else
        return varm(v, wv, mean; corrected=corrected)
    end
end

## var along dim

Base.varm!(R::AbstractArray, A::RealArray, wv::AbstractWeights, M::RealArray, dim::Int; corrected=false) =
    scale!(_wsum_centralize!(R, abs2, A, values(wv), M, dim, true), varcorrection(wv, corrected))

function var!(R::AbstractArray, A::RealArray, wv::AbstractWeights, dim::Int; mean=nothing, corrected=false)
    if mean == 0
        Base.varm!(
            R, A, wv, Base.reducedim_initarray(A, dim, 0, eltype(R)), dim;
            corrected=corrected
        )
    elseif mean == nothing
        Base.varm!(
            R, A, wv, Base.mean(A, wv, dim), dim;
            corrected=corrected
        )
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
        Base.varm!(R, A, wv, mean, dim; corrected=corrected)
    end
end

function Base.varm(A::RealArray, wv::AbstractWeights, M::RealArray, dim::Int; corrected=false)
    @static if VERSION < v"0.6.0-dev.1121"
        return Base.varm!(
            similar(A, Float64, Base.reduced_dims(size(A), dim)),
            A, wv, M, dim; corrected=corrected
        )
    else
        return Base.varm!(
            similar(A, Float64, Base.reduced_indices(indices(A), dim)),
            A, wv, M, dim; corrected=corrected
        )
    end
end

function Base.var(A::RealArray, wv::AbstractWeights, dim::Int; mean=nothing, corrected=false)
    @static if VERSION < v"0.6.0-dev.1121"
        return var!(similar(A, Float64, Base.reduced_dims(size(A), dim)), A, wv, dim;
            mean=mean, corrected=corrected)
    else
        return var!(similar(A, Float64, Base.reduced_indices(indices(A), dim)), A, wv, dim;
            mean=mean, corrected=corrected)
    end
end

## std
"""
    stdm(v, wv::AbstractWeights, m, [dim])

Return the standard deviation of a real-valued array `v` with a known mean `m`,
optionally over a dimension `dim`. The weighting vector `wv` specifies frequency
weights (also called case weights) for the estimate.
"""
Base.stdm(v::RealArray, wv::AbstractWeights, m::Real; corrected=false) =
    sqrt(varm(v, wv, m; corrected=corrected))

"""
    std(v, wv::AbstractWeights, [dim]; mean=nothing)

Return the standard deviation of a real-valued array `v`, optionally over a
dimension `dim`. The weighting vector `wv` specifies frequency weights (also
called case weights) for the estimate.
"""
Base.std(v::RealArray, wv::AbstractWeights; mean=nothing, corrected=false) =
    sqrt.(var(v, wv; mean=mean, corrected=corrected))

Base.stdm(v::RealArray, m::RealArray, dim::Int; corrected=false) =
    Base.sqrt!(varm(v, m, dim; corrected=corrected))

Base.stdm(v::RealArray, wv::AbstractWeights, m::RealArray, dim::Int; corrected=false) =
    sqrt.(varm(v, wv, m, dim; corrected=corrected))

Base.std(v::RealArray, wv::AbstractWeights, dim::Int; mean=nothing, corrected=false) =
    sqrt.(var(v, wv, dim; mean=mean, corrected=corrected))

##### Fused statistics
"""
    mean_and_var(x, [wv::AbstractWeights], [dim]) -> (mean, var)

Return the mean and variance of a real-valued array `x`, optionally over a dimension
`dim`, as a tuple. A weighting vector `wv` can be specified to weight the estimates.
The weights are assumed to be frequency weights, also called case weights.
"""
function mean_and_var(A::RealArray; corrected=false)
    m = mean(A)
    v = varm(A, m; corrected=corrected)
    m, v
end

"""
    mean_and_std(x, [wv::AbstractWeights], [dim]) -> (mean, std)

Return the mean and standard deviation of a real-valued array `x`, optionally
over a dimension `dim`, as a tuple. A weighting vector `wv` can be specified
to weight the estimates. The weights are assumed to be frequency weights, also
called case weights.
"""
function mean_and_std(A::RealArray; corrected=false)
    m = mean(A)
    s = stdm(A, m; corrected=corrected)
    m, s
end

function mean_and_var(A::RealArray, wv::AbstractWeights; corrected=false)
    m = mean(A, wv)
    v = varm(A, wv, m; corrected=corrected)
    m, v
end

function mean_and_std(A::RealArray, wv::AbstractWeights; corrected=false)
    m = mean(A, wv)
    s = stdm(A, wv, m; corrected=corrected)
    m, s
end

function mean_and_var(A::RealArray, dim::Int; corrected=false)
    m = mean(A, dim)
    v = varm(A, m, dim; corrected=corrected)
    m, v
end

function mean_and_std(A::RealArray, dim::Int; corrected=false)
    m = mean(A, dim)
    s = stdm(A, m, dim; corrected=corrected)
    m, s
end

function mean_and_var(A::RealArray, wv::AbstractWeights, dim::Int; corrected=false)
    m = mean(A, wv, dim)
    v = varm(A, wv, m, dim; corrected=corrected)
    m, v
end

function mean_and_std(A::RealArray, wv::AbstractWeights, dim::Int; corrected=false)
    m = mean(A, wv, dim)
    s = stdm(A, wv, m, dim; corrected=corrected)
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
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z * z) * w[i]
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
    varcorrection(n, false) * s
end

function _moment3(v::RealArray, wv::AbstractWeights, m::Real)
    n = length(v)
    s = 0.0
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z * z * z) * w[i]
    end
    varcorrection(wv, false) * s
end

function _moment4(v::RealArray, m::Real; corrected=false)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        s += abs2(z * z)
    end
    varcorrection(n, false) * s
end

function _moment4(v::RealArray, wv::AbstractWeights, m::Real)
    n = length(v)
    s = 0.0
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += abs2(z * z) * w[i]
    end
    varcorrection(wv, false) * s
end

function _momentk(v::RealArray, k::Int, m::Real)
    n = length(v)
    s = 0.0
    for i = 1:n
        @inbounds z = v[i] - m
        s += (z ^ k)
    end
    varcorrection(n, false) * s
end

function _momentk(v::RealArray, k::Int, wv::AbstractWeights, m::Real)
    n = length(v)
    s = 0.0
    w = values(wv)
    for i = 1:n
        @inbounds z = v[i] - m
        @inbounds s += (z ^ k) * w[i]
    end
    varcorrection(wv, false) * s
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

moment(v::RealArray, k::Int; corrected=true) = moment(v, k, mean(v))
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
    cf = varcorrection(n, false)
    cm3 *= cf
    cm2 *= cf
    return cm3 / sqrt(cm2 * cm2 * cm2)  # this is much faster than cm2^1.5
end

function skewness(v::RealArray, wv::AbstractWeights, m::Real)
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
    cf = varcorrection(wv, false)
    cm3 *= cf
    cm2 *= cf
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
    cf = varcorrection(n, false)
    cm4 *= cf
    cm2 *= cf
    return (cm4 / (cm2 * cm2)) - 3.0
end

function kurtosis(v::RealArray, wv::AbstractWeights, m::Real)
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
    cf = varcorrection(wv, false)
    cm4 *= cf
    cm2 *= cf
    return (cm4 / (cm2 * cm2)) - 3.0
end

kurtosis(v::RealArray) = kurtosis(v, mean(v))
kurtosis(v::RealArray, wv::AbstractWeights) = kurtosis(v, wv, mean(v, wv))
