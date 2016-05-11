##### Weighted var & std

## var

function Base.varzm(v::RealArray, wv::WeightVec) 
    n = length(v)
    length(wv) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    w = values(wv)
    s = 0.0
    for i = 1:n
        @inbounds s += w[i] * abs2(v[i])
    end
    return s / sum(wv)
end

Base.varm(v::RealArray, wv::WeightVec, m::Real) = _moment2(v, wv, m)

function Base.var(v::RealArray, wv::WeightVec; mean=nothing) 
    mean == 0 ? Base.varzm(v, wv) :
    mean == nothing ? varm(v, wv, Base.mean(v, wv)) :
    varm(v, wv, mean)
end

## var along dim

Base.varzm!(R::AbstractArray, A::RealArray, wv::WeightVec, dim::Int) =
    scale!(_wsum_general!(R, @functorize(abs2), A, values(wv), dim, true), inv(sum(wv)))

Base.varm!(R::AbstractArray, A::RealArray, wv::WeightVec, M::RealArray, dim::Int) =
    scale!(_wsum_centralize!(R, @functorize(abs2), A, values(wv), M, dim, true), inv(sum(wv)))

function var!(R::AbstractArray, A::RealArray, wv::WeightVec, dim::Int; mean=nothing)
    if mean == 0
        Base.varzm!(R, A, wv, dim)
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
    Base.varm!(Array(Float64, Base.reduced_dims(A, dim)), A, wv, M, dim)

Base.var(A::RealArray, wv::WeightVec, dim::Int; mean=nothing) = 
    var!(Array(Float64, Base.reduced_dims(A, dim)), A, wv, dim; mean=mean)

## std

Base.stdm(v::RealArray, wv::WeightVec, m::Real) = sqrt(varm(v, wv, m))
Base.std(v::RealArray, wv::WeightVec; mean=nothing) = sqrt(var(v, wv; mean=mean))

Base.stdm(v::RealArray, m::RealArray, dim::Int) = Base.sqrt!(varm(v, m, dim))
Base.stdm(v::RealArray, wv::WeightVec, m::RealArray, dim::Int) = sqrt(varm(v, wv, m, dim))
Base.std(v::RealArray, wv::WeightVec, dim::Int; mean=nothing) = sqrt(var(v, wv, dim; mean=mean))


##### Fused statistics

mean_and_var(A::RealArray) = (m = mean(A); (m, varm(A, m)))
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

