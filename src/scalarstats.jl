# Statistics for an array of scalars


#############################
#
#   Moments
#
#############################

# Skewness
# This is Type 1 definition according to Joanes and Gill (1998)
function skewness{T<:Real}(v::AbstractVector{T}, m::Real)
    n = length(v)
    cm2 = 0.0   # empirical 2nd centered moment (variance)
    cm3 = 0.0   # empirical 3rd centered moment
    for x_i in v
        z = x_i - m
        z2 = z * z

        cm2 += z2
        cm3 += z2 * z
    end
    cm3 /= n
    cm2 /= n
    return cm3 / sqrt(cm2 * cm2 * cm2)  # this is much faster than cm2^1.5
end

skewness{T<:Real}(v::AbstractVector{T}) = skewness(v, mean(v))

# (excessive) Kurtosis
# This is Type 1 definition according to Joanes and Gill (1998)
function kurtosis{T<:Real}(v::AbstractVector{T}, m::Real)
    n = length(v)
    cm2 = 0.0  # empirical 2nd centered moment (variance)
    cm4 = 0.0  # empirical 4th centered moment
    for x_i in v
        z = x_i - m
        z2 = z * z
        cm2 += z2
        cm4 += z2 * z2
    end
    cm4 /= n
    cm2 /= n
    return (cm4 / (cm2 * cm2)) - 3.0
end

kurtosis{T<:Real}(v::AbstractVector{T}) = kurtosis(v, mean(v))


#############################
#
#   Variability measurements
#
#############################

# Variation: std / mean
variation{T<:Real}(x::AbstractArray{T}, m::Real) = stdm(x, m) / m
variation{T<:Real}(x::AbstractArray{T}) = variation(x, mean(x))

# Standard error of the mean: std(a)
sem{T<:Real}(a::AbstractArray{T}) = sqrt(var(a) / length(a))

# Median absolute deviation
mad{T<:Real}(v::AbstractArray{T}, center::Real) = 1.4826 * median!(abs(v-center))

function mad!{T<:Real}(v::AbstractArray{T}, center::Real)
    for i in 1:length(v)
        v[i] = abs(v[i]-center)
    end
    1.4826 * median!(v, checknan=false)
end

mad!{T<:Real}(v::AbstractArray{T}) = mad!(v, median!(v))

mad{T<:Real}(v::AbstractArray{T}) = mad!(copy(v))
mad{T<:Real}(v::Range1{T}) = mad!([v])


#############################
#
#   min/max related
#
#############################

# Minimum and maximum
function minmax{T<:Real}(x::AbstractArray{T})
    isempty(x) && error("minmax: x cannot be empty.")

    # find the first non-NaN value
    x0 = x[1]
    n = length(x)
    i = 1
    while i < n && (x0 != x0)
        i += 1
        @inbounds x0 = x[i] 
    end
    xmin = xmax = x0

    # deal with the remaining
    while i < n
        i += 1
        @inbounds xi = x[i]
        if xi < xmin
            xmin = xi
        elseif xi > xmax
            xmax = xi
        end
    end

    return xmin, xmax
end

# Mid-range
function midrange{T<:Real}(a::AbstractArray{T})
    xmin, xmax = minmax(a)
    return xmin + (xmax - xmin) / 2
end

# Range: xmax - xmin
function range{T<:Real}(a::AbstractArray{T})
    xmin, xmax = minmax(a)
    return xmax - xmin
end


#############################
#
#   quantile and friends
#
#############################

quantile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.0, .25, .5, .75, 1.0])
percentile{T<:Real}(v::AbstractVector{T}) = quantile(v, [1:99] / 100)
quartile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.25, .5, .75])
quintile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.2, .4, .6, .8])
decile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.1, .2, .3, .4, .5, .6, .7, .8, .9])
iqr{T<:Real}(v::AbstractVector{T}) = quantile(v, [.25, .75])


#############################
#
#   table & mode
#
#############################

function table{T}(a::AbstractArray{T})
    counts = Dict{T, Int}()
    for i = 1:length(a)
        tmp = a[i]
        counts[tmp] = get(counts, tmp, 0) + 1
    end
    return counts
end

function mode{T}(a::AbstractArray{T})
    if isempty(a)
        throw(ArgumentError("mode: a cannot be empty."))
    end
    tab = table(a)
    m = maximum(values(tab))
    for (k, v) in tab
        if v == m
            return k
        end
    end
end

function modes{T}(a::AbstractArray{T})
    if isempty(a)
        throw(ArgumentError("mode: a cannot be empty."))
    end
    res = Array(T, 0)
    tab = table(a)
    m = maximum(values(tab))
    for (k, v) in tab
        if v == m
            push!(res, k)
        end
    end
    return res
end



#############################
#
#   summary
#
#############################

immutable SummaryStats{T<:FloatingPoint}
    mean::T
    min::T
    q25::T    
    median::T    
    q75::T
    max::T
end

function summarystats{T<:Real}(a::AbstractArray{T})
    m = mean(a)
    qs = quantile(a, [0.00, 0.25, 0.50, 0.75, 1.00])    
    R = typeof(convert(FloatingPoint, zero(T)))
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

describe{T<:Real}(a::AbstractArray{T}) = show(summarystats(a))

