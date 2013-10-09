# Statistics for an array of scalars


# Variation: std / mean
variation{T<:Real}(x::AbstractArray{T}, m::Real) = stdm(x, m) / m
variation{T<:Real}(x::AbstractArray{T}) = variation(x, mean(x))


# Standard error of the mean: std(a)
sem{T<:Real}(a::AbstractArray{T}) = sqrt(var(a) / length(a))


# Median absolute deviation
mad{T<:Real}(v::AbstractArray{T}, center::Real) = 1.4826 * median!(abs(v-center))

function mad{T<:Real}(v::AbstractArray{T})
    v = copy(v)
    center = median!(v)
    for i in 1:length(v)
        v[i] = abs(v[i]-center)
    end
    1.4826 * median!(v, checknan=false)
end


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
    return cm3 / (cm2^1.5)
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
    return (cm4 / (cm2^2)) - 3.0
end

kurtosis{T<:Real}(v::AbstractVector{T}) = kurtosis(v, mean(v))


# Minimum and maximum
function minmax{T<:Real}(x::AbstractArray{T})
    if isempty(x)
        error("minmax: x cannot be empty.")
    end

    xmin = xmax = x[1]
    for i = 2:length(x)
        x_i = x[i]
        if x_i < xmin
            xmin = x_i
        elseif x_i > xmax
            xmax = x_i
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


# order (aka. rank), resolving ties using the mean rank
function tiedrank{T<:Real}(v::AbstractArray{T})
    n     = length(v)
    place = sortperm(v)
    ord   = Array(Float64, n)

    i = 1
    while i <= n
        j = i
        while j + 1 <= n && v[place[i]] == v[place[j + 1]]
            j += 1
        end

        if j > i
            m = sum(i:j) / (j - i + 1)
            for k = i:j
                ord[place[k]] = m
            end
        else
            ord[place[i]] = i
        end

        i = j + 1
    end

    return ord
end

tiedrank{T<:Real}(X::AbstractMatrix{T}) = tiedrank(reshape(X, length(X)))

# quantile and friends

quantile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.0, .25, .5, .75, 1.0])
percentile{T<:Real}(v::AbstractVector{T}) = quantile(v, [1:99] / 100)
quartile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.25, .5, .75])
quintile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.2, .4, .6, .8])
decile{T<:Real}(v::AbstractVector{T}) = quantile(v, [.1, .2, .3, .4, .5, .6, .7, .8, .9])
iqr{T<:Real}(v::AbstractVector{T}) = quantile(v, [.25, .75])


# table & mode

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
    m = max(values(tab))
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
    m = max(values(tab))
    for (k, v) in tab
        if v == m
            push!(res, k)
        end
    end
    return res
end


# Print out basic summary statistics
function describe{T<:Real}(a::AbstractArray{T})
    q00, q25, q50, q75, q10 = quantile(a, [0.00, 0.25, 0.50, 0.75, 1.00])
    @printf "Min:          %.6f\n" q00
    @printf "1st Quartile: %.6f\n" q25
    @printf "Median:       %.6f\n" q50
    @printf "Mean:         %.6f\n" mean(a)
    @printf "3rd Quartile: %.6f\n" q75
    @printf "Max:          %.6f\n" q10
    return
end

