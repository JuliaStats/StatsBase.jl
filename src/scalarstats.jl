# Descriptive Statistics


#############################
#
#   Location
#
#############################

# Geometric mean
function geomean(a::RealArray)
    s = 0.0
    n = length(a)
    for i = 1 : n
        @inbounds s += log(a[i])
    end
    return exp(s / n)
end

# Harmonic mean
function harmmean(a::RealArray)
    s = 0.0
    n = length(a)
    for i in 1 : n
        @inbounds s += inv(a[i])
    end
    return n / s
end

# Trimmed mean
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


#############################
#
#   quantile and friends
#
#############################

percentile{T<:Real}(v::AbstractArray{T}, p) = quantile(v, p * 0.01)

quantile{T<:Real}(v::AbstractArray{T}) = quantile(v, [.0, .25, .5, .75, 1.0])
nquantile{T<:Real}(v::AbstractArray{T}, n::Integer) = quantile(v, (0:n)/n)


#############################
#
#   Dispersion
#
#############################

# span, i.e. the range minimum(x):maximum(x)
span{T<:Integer}(x::AbstractArray{T}) = ((a, b) = extrema(x); a:b)

# Variation coefficient: std / mean
variation{T<:Real}(x::AbstractArray{T}, m::Real) = stdm(x, m) / m
variation{T<:Real}(x::AbstractArray{T}) = variation(x, mean(x))

# Standard error of the mean: std(a) / sqrt(len)
sem{T<:Real}(a::AbstractArray{T}) = sqrt(var(a) / length(a))

# Median absolute deviation
mad{T<:Real}(v::AbstractArray{T}, args...;arg...) = mad!(copy(v), args...;arg...)
mad{T<:Real}(v::Range{T}, args...;arg...) = mad!([v;], args...;arg...)

function mad!{T<:Real}(v::AbstractArray{T}, center::Real=median!(v); constant::Real=1.4826)
    for i in 1:length(v)
        @inbounds v[i] = abs(v[i]-center)
    end
    constant * median!(v)
end

# Interquartile range
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

#############################
#
#   summary
#
#############################
immutable SummaryStats{T<:AbstractFloat} 
    mean::T
    std::T
    min::T
    p25::T
    median::T
    p75::T
    max::T
end

function describe{T<:Real}(a::AbstractArray{T}, ::Type{Val{false}})
    m = mean(a)
    std = stdm(a, m)
    R = typeof(zero(T)/1)
    p = [0.00, 0.25, 0.50, 0.75, 1.00]
    qs = quantile(a, p)
    return SummaryStats{R}(
        convert(R, m),
        convert(R, std),
        convert(Vector{R}, qs)...
        )
end

immutable DetailedSummaryStats{T<:AbstractFloat}  
    mean::T
    std::T
    skewness::T
    kurtosis::T
    min::T
    p1::T
    p5::T
    p10::T
    p25::T
    median::T
    p75::T
    p90::T
    p95::T
    p99::T
    max::T
end

function describe{T<:Real}(a::AbstractArray{T}, ::Type{Val{true}})
    m = mean(a)
    std = stdm(a, m)
    R = typeof(zero(T)/1)
    sk = skewness(a, m)
    ku = kurtosis(a, m)
    p = [0.00, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.00]
    qs = quantile(a, p)
    return DetailedSummaryStats{R}(
        convert(R, m), 
        convert(R, std),
        convert(R, sk),
        convert(R, ku),
        convert(Vector{R}, qs)...
        )
end

# Weighted
immutable WeightedSummaryStats{T<:AbstractFloat, W<:Real} 
    mean::T
    std::T
    min::T
    p25::T
    median::T
    p75::T
    max::T 
    wsum::W
end

function describe{T<:Real, W<:Real}(a::AbstractVector{T}, w::WeightVec{W}, detail::Type{Val{false}})
    m = mean(a, w)
    std = stdm(a, m, w)
    R = typeof(zero(T)/1)
    p = [0.00, 0.25, 0.50, 0.75, 1.00]
    qs = quantile(a, w, p)
    return WeightedSummaryStats{R, W}(
        convert(R, m),
        convert(R, std),
        convert(Vector{R}, qs)...,
        w.sum
        )
end

immutable WeightedDetailedSummaryStats{T<:AbstractFloat, W<:Real} 
    mean::T
    std::T
    skewness::T
    kurtosis::T
    min::T
    p1::T
    p5::T
    p10::T
    p25::T
    median::T
    p75::T
    p90::T
    p95::T
    p99::T
    max::T
    wsum::W
end

function describe{T<:Real, W<:Real}(a::AbstractVector{T}, w::WeightVec{W}, detail::Type{Val{true}})
    m = mean(a, w)
    std = stdm(a, m, w)
    R = typeof(zero(T)/1)
    sk = skewness(a, w, m)
    ku = kurtosis(a, w, m)
    p = [0.00, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.00]
    qs = quantile(a, w, p)
    return WeightedDetailedSummaryStats{R, W}(
        convert(R, m),
        convert(R, std),
        convert(R, sk),
        convert(R, ku),
        convert(Vector{R}, qs)...,
        w.sum
        )  
end  

# value vs type for detail
function describe{T<:Real}(a::AbstractArray{T}; detail = false) 
    return describe(a, Val{detail})
end

function describe{T<:Real, W<:Real}(a::AbstractVector{T}, w::WeightVec{W}; detail = false) 
    return describe(a, w, Val{detail})
end

# show methods
function Base.show(io::IO, ss::@compat(Union{SummaryStats, WeightedSummaryStats}))
    println(io, "Summary Stats:")
    @printf(io, "Mean:           %.6f\n", ss.mean)
    @printf(io, "Std:            %.6f\n", ss.std)
    @printf(io, "Minimum:        %.6f\n", ss.min)
    @printf(io, "1st Quartile:   %.6f\n", ss.p25)
    @printf(io, "Median:         %.6f\n", ss.median)
    @printf(io, "3rd Quartile:   %.6f\n", ss.p75)
    @printf(io, "Maximum:        %.6f\n", ss.max)
    show_sumw(io, ss)
end

function Base.show(io::IO, ss::@compat(Union{DetailedSummaryStats, WeightedDetailedSummaryStats}))
    println(io, "Detailed Summary Stats:")
    @printf(io, "Mean:           %.6f\n", ss.mean)
    @printf(io, "Std:            %.6f\n", ss.std)
    @printf(io, "Skewness:       %.6f\n", ss.skewness)
    @printf(io, "Kurtosis:       %.6f\n", ss.kurtosis)
    print("\n")  
    @printf(io, "Minimum:        %.6f\n", ss.min)
    @printf(io, "p1:             %.6f\n", ss.p1)
    @printf(io, "p5:             %.6f\n", ss.p5)
    @printf(io, "p10:            %.6f\n", ss.p10)
    @printf(io, "p25:            %.6f\n", ss.p25)
    @printf(io, "Median:         %.6f\n", ss.median)
    @printf(io, "p75             %.6f\n", ss.p75)
    @printf(io, "p90:            %.6f\n", ss.p90)
    @printf(io, "p95:            %.6f\n", ss.p95)
    @printf(io, "p99:            %.6f\n", ss.p99)
    @printf(io, "Maximum:        %.6f\n", ss.max)
    show_sumw(io, ss)
end

show_sumw(io::IO, ss::@compat(Union{SummaryStats, DetailedSummaryStats})) = nothing
function show_sumw(io::IO, ss::@compat(Union{WeightedSummaryStats, WeightedDetailedSummaryStats}))
    print("\n")  
    @printf(io, "Sum of weights: %.6f\n", ss.wsum)
end