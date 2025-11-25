# Computing deviation in a variety of ways

## count the number of equal/non-equal pairs

"""
    counteq(a, b)

Count the number of indices at which the elements of the arrays
`a` and `b` are equal.
"""
function counteq(a::AbstractArray, b::AbstractArray)
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    c = 0
    for i in eachindex(a, b)
        if a[i] == b[i]
            c += 1
        end
    end
    return c
end


"""
    countne(a, b)

Count the number of indices at which the elements of the arrays
`a` and `b` are not equal.
"""
function countne(a::AbstractArray, b::AbstractArray)
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    c = 0
    for i in eachindex(a, b)
        if a[i] != b[i]
            c += 1
        end
    end
    return c
end


"""
    sqL2dist(a, b)

Compute the squared L2 distance between two arrays: ``\\sum_{i=1}^n |a_i - b_i|^2``.
Efficient equivalent of `sum(abs2, a - b)`.
"""
function sqL2dist(a::AbstractArray{<:Number}, b::AbstractArray{<:Number})
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    if iszero(n)
        r = zero(abs2(zero(eltype(a)) - zero(eltype(b))))
    else
        broadcasted = Broadcast.broadcasted((ai, bi) -> abs2(ai - bi), vec(a), vec(b))
        r = sum(Broadcast.instantiate(broadcasted))
    end
    return r
end


# L2 distance
"""
    L2dist(a, b)

Compute the L2 distance between two arrays: ``\\sqrt{\\sum_{i=1}^n |a_i - b_i|^2}``.
Efficient equivalent of `sqrt(sum(abs2, a - b))`.
"""
L2dist(a::AbstractArray{<:Number}, b::AbstractArray{<:Number}) = sqrt(sqL2dist(a, b))


# L1 distance
"""
    L1dist(a, b)

Compute the L1 distance between two arrays: ``\\sum_{i=1}^n |a_i - b_i|``.
Efficient equivalent of `sum(abs, a - b)`.
"""
function L1dist(a::AbstractArray{<:Number}, b::AbstractArray{<:Number})
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    if iszero(n)
        r = zero(abs(zero(eltype(a)) - zero(eltype(b))))
    else
        broadcasted = Broadcast.broadcasted((ai, bi) -> abs(ai - bi), vec(a), vec(b))
        r = sum(Broadcast.instantiate(broadcasted))
    end
    return r
end


# Linf distance
"""
    Linfdist(a, b)

Compute the L∞ distance, also called the Chebyshev distance, between
two arrays: ``\\max_{1≤i≤n} |a_i - b_i|``.
Efficient equivalent of `maxabs(a - b)`.
"""
function Linfdist(a::AbstractArray{<:Number}, b::AbstractArray{<:Number})
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    if iszero(n)
        r = zero(abs(zero(eltype(a)) - zero(eltype(b))))
    else
        broadcasted = Broadcast.broadcasted((ai, bi) -> abs(ai - bi), vec(a), vec(b))
        r = maximum(Broadcast.instantiate(broadcasted))
    end
    return r
end


# Generalized KL-divergence
"""
    gkldiv(a, b)

Compute the generalized Kullback-Leibler divergence between two arrays:
``\\sum_{i=1}^n (a_i \\log(a_i/b_i) - a_i + b_i)``.
Efficient equivalent of `sum(a*log(a/b)-a+b)`.
"""
function gkldiv(a::AbstractArray{<:Real}, b::AbstractArray{<:Real})
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    if iszero(n)
        za = zero(eltype(a))
        zb = zero(eltype(b))
        r = zero(xlogy(za, za / zb) + (zb - za))
    else
        broadcasted = Broadcast.broadcasted(vec(a), vec(b)) do ai, bi
            return xlogy(ai, ai / bi) + (bi - ai)
        end
        return sum(Broadcast.instantiate(broadcasted))
    end
    return r
end


# MeanAD: mean absolute deviation
"""
    meanad(a, b)

Return the mean absolute deviation between two arrays: `mean(abs, a - b)`.
"""
meanad(a::AbstractArray{<:Number}, b::AbstractArray{<:Number}) = L1dist(a, b) / length(a)


# MaxAD: maximum absolute deviation
"""
    maxad(a, b)

Return the maximum absolute deviation between two arrays: `maxabs(a - b)`.
"""
maxad(a::AbstractArray{<:Number}, b::AbstractArray{<:Number}) = Linfdist(a, b)


# MSD: mean squared deviation
"""
    msd(a, b)

Return the mean squared deviation between two arrays: `mean(abs2, a - b)`.
"""
msd(a::AbstractArray{<:Number}, b::AbstractArray{<:Number}) = sqL2dist(a, b) / length(a)


# RMSD: root mean squared deviation
"""
    rmsd(a, b; normalize=false)

Return the root mean squared deviation between two optionally
normalized arrays. The root mean squared deviation is computed
as `sqrt(msd(a, b))`.
"""
function rmsd(a::AbstractArray{<:Number}, b::AbstractArray{<:Number}; normalize::Bool=false)
    v = sqrt(msd(a, b))
    if normalize
        amin, amax = isempty(a) ? (zero(eltype(a)), zero(eltype(a))) : extrema(a)
        return v / (amax - amin)
    else
        return v
    end
end


# PSNR: peak signal-to-noise ratio
"""
    psnr(a, b, maxv)

Compute the peak signal-to-noise ratio between two arrays `a` and `b`.
`maxv` is the maximum possible value either array can take. The PSNR
is computed as `10 * log10(maxv^2 / msd(a, b))`.
"""
function psnr(a::AbstractArray{<:Real}, b::AbstractArray{<:Real}, maxv::Real)
    msd_a_b, _maxv = promote(msd(a, b), maxv)
    return 20 * log10(_maxv) - 10 * log10(msd_a_b)
end
