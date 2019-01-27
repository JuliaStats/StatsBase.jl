# Computing deviation in a variety of ways

## count the number of equal/non-equal pairs

"""
    counteq(a, b)

Count the number of indices at which the elements of the arrays
`a` and `b` are equal.
"""
function counteq(a::AbstractArray, b::AbstractArray)
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent lengths."))
    c = 0
    for i = 1:n
        @inbounds if a[i] == b[i]
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
    length(b) == n || throw(DimensionMismatch("Inconsistent lengths."))
    c = 0
    for i = 1:n
        @inbounds if a[i] != b[i]
            c += 1
        end
    end
    return c
end


"""
    sqL2dist(a, b)

Compute the squared L2 distance between `a` and `b`: ``\\sum_{i=1}^n |a_i - b_i|^2``.
"""
sqL2dist(a, b) = sum(abs2.(a .- b))

# L2 distance
"""
    L2dist(a, b)

Compute the L2 distance between `a` and b`: ``\\sqrt{\\sum_{i=1}^n |a_i - b_i|^2}``.
"""
L2dist(a, b) = sqrt(sqL2dist(a, b))


# L1 distance
"""
    L1dist(a, b)

Compute the L1 distance between two arrays: ``\\sum_{i=1}^n |a_i - b_i|``.
"""
L1dist(a, b) = sum(abs.(a .- b))


# Linf distance
"""
    Linfdist(a, b)

Compute the L∞ distance, also called the Chebyshev distance, between
two arrays: ``\\max_{i\\in1:n} |a_i - b_i|``.
Efficient equivalent of `maxabs(a - b)`.
"""
function Linfdist(a, b) = maximum(abs.(a .- b))

# Generalized KL-divergence
"""
    gkldiv(a, b)

Compute the generalized Kullback-Leibler divergence between two arrays:
``\\sum_{i=1}^n (a_i \\log(a_i/b_i) - a_i + b_i)``.
Efficient equivalent of `sum(a*log(a/b)-a+b)`.
"""
function gkldiv(a::AbstractArray{T}, b::AbstractArray{T}) where T<:AbstractFloat
    n = length(a)
    r = 0.0
    for i = 1:n
        @inbounds ai = a[i]
        @inbounds bi = b[i]
        if ai > 0
            r += (ai * log(ai / bi) - ai + bi)
        else
            r += bi
        end
    end
    return r::Float64
end


# MeanAD: mean absolute deviation
"""
    meanad(a, b)

Return the mean absolute deviation between two arrays: `mean(abs(a - b))`.
"""
meanad(a, b) = mean(abs.(a .- b))

# MaxAD: maximum absolute deviation
"""
    maxad(a, b)

Return the maximum absolute deviation between `a` and `b`.
Equivelent to `Linfdist`
"""
maxad(a, b) = Linfdist(a, b)


# MSD: mean squared deviation
"""
    msd(a, b)

Return the mean squared deviation between two arrays: `mean(abs2(a - b))`.
"""
msd(a, b) = mean(abs2.(a .- b))

# RMSD: root mean squared deviation
"""
    rmsd(a, b; normalize=false)

Return the root mean squared deviation between two optionally
normalized arrays. The root mean squared deviation is computed
as `sqrt(msd(a, b))`.
"""
function rmsd(a::AbstractArray{T}, b::AbstractArray{T}; normalize::Bool=false) where T<:Number
    v = sqrt(msd(a, b))
    if normalize
        amin, amax = extrema(a)
        v /= (amax - amin)
    end
    return v
end


# PSNR: peak signal-to-noise ratio
"""
    psnr(a, b, maxv)

Compute the peak signal-to-noise ratio between `a` and `b`.
`maxv` is the maximum possible value either array can take. The PSNR
is computed as `10 * log10(maxv^2 / msd(a, b))`.
"""
psnr(a, b, maxv) = 20log10(maxv) - 10log10(msd(a, b))
