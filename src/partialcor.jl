# Partial correlation

"""
    partialcor(x, y, Z)

Compute the partial correlation of the vectors `x` and `y` given `Z`, which can be
a vector or matrix.
"""
function partialcor(x::AbstractVector, y::AbstractVector, Z::AbstractVecOrMat)
    length(x) == length(y) == size(Z, 1) ||
        throw(DimensionMismatch("Inputs must have the same number of observations"))
    length(x) > 0 || throw(ArgumentError("Inputs must be non-empty"))
    return Statistics.clampcor(_partialcor(x, mean(x), y, mean(y), Z))
end

function _partialcor(x::AbstractVector, μx, y::AbstractVector, μy, Z::AbstractMatrix)
    p = size(Z, 2)
    p == 1 && return _partialcor(x, μx, y, μy, vec(Z))
    z₀   = view(Z, :, 1)
    Zmz₀ = view(Z, :, 2:p)
    μz₀ = mean(z₀)
    rxz = _partialcor(x,  μx,  z₀, μz₀, Zmz₀)
    rzy = _partialcor(z₀, μz₀, y,  μy,  Zmz₀)
    rxy = _partialcor(x,  μx,  y,  μy,  Zmz₀)::typeof(rxz)
    return (rxy - rxz * rzy) / (sqrt(1 - rxz^2) * sqrt(1 - rzy^2))
end

function _partialcor(x::AbstractVector, μx, y::AbstractVector, μy, z::AbstractVector)
    μz = mean(z)

    # Initialize all of the accumulators to 0 of the appropriate types
    Σxx = abs2(zero(eltype(x)) - zero(μx))
    Σyy = abs2(zero(eltype(y)) - zero(μy))
    Σzz = abs2(zero(eltype(z)) - zero(μz))
    Σxy = zero(Σxx * Σyy)
    Σxz = zero(Σxx * Σzz)
    Σzy = zero(Σzz * Σyy)

    # We only want to make one pass over all of the arrays
    @inbounds begin
        @simd for i in eachindex(x, y, z)
            xi = x[i] - μx
            yi = y[i] - μy
            zi = z[i] - μz

            Σxx += abs2(xi)
            Σyy += abs2(yi)
            Σzz += abs2(zi)

            Σxy += xi * yi
            Σxz += xi * zi
            Σzy += zi * yi
        end
    end

    # Individual pairwise correlations
    rxy = Σxy / sqrt(Σxx * Σyy)
    rxz = Σxz / sqrt(Σxx * Σzz)
    rzy = Σzy / sqrt(Σzz * Σyy)

    return (rxy - rxz * rzy) / (sqrt(1 - rxz^2) * sqrt(1 - rzy^2))
end
