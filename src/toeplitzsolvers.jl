# Symmetric Toeplitz solver
function durbin!(r::AbstractVector{T}, y::AbstractVector{T}) where T<:BlasReal
    n = length(r)
    n <= length(y) || throw(DimensionMismatch("Auxiliary vector cannot be shorter than data vector"))
    y[1] = -r[1]
    β = one(T)
    α = -r[1]
    for k = 1:n-1
        β *= one(T) - α*α
        α = -r[k+1]
        for j = 1:k
            α -= r[k-j+1]*y[j]
        end
        α /= β
        for j = 1:div(k,2)
            tmp = y[j]
            y[j] += α*y[k-j+1]
            y[k-j+1] += α*tmp
        end
        if isodd(k) y[div(k,2)+1] *= one(T) + α end
        y[k+1] = α
    end
    return y
end
durbin(r::AbstractVector{T}) where {T<:BlasReal} = durbin!(r, zeros(T, length(r)))

function levinson!(r::AbstractVector{T}, b::AbstractVector{T}, x::AbstractVector{T}) where T<:BlasReal
    n = length(b)
    n == length(r) || throw(DimensionMismatch("Vectors must have same length"))
    n <= length(x) || throw(DimensionMismatch("Auxiliary vector cannot be shorter than data vector"))
    x[1] = b[1]
    b[1] = -r[2]/r[1]
    β = one(T)
    α = -r[2]/r[1]
    for k = 1:n-1
        β *= one(T) - α*α
        μ = b[k+1]
        for j = 2:k+1
            μ -= r[j]/r[1]*x[k-j+2]
        end
        μ /= β
        for j = 1:k
            x[j] += μ*b[k-j+1]
        end
        x[k+1] = μ
        if k < n - 1
            α = -r[k+2]
            for j = 2:k+1
                α -= r[j]*b[k-j+2]
            end
            α /= β*r[1]
            for j = 1:div(k,2)
                tmp = b[j]
                b[j] += α*b[k-j+1]
                b[k-j+1] += α*tmp
            end
            if isodd(k) b[div(k,2)+1] *= one(T) + α end
            b[k+1] = α
        end
    end
    for i = 1:n
        x[i] /= r[1]
    end
    return x
end
levinson(r::AbstractVector{T}, b::AbstractVector{T}) where {T<:BlasReal} = levinson!(r, copy(b), zeros(T, length(b)))
