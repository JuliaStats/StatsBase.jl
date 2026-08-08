"""
    durbin!(r::AbstractVector, y::AbstractVector, p::AbstractVector) -> Vector

Solve Yule-Walker equations using Durbin algorithm.

Solution of N×N system is obtained iteratively, by solving 1×1, 2×2, ...
in succession. For use in computing partial autocorrelation matrix,
return the last elements of the successive solutions in vector `p`.

Reference: Golub, G. H., and C. F. Van Loan. "Matrix computations 4th edition the johns hopkins university press." Baltimore, MD (2013), section 4.7

# Arguments
- `r::AbstractVector`: first column of a Herimitian positive definite
    Toeplitz matrix, excluding the diagonal element (equal to one).
- `y::AbstractVector`: work vector for solution, should have the same length
    as `r`
- `p::AbstractVector`: last elements of the successive solutions, should
     have the same length  as `r`

# Returns
- `y::AbstractVector`: return the solution vector (same as the second argument)
"""
function durbin!(r::AbstractVector{T}, y::AbstractVector{T}, p::AbstractVector{T}) where {T}
    n = length(r)
    n <= length(p) || n <= length(y) || throw(DimensionMismatch("Auxiliary vectors cannot be shorter than data vector"))
    y[1] = -r[1]
    p[1] = -r[1]
    β = one(T)
    α = -r[1]
    for k = 1:n-1
        β *= one(T) - abs2(α)
        α = -r[k+1]
        for j = 1:k
            α -= r[k-j+1]*y[j]
        end
        α /= β
        for j = 1:div(k,2)
            tmp = y[j]
            y[j] += α*conj(y[k-j+1])
            y[k-j+1] += α*conj(tmp)
        end
        if isodd(k)
            y[div(k,2)+1] += α*conj(y[div(k,2)+1]) 
        end
        y[k+1] = α
        p[k+1] = α
    end
    return y
end
durbin(r::AbstractVector{T}) where {T} = durbin!(r, zeros(T, length(r)), zeros(T, length(r)))

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
