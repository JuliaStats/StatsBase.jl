# Rank-based correlations
#
# - Spearman's correlation
# - Kendall's correlation
#

#######################################
#
#   Spearman correlation
#
#######################################

"""
    corspearman(x, y=x)

Compute Spearman's rank correlation coefficient. If `x` and `y` are vectors, the
output is a float, otherwise it's a matrix corresponding to the pairwise correlations
of the columns of `x` and `y`.
"""
corspearman(x::RealVector, y::RealVector) = cor(tiedrank(x), tiedrank(y))

corspearman(X::RealMatrix, Y::RealMatrix) =
    cor(mapslices(tiedrank, X, dims=1), mapslices(tiedrank, Y, dims=1))
corspearman(X::RealMatrix, y::RealVector) = cor(mapslices(tiedrank, X, dims=1), tiedrank(y))
corspearman(x::RealVector, Y::RealMatrix) = cor(tiedrank(x), mapslices(tiedrank, Y, dims=1))

corspearman(X::RealMatrix) = (Z = mapslices(tiedrank, X, dims=1); cor(Z, Z))


#######################################
# 
#   Kendall correlation
# 
#######################################

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833. Accessed 15 Jan. 2021.
function corkendall!(x::RealVector, y::RealVector, permx=sortperm(x))
    if any(isnan, x) || any(isnan, y) return NaN end
    n = length(x)
    if n ≠ length(y) error("Vectors must have same length") end

    # Initial sorting
    x[:] = x[permx]
    y[:] = y[permx]

    npairs = div(n * (n - 1), 2)
    ntiesx, ntiesy, ndoubleties, k, nswaps = 0, 0, 0, 0, 0

    @inbounds for i ∈ 2:n
        if x[i - 1] == x[i]
            k += 1
        elseif k > 0
            # Sort the corresponding chunk of y, so the rows of hcat(x,y) are 
            # sorted first on x, then (where x values are tied) on y. Hence 
            # double ties can be counted by calling countties.
            sort!(view(y, (i - k - 1):(i - 1)))
            ntiesx += div(k * (k + 1), 2)
            ndoubleties += countties(y,  i - k - 1, i - 1)
            k = 0
        end
    end
    if k > 0
        sort!(view(y, ((n - k):n)))
        ntiesx += div(k * (k + 1), 2)
        ndoubleties += countties(y,  n - k, n)
    end

    nswaps = mergesort!(y, 1, n)
    ntiesy = countties(y, 1, n)

    (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
     sqrt((npairs - ntiesx) * (npairs - ntiesy))
end

"""
    countties(x::RealVector,lo::Int64,hi::Int64)

Assumes `x` is sorted. Returns the number of ties within `x[lo:hi]`.
"""
function countties(x::AbstractVector, lo::Integer, hi::Integer)
    thistiecount, result = 0, 0
    (lo < 1 || hi > length(x)) && error("Bounds error")
    @inbounds for i ∈ (lo + 1):hi
        if x[i] == x[i - 1]
            thistiecount += 1
        elseif thistiecount > 0
            result += div(thistiecount * (thistiecount + 1), 2)
            thistiecount = 0
        end
    end

    if thistiecount > 0
        result += div(thistiecount * (thistiecount + 1), 2)
    end
    result
end

"""
    corkendall(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors.
"""
corkendall(x::RealVector, y::RealVector) = corkendall!(float(copy(x)), float(copy(y)))

corkendall(X::RealMatrix, y::RealVector) = (permy = sortperm(y);Float64[corkendall!(float(copy(y)), float(X[:,i]), permy) for i in 1:size(X, 2)])

corkendall(x::RealVector, Y::RealMatrix) = (n = size(Y, 2); permx = sortperm(x); reshape(Float64[corkendall!(float(copy(x)), float(Y[:,i]), permx) for i in 1:n], 1, n))

function corkendall(X::RealMatrix)
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra
    @inbounds for j ∈ 2:n
        permx = sortperm(X[:,j])
        for i ∈ 1:j - 1
            C[j,i] = corkendall!(X[:,j], X[:,i], permx)
            C[i,j] = C[j,i]
        end
    end
    return C
end

function corkendall(X::RealMatrix, Y::RealMatrix)
    nr = size(X, 2)
    nc = size(Y, 2)
    C = zeros(float(eltype(X)), nr, nc)
    @inbounds for j ∈ 1:nr
        permx = sortperm(X[:,j])
        for i ∈ 1:nc
            C[j,i] = corkendall!(X[:,j], Y[:,i], permx)
        end
    end
    return C
end

# Auxilliary functions for Kendall's rank correlation

const SMALL_THRESHOLD  = 64

# Copy was from https://github.com/JuliaLang/julia/commit/28330a2fef4d9d149ba0fd3ffa06347b50067647 dated 20 Sep 2020
"""
    mergesort!(v::AbstractVector, lo::Integer, hi::Integer, t=similar(v, 0))    

Mutates `v` by sorting elements `x[lo:hi]` using the mergesort algorithm. 
This method is a copy-paste-edit of sort! in base/sort.jl (the method specialised on MergeSortAlg),
but amended to return the bubblesort distance.
"""
function mergesort!(v::AbstractVector, lo::Integer, hi::Integer, t=similar(v, 0))
    nswaps = 0
    @inbounds if lo < hi
        hi - lo <= SMALL_THRESHOLD && return insertionsort!(v, lo, hi)

        m = midpoint(lo, hi)
        (length(t) < m - lo + 1) && resize!(t, m - lo + 1)

        nswaps = mergesort!(v, lo,  m, t)
        nswaps += mergesort!(v, m + 1, hi, t)

        i, j = 1, lo
        while j <= m
            t[i] = v[j]
            i += 1
            j += 1
        end

        i, k = 1, lo
        while k < j <= hi
            if v[j] < t[i]
                v[k] = v[j]
                j += 1
                nswaps += m - lo + 1 - (i - 1)
            else
                v[k] = t[i]
                i += 1
            end
            k += 1
        end
        while k < j
            v[k] = t[i]
            k += 1
            i += 1
        end
    end
    return nswaps
end

# This implementation of `midpoint` is performance-optimized but safe only if `lo <= hi`.
# This function is also copied from base/sort.jl
midpoint(lo::T, hi::T) where T <: Integer = lo + ((hi - lo) >>> 0x01)
midpoint(lo::Integer, hi::Integer) = midpoint(promote(lo, hi)...)

# Copy was from https://github.com/JuliaLang/julia/commit/28330a2fef4d9d149ba0fd3ffa06347b50067647 dated 20 Sep 2020
"""
    insertionsort!(v::AbstractVector, lo::Integer, hi::Integer)

Mutates `v` by sorting elements `x[lo:hi]` using the insertionsort algorithm. 
This method is a copy-paste-edit of sort! in base/sort.jl (the method specialised on InsertionSortAlg),
amended to return the bubblesort distance.
"""
function insertionsort!(v::AbstractVector, lo::Integer, hi::Integer)
    if lo == hi return 0 end
    nswaps = 0
    @inbounds for i = lo + 1:hi
        j = i
        x = v[i]
        while j > lo
            if x < v[j - 1]
                nswaps += 1
                v[j] = v[j - 1]
                j -= 1
                continue
            end
            break
        end
        v[j] = x
    end
    return nswaps
end

