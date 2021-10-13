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
function corspearman(x::RealVector, y::RealVector)
    n = length(x)
    n == length(y) || throw(DimensionMismatch("vectors must have same length"))
    (any(isnan, x) || any(isnan, y)) && return NaN
    return cor(tiedrank(x), tiedrank(y))
end

function corspearman(X::RealMatrix, y::RealVector)
    size(X, 1) == length(y) ||
        throw(DimensionMismatch("X and y have inconsistent dimensions"))
    n = size(X, 2)
    C = Matrix{Float64}(I, n, 1)
    any(isnan, y) && return fill!(C, NaN)
    yrank = tiedrank(y)
    for j = 1:n
        Xj = view(X, :, j)
        if any(isnan, Xj)
            C[j,1] = NaN
        else
            Xjrank = tiedrank(Xj)
            C[j,1] = cor(Xjrank, yrank)
        end
    end
    return C
end

function corspearman(x::RealVector, Y::RealMatrix)
    size(Y, 1) == length(x) ||
        throw(DimensionMismatch("x and Y have inconsistent dimensions"))
    n = size(Y, 2)
    C = Matrix{Float64}(I, 1, n)
    any(isnan, x) && return fill!(C, NaN)
    xrank = tiedrank(x)
    for j = 1:n
        Yj = view(Y, :, j)
        if any(isnan, Yj)
            C[1,j] = NaN
        else
            Yjrank = tiedrank(Yj)
            C[1,j] = cor(xrank, Yjrank)
        end
    end
    return C
end

function corspearman(X::RealMatrix)
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)
    anynan = Vector{Bool}(undef, n)
    for j = 1:n
        Xj = view(X, :, j)
        anynan[j] = any(isnan, Xj)
        if anynan[j]
            C[:,j] .= NaN
            C[j,:] .= NaN
            C[j,j] = 1
            continue
        end
        Xjrank = tiedrank(Xj)
        for i = 1:(j-1)
            Xi = view(X, :, i)
            if anynan[i]
                C[i,j] = C[j,i] = NaN
            else
                Xirank = tiedrank(Xi)
                C[i,j] = C[j,i] = cor(Xjrank, Xirank)
            end
        end
    end
    return C
end

function corspearman(X::RealMatrix, Y::RealMatrix)
    size(X, 1) == size(Y, 1) ||
        throw(ArgumentError("number of rows in each array must match"))
    nr = size(X, 2)
    nc = size(Y, 2)
    C = Matrix{Float64}(undef, nr, nc)
    for j = 1:nr
        Xj = view(X, :, j)
        if any(isnan, Xj)
            C[j,:] .= NaN
            continue
        end
        Xjrank = tiedrank(Xj)
        for i = 1:nc
            Yi = view(Y, :, i)
            if any(isnan, Yi)
                C[j,i] = NaN
            else
                Yirank = tiedrank(Yi)
                C[j,i] = cor(Xjrank, Yirank)
            end
        end
    end
    return C
end


#######################################
#
#   Kendall correlation
#
#######################################

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.
function corkendall!(x::RealVector, y::RealVector, permx::AbstractVector{<:Integer}=sortperm(x))
    if any(isnan, x) || any(isnan, y) return NaN end
    n = length(x)
    if n != length(y) error("Vectors must have same length") end

    # Initial sorting
    permute!(x, permx)
    permute!(y, permx)

    # Use widen to avoid overflows on both 32bit and 64bit
    npairs = div(widen(n) * (n - 1), 2)
    ntiesx = ndoubleties = nswaps = widen(0)
    k = 0

    @inbounds for i = 2:n
        if x[i - 1] == x[i]
            k += 1
        elseif k > 0
            # Sort the corresponding chunk of y, so the rows of hcat(x,y) are 
            # sorted first on x, then (where x values are tied) on y. Hence 
            # double ties can be counted by calling countties.
            sort!(view(y, (i - k - 1):(i - 1)))
            ntiesx += div(widen(k) * (k + 1), 2) # Must use wide integers here
            ndoubleties += countties(y,  i - k - 1, i - 1)
            k = 0
        end
    end
    if k > 0
        sort!(view(y, (n - k):n))
        ntiesx += div(widen(k) * (k + 1), 2)
        ndoubleties += countties(y, n - k, n)
    end

    nswaps = merge_sort!(y, 1, n)
    ntiesy = countties(y, 1, n)

    # Calls to float below prevent possible overflow errors when
    # length(x) exceeds 77_936 (32 bit) or 5_107_605_667 (64 bit)
    (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
     sqrt(float(npairs - ntiesx) * float(npairs - ntiesy))
end

"""
    corkendall(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors.
"""
corkendall(x::RealVector, y::RealVector) = corkendall!(copy(x), copy(y))

function corkendall(X::RealMatrix, y::RealVector)
    permy = sortperm(y)
    return([corkendall!(copy(y), X[:,i], permy) for i in 1:size(X, 2)])
end

function corkendall(x::RealVector, Y::RealMatrix)
    n = size(Y, 2)
    permx = sortperm(x)
    return(reshape([corkendall!(copy(x), Y[:,i], permx) for i in 1:n], 1, n))
end

function corkendall(X::RealMatrix)
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)
    for j = 2:n
        permx = sortperm(X[:,j])
        for i = 1:j - 1
            C[j,i] = corkendall!(X[:,j], X[:,i], permx)
            C[i,j] = C[j,i]
        end
    end
    return C
end

function corkendall(X::RealMatrix, Y::RealMatrix)
    nr = size(X, 2)
    nc = size(Y, 2)
    C = Matrix{Float64}(undef, nr, nc)
    for j = 1:nr
        permx = sortperm(X[:,j])
        for i = 1:nc
            C[j,i] = corkendall!(X[:,j], Y[:,i], permx)
        end
    end
    return C
end

# Auxilliary functions for Kendall's rank correlation

"""
    countties(x::RealVector, lo::Integer, hi::Integer)

Return the number of ties within `x[lo:hi]`. Assumes `x` is sorted. 
"""
function countties(x::AbstractVector, lo::Integer, hi::Integer)
    # Use of widen below prevents possible overflow errors when
    # length(x) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    thistiecount = result = widen(0)
    checkbounds(x, lo:hi)
    @inbounds for i = (lo + 1):hi
        if x[i] == x[i - 1]
            thistiecount += 1
        elseif thistiecount > 0
            result += div(thistiecount * (thistiecount + 1), 2)
            thistiecount = widen(0)
        end
    end

    if thistiecount > 0
        result += div(thistiecount * (thistiecount + 1), 2)
    end
    result
end

# Tests appear to show that a value of 64 is optimal,
# but note that the equivalent constant in base/sort.jl is 20.
const SMALL_THRESHOLD = 64

# merge_sort! copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
"""
    merge_sort!(v::AbstractVector, lo::Integer, hi::Integer, t::AbstractVector=similar(v, 0))    

Mutates `v` by sorting elements `x[lo:hi]` using the merge sort algorithm. 
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort distance.
"""
function merge_sort!(v::AbstractVector, lo::Integer, hi::Integer, t::AbstractVector=similar(v, 0))
    # Use of widen below prevents possible overflow errors when
    # length(v) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    nswaps = widen(0)
    @inbounds if lo < hi
        hi - lo <= SMALL_THRESHOLD && return insertion_sort!(v, lo, hi)

        m = midpoint(lo, hi)
        (length(t) < m - lo + 1) && resize!(t, m - lo + 1)

        nswaps = merge_sort!(v, lo,  m, t)
        nswaps += merge_sort!(v, m + 1, hi, t)

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

# insertion_sort! and midpoint copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
midpoint(lo::T, hi::T) where T <: Integer = lo + ((hi - lo) >>> 0x01)
midpoint(lo::Integer, hi::Integer) = midpoint(promote(lo, hi)...)

"""
    insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)

Mutates `v` by sorting elements `x[lo:hi]` using the insertion sort algorithm. 
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort distance.
"""
function insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)
    if lo == hi return widen(0) end
    nswaps = widen(0)
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
