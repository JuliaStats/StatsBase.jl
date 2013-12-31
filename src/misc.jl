# Miscelleneous stuff

# run-length encoding
function rle{T}(v::Vector{T})
    n = length(v)
    vals = T[]
    lens = Int[]

    cv = v[1]
    cl = 1

    i = 2
    @inbounds while i <= n
        vi = v[i]
        if vi == cv
            cl += 1
        else
            push!(vals, cv)
            push!(lens, cl)
            cv = vi
            cl = 1
        end
        i += 1
    end

    # the last section
    push!(vals, cv)
    push!(lens, cl)

    return (vals, lens)
end

# inverse run-length encoding
function inverse_rle{T}(vals::AbstractVector{T}, lens::IntegerVector)
    m = length(vals)
    length(lens) == m || raise_dimerror()

    r = Array(T, sum(lens))
    p = 0
    @inbounds for i = 1 : m
        j = lens[i]
        v = vals[i]
        while j > 0
            r[p+=1] = v
            j -=1 
        end
    end
    return r
end


# findat (get positions (within a) for elements in b)

function indexmap{T}(a::AbstractArray{T})
    d = (T=>Int)[]
    for i = 1 : length(a)
        @inbounds k = a[i]
        if !haskey(d, k)
            d[k] = i
        end
    end
    return d
end

function findat!{T}(r::IntegerArray, a::AbstractArray{T}, b::AbstractArray{T})
    length(r) == length(b) || raise_dimerror()
    d = indexmap(a)
    @inbounds for i = 1 : length(b)
        r[i] = get(d, b[i], 0)
    end
    return r
end

findat(a::AbstractArray, b::AbstractArray) = findat!(Array(Int, size(b)), a, b)


# indicatormat

# x: input elements, 
# c: categories
# k: the maximum integer in x

function indicatormat(x::IntegerArray, k::Integer; sparse::Bool=false)
    sparse ? _indicatormat_sparse(x, k) : _indicatormat_dense(x, k)
end

function indicatormat(x::AbstractArray, c::AbstractArray; sparse::Bool=false)
    sparse ? _indicatormat_sparse(x, c) : _indicatormat_dense(x, c)
end

indicatormat{T<:Union(Real,String)}(x::AbstractArray{T}; sparse::Bool=false) = indicatormat(x, sort!(unique(x)); sparse=sparse)
indicatormat(x::AbstractArray; sparse::Bool=false) = indicatormat(x, unique(x); sparse=sparse)


function _indicatormat_dense(x::IntegerArray, k::Integer)
    n = length(x)
    r = zeros(Bool, k, n)
    for i = 1 : n
        r[x[i], i] = true
    end
    return r
end

function _indicatormat_dense{T}(x::AbstractArray{T}, c::AbstractArray{T})
    d = indexmap(c)
    m = length(c)
    n = length(x)
    r = zeros(Bool, m, n)
    o = 0
    @inbounds for i = 1 : n
        xi = x[i]
        r[o + d[xi]] = true
        o += m
    end
    return r
end

_indicatormat_sparse(x::IntegerArray, k::Integer) = (n = length(x); sparse(x, 1:n, true, k, n))

function _indicatormat_sparse{T}(x::AbstractArray{T}, c::AbstractArray{T})
    d = indexmap(c)
    m = length(c)
    n = length(x)

    rinds = Array(Int, n)
    @inbounds for i = 1 : n
        rinds[i] = d[x[i]]
    end
    return sparse(rinds, 1:n, true, m, n)
end

# Principal Component Analysis
# Performs PCA on the NxM data matrix X and returns the principal component coefficients and scores.
# Rows of X correspond to observations, columns to variables. This function will subtract the column 
# means of X, but will not normalize the variances. To get PCA with "standardized variables", use 
# pca(zscore(X))
function pca(x)
    x .-= mean(x, 1)
    # compute by singular value decomposition.
    u,s,v = svd(x, true)
    scores = scale(u, s) # equivalent to: x * v, but more efficient
    v, scores
end

# zscore(x) computes the centered, scaled version of x. When x is a vector, this is just 
# z = (x - mean(x)) ./ std(x). For array inputs, the mean and standard deviation are computed
# along the given dimension.
zscore(x, dim=1) =  (x .- mean(x, dim)) ./ std(x, dim)
