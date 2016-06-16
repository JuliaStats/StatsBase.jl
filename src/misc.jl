# Miscelleneous stuff

# test whether a contains no repeated elements
function norepeat(a::AbstractArray)
    sa = sort(a)
    for i = 2:length(a)
        if a[i] == a[i-1]
            return false
        end
    end
    return true
end

# run-length encoding
"""
    rle(v) -> (vals, lens)

Return the run-length encoding of a vector as a tuple. The
first element of the tuple is a vector of values of the input
and the second is the number of consecutive occurrences of each
element.
"""
function rle{T}(v::Vector{T})
    n = length(v)
    vals = T[]
    lens = Int[]

    n>0 || return (vals,lens)

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
"""
    inverse_rle(vals, lens)

Reconstruct a vector from its run-length encoding. `vals` is a
vector of the values and `lens` is a vector of the corresponding
run lengths.
"""
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
"""
    indexmap(a)

Construct a dictionary that maps each unique value in `a` to
the index of its first occurrence in `a`.
"""
function indexmap{T}(a::AbstractArray{T})
    d = Dict{T,Int}()
    for i = 1 : length(a)
        @inbounds k = a[i]
        if !haskey(d, k)
            d[k] = i
        end
    end
    return d
end


"""
    levelsmap(a)

Construct a dictionary that maps each of the `n` unique values
in `a` to a number between 1 and `n`.
"""
function levelsmap{T}(a::AbstractArray{T})
    d = Dict{T,Int}()
    index = 1
    for i = 1 : length(a)
        @inbounds k = a[i]
        if !haskey(d, k)
            d[k] = index
            index += 1
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


"""
    findat(a, b)

For each element in `b`, find its first index in `a`. If the value does
not occur in `a`, the corresponding index is 0.
"""
findat(a::AbstractArray, b::AbstractArray) = findat!(Array(Int, size(b)), a, b)


# indicatormat

# x: input elements,
# c: categories
# k: the maximum integer in x

"""
    indicatormat(x, k::Integer; sparse=false)

Construct a boolean matrix `I` of size `(k, length(x))` such that
`I[x[i], i] = true` and all other elements are set to `false`.
If `sparse` is `true`, the output will be a sparse matrix, otherwise
it will be dense (default).
"""
function indicatormat(x::IntegerArray, k::Integer; sparse::Bool=false)
    sparse ? _indicatormat_sparse(x, k) : _indicatormat_dense(x, k)
end


"""
    indicatormat(x, [c=sort(unique(x))]; sparse=false)

Construct a boolean matrix `I` of size `(length(c), length(x))`.
Let `ci` be the index of `x[i]` in `c`. Then `I[ci, i] = true` and
all other elements are `false`.
"""
function indicatormat(x::AbstractArray, c::AbstractArray; sparse::Bool=false)
    sparse ? _indicatormat_sparse(x, c) : _indicatormat_dense(x, c)
end

indicatormat(x::AbstractArray; sparse::Bool=false) =
    indicatormat(x, sort!(unique(x)); sparse=sparse)


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

