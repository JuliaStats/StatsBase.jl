# Robust Statistics

#############################
#
#   Trimming outliers
#
#############################

# Trimmed set
"Return the upper and lower bound elements used by `trim` and `winsor`"
function uplo(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    n = length(x)
    n > 0 || throw(ArgumentError("x can not be empty."))

    if count == 0
        0 <= prop < 0.5 || throw(ArgumentError("prop must satisfy 0 ≤ prop < 0.5."))
        count = floor(Int, n * prop)
    else
        prop == 0 || throw(ArgumentError("prop and count can not both be > 0."))
        0 <= count < n/2 || throw(ArgumentError("count must satisfy 0 ≤ count < length(x)/2."))
    end

    # indices for lowest count values
    x2 = Base.copymutable(x)
    lo = partialsort!(x2, count+1)
    up = partialsort!(x2, n-count)

    up, lo
end

"""
    trim(x::AbstractVector; prop=0.0, count=0)

Return an iterator of all elements of `x` that omits either `count` or proportion
`prop` of the highest and lowest elements.

The number of trimmed elements could be smaller than specified if several
elements equal the lower or upper bound.

To compute the trimmed mean of `x` use `mean(trim(x))`;
to compute the variance use `trimvar(x)` (see [`trimvar`](@ref)).

# Example
```julia
julia> collect(trim([5,2,4,3,1], prop=0.2))
3-element Array{Int64,1}:
 2
 4
 3
```
"""
function trim(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    up, lo = uplo(x; prop=prop, count=count)

    (xi for xi in x if lo <= xi <= up)
end

"""
    trim!(x::AbstractVector; prop=0.0, count=0)

A variant of [`trim`](@ref) that modifies `x` in place.
"""
function trim!(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    up, lo = uplo(x; prop=prop, count=count)
    ix = (i for (i,xi) in enumerate(x) if lo > xi || xi > up)
    deleteat!(x, ix)
    return x
end

"""
    winsor(x::AbstractVector; prop=0.0, count=0)

Return an iterator of all elements of `x` that replaces either `count` or
proportion `prop` of the highest elements with the previous-highest element
and an equal number of the lowest elements with the next-lowest element.
                        
The number of replaced elements could be smaller than specified if several
elements equal the lower or upper bound.
                        
To compute the Winsorized mean of `x` use `mean(winsor(x))`.

# Example
```julia
julia> collect(winsor([5,2,3,4,1], prop=0.2))
5-element Array{Int64,1}:
 4
 2
 3
 4
 2
```
"""
function winsor(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    up, lo = uplo(x; prop=prop, count=count)

    (clamp(xi, lo, up) for xi in x)
end

"""
    winsor!(x::AbstractVector; prop=0.0, count=0)

A variant of [`winsor`](@ref) that modifies vector `x` in place.
"""
function winsor!(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    copyto!(x, winsor(x; prop=prop, count=count))
    return x
end


#############################
#
#   Other
#
#############################

# Variance of a trimmed set.
"""
    trimvar(x; prop=0.0, count=0)

Compute the variance of the trimmed mean of `x`. This function uses
the Winsorized variance, as described in Wilcox (2010).
"""
function trimvar(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    n = length(x)
    n > 0 || throw(ArgumentError("x can not be empty."))

    if count == 0
        0 <= prop < 0.5 || throw(ArgumentError("prop must satisfy 0 ≤ prop < 0.5."))
        count = floor(Int, n * prop)
    else
        0 <= count < n/2 || throw(ArgumentError("count must satisfy 0 ≤ count < length(x)/2."))
        prop = count/n
    end

    return var(winsor(x, count=count)) / (n * (1 - 2prop)^2)
end
