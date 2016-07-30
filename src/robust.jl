# Robust Statistics

#############################
#
#   Trimming outliers
#
#############################

# Trimmed set
"""
    trim(x, [prop=0.0, count=0])


Remove either `count` or proportion `prop` of the highest and lowest
elements in `x` and return the result. To compute the trimmed mean of
`x`, use `mean(trim(x, prop=prop))`; to compute the variance of this
quantity, use `trimvar(x, prop=prop)`.

# Example
```julia
julia> trim([1,2,3,4,5], prop=0.2)
3-element Array{Int64,1}:
 2
 3
 4
```
"""
function trim(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    trim!(copy(x); prop=prop, count=count)
end

"""
    trim!(x, [prop=0.0, count=0])

A variant of `trim` that modifies `x` in place.
"""
function trim!(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    n = length(x)
    n > 0 || throw(ArgumentError("x can not be empty."))

    if count == 0
        0 <= prop < 0.5 || throw(ArgumentError("prop must satisfy 0 ≤ prop < 0.5."))
        count = floor(Int, n * prop)
    else
        0 <= count < n/2 || throw(ArgumentError("count must satisfy 0 ≤ count < length(x)/2."))
    end

    select!(x, (n-count+1):n)
    select!(x, 1:count)
    deleteat!(x, (n-count+1):n)
    deleteat!(x, 1:count)

    return x
end

"""
    winsor(x, [prop=0.0, count=0])

Return a Winsorized version of vector `x`; i.e. replace either `count` or
proportion `prop` of the lowest elements of `x` with the next-lowest, and
replace an equal number of the highest elements with the previous-highest.
This function is used to compute the Winsorized mean of `x`. It is also
used by `trimvar` to compute the variance of the trimmed mean of `x`.

# Example
```julia
julia> winsor([1,2,3,4,5], prop=0.2)
5-element Array{Int64,1}:
 2
 2
 3
 4
 4
```
"""
function winsor(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    winsor!(copy(x); prop=prop, count=count)
end

"""
    winsor!(x, [prop=0.0, count=0])

A variant of `winsor` that modifies vector `x` in place.
"""
function winsor!(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    n = length(x)
    n > 0 || throw(ArgumentError("x can not be empty."))

    if count == 0
        0 <= prop < 0.5 || throw(ArgumentError("prop must satisfy 0 ≤ prop < 0.5."))
        count = floor(Int, n * prop)
    else
        0 <= count < n/2 || throw(ArgumentError("count must satisfy 0 ≤ count < length(x)/2."))
    end

    select!(x, (n-count+1):n)
    select!(x, 1:count)
    x[1:count] = x[count+1]
    x[n-count+1:end] = x[n-count]

    return x
end


#############################
#
#   Other
#
#############################

# Variance of a trimmed set.
"""
    trimvar(x, [prop=0.0, count=0])

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

