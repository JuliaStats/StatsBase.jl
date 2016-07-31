# Robust Statistics

#############################
#
#   Trimming outliers
#
#############################

# Trimmed set
"""
    trim(x, p=0.2)

Remove a proportion `p` of the highest elements, and `p` of the lowest
elements of vector `x` and return the result. To compute the trimmed
mean of 'x', use `mean(trim(x, p))`; to compute the variance of this
quantity, use `trimvar(x, p)`.


If you want to compute
the variance of `mean(trim(x,p))`, use `trimvar`.

# Example
```julia
julia> trim([1,2,3,4,5], 0.2)
 2
 3
 4
```
"""
trim(x::AbstractVector, p::Real=0.2) = trim!(copy(x), p)

"""
    trim!(x, p=0.2)

A variant of `trim` that partially sorts `x` in place while removing
the proportion `p` of the highest and lowest elements of `x`.
"""
function trim!(x::AbstractVector, p::Real=0.2)
    n = length(x)
    n > 0 || error("x can not be empty.")
    0 <= p < 0.5 || error("p must be in 0 <= p < 0.5.")
    g = floor(Int, n * p)
    
    select!(x, 1:g)
    select!(x, (n-g+1):n)
    
    return x[(g+1):(n-g)]
end

"""
    winsor(x, p=0.2)

Return a Winsorized version of vector `x`; i.e. replace the proportion
`p` of the lowest elements of `x` with the next-lowest, and replace the
proportion `p` of the highest elements with the previous-highest. This
function is used to compute the Winsorized mean of `x`. It is also used
by `trimvar` to compute the variance of the trimmed mean of `x`.

# Example
```julia
julia> winsor([1,2,3,4,5], 0.2)
 2
 2
 3
 4
 4
```
"""
winsor(x::AbstractVector, p::Real=0.2) = winsor!(copy(x), p)

"""
    winsor!(x, p=0.2)

A variant of `winsor` that modifies vector `x` in place.
"""
function winsor!(x::AbstractVector, p::Real=0.2)
    n = length(x)
    n > 0 || error("x can not be empty.")
    0 <= p < 0.5 || error("p must be in 0 <= p < 0.5.")
    g = floor(Int, n * p)
    
    select!(x, 1:g)
    select!(x, (n-g+1):n)
    x[1:g] = x[g+1]
    x[n-g+1:end] = x[n-g]
    
    return x
end


#############################
#
#   Other
#
#############################

# Variance of a trimmed set.
"""
    trimvar(x, p=0.2)

Compute the variance of the trimmed mean of `x`. This function uses
the Winsorized variance, as described in Wilcox (2010).
"""
function trimvar(x::RealArray, p::Real=0.2)
    n = length(x)
    n > 0 || error("x can not be empty.")
    0 <= p < 0.5 || error("p must be in 0 <= p < 0.5.")
    
    return var(winsor(x,p)) / (n * (1 - 2p)^2)
end