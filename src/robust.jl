# Robust Statistics

#############################
#
#   Trimming outliers
#
#############################

# Trimmed set
"""
    trim(x, p=0.2)

Remove a proportion `p` of its highest elements, and `p` of its lowest elements and return the result. If you want to compute the variance of `mean(trim(x,p))`, use `trimvar`.

# Example
```julia
julia> trim([1,2,3,4,5] , 0.2)
 2
 3
 4
```
"""
trim(x,p::Real=0.2) = trim!(copy(x),p)
trim!(x,p::Real=0.2) = trim!(collect(x), p)

"""
    trim!(x, p=0.2)

A variant of `trim` that modifies vector `x` in place. See also `trimvar`.
"""
function trim!(x::RealArray, p::Real=0.2)
    n = length(x)
    n > 0 || error("x can not be empty.")
    0 <= p < 0.5 || error("p must be in 0 <= p < 0.5.")
    g = floor(Int, n * p)
    
    return select!(x, (g+1):(n-g))
end

"""
    winsor(x, p=0.2)

Replace  the proportion `p` of the lowest elements with the next-lowest element, and the proportion `p` of the highest elements with the next highest elements. This function is typically used to compute the variance of `mean(trim(x,p))`.

# Example
```julia
julia> winsor([1,2,3,4,5] , 0.2)
 2
 2
 3
 4
 4
```
"""
winsor(x,p::Real=0.2) = winsor!(copy(x),p)
winsor!(x,p::Real=0.2) = winsor!(collect(x), p)

"""
    winsor!(x, p=0.2)

A variant of `winsor` that modifies vector `x` in place.
"""
function winsor!(x::RealArray, p::Real=0.2)
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
    trimvar(x,p=0.2)

Compute the variance of `mean(trim(x,p))` using the Winsorized variance, as described in Wilcox (2010).
"""
function trimvar(x::RealArray,p::Real=0.2)
    n = length(x)
    n > 0 || error("x can not be empty.")
    0 <= p < 0.5 || error("p must be in 0 <= p < 0.5.")
    
    return var(winsor(x,p)) / (n * (1 - 2p)^2)
end


