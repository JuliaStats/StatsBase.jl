# a variety of rankings
#
# Please refer to http://en.wikipedia.org/wiki/Ranking#Strategies_for_assigning_rankings
# to see the definitions of a variety of ranking strategies
#
# The implementations here follow this wikipedia page.
#

function _check_randparams(rks, x, p)
    n = length(rks)
    length(x) == length(p) == n || raise_dimerror()
    return n
end

# ranking helper function: calls sortperm(x) and then ranking method f!
function _rank(f!, x::AbstractArray, R::Type=Int; sortkwargs...)
    rks = similar(x, R)
    ord = reshape(sortperm(vec(x); sortkwargs...), size(x))
    return f!(rks, x, ord)
end

# ranking helper function for arrays with missing values
function _rank(f!, x::AbstractArray{>: Missing}, R::Type=Int; sortkwargs...)
    inds = findall(!ismissing, vec(x))
    isempty(inds) && return missings(R, size(x))
    xv = disallowmissing(view(vec(x), inds))
    ordv = sortperm(xv; sortkwargs...)
    rks = missings(R, size(x))
    f!(view(rks, inds), xv, ordv)
    return rks
end

# Ordinal ranking ("1234 ranking") -- use the literal order resulted from sort
function _ordinalrank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    _check_randparams(rks, x, p)
    @inbounds for i in eachindex(p)
        rks[p[i]] = i
    end
    return rks
end


"""
    ordinalrank(x; lt=isless, by=identity, rev::Bool=false, ...)

Return the [ordinal ranking](https://en.wikipedia.org/wiki/Ranking#Ordinal_ranking_.28.221234.22_ranking.29)
("1234" ranking) of an array. Supports the same keyword arguments as the `sort` function.
All items in `x` are given distinct, successive ranks based on their position
in the sorted vector.
Missing values are assigned rank `missing`.
"""
ordinalrank(x::AbstractArray; sortkwargs...) =
    _rank(_ordinalrank!, x; sortkwargs...)


# Competition ranking ("1224" ranking) -- resolve tied ranks using min
function _competerank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    @inbounds if n > 0
        p1 = p[1]
        v = x[p1]
        rks[p1] = k = 1

        for i in 2:n
            pi = p[i]
            xi = x[pi]
            if xi != v
                v = xi
                k = i
            end
            rks[pi] = k
        end
    end

    return rks
end


"""
    competerank(x; lt=isless, by=identity, rev::Bool=false, ...)

Return the [standard competition ranking](http://en.wikipedia.org/wiki/Ranking#Standard_competition_ranking_.28.221224.22_ranking.29)
("1224" ranking) of an array. Supports the same keyword arguments as the `sort` function.
Equal (*"tied"*) items are given the same rank, and the next rank comes after a gap
that is equal to the number of tied items - 1.
Missing values are assigned rank `missing`.
"""
competerank(x::AbstractArray; sortkwargs...) =
    _rank(_competerank!, x; sortkwargs...)


# Dense ranking ("1223" ranking) -- resolve tied ranks using min
function _denserank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    @inbounds if n > 0
        p1 = p[1]
        v = x[p1]
        rks[p1] = k = 1

        for i in 2:n
            pi = p[i]
            xi = x[pi]
            if xi != v
                v = xi
                k += 1
            end
            rks[pi] = k
        end
    end

    return rks
end


"""
    denserank(x; lt=isless, by=identity, rev::Bool=false, ...)

Return the [dense ranking](http://en.wikipedia.org/wiki/Ranking#Dense_ranking_.28.221223.22_ranking.29)
("1223" ranking) of an array. Supports the same keyword arguments as the `sort` function.
Equal items receive the same rank, and the next subsequent rank is
assigned with no gap.
Missing values are assigned rank `missing`.
"""
denserank(x::AbstractArray; sortkwargs...) =
    _rank(_denserank!, x; sortkwargs...)


# Tied ranking ("1 2.5 2.5 4" ranking) -- resolve tied ranks using average
function _tiedrank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    @inbounds if n > 0
        v = x[p[1]]

        s = 1  # starting index of current range
        for e in 2:n # e is pass-by-end index of current range
            cx = x[p[e]]
            if cx != v
                # fill average rank to s : e-1
                ar = (s + e - 1) / 2
                for i = s : e-1
                    rks[p[i]] = ar
                end
                # switch to next range
                s = e
                v = cx
            end
        end

        # the last range
        ar = (s + n) / 2
        for i = s : n
            rks[p[i]] = ar
        end
    end

    return rks
end

"""
    tiedrank(x; lt=isless, by=identity, rev::Bool=false, ...)

Return the [tied ranking](http://en.wikipedia.org/wiki/Ranking#Fractional_ranking_.28.221_2.5_2.5_4.22_ranking.29),
also called fractional or "1 2.5 2.5 4" ranking,
of an array. Supports the same keyword arguments as the `sort` function.
Equal (*"tied"*) items receive the mean of the ranks they would
have been assigned under the ordinal ranking (see [`ordinalrank`](@ref)).
Missing values are assigned rank `missing`.
"""
tiedrank(x::AbstractArray; sortkwargs...) =
    _rank(_tiedrank!, x, Float64; sortkwargs...)

"""
    quantilerank(v, value; method=:inc)

Compute the quantile(s)-position in [0-1] of a `value` relative to a collection `v`, e.g. a 
quantile rank of x means that (x*100)% of the elements in `v` are lesser (strict) or 
lesser-equal (weak) the given `value`. 

The function gets the counts of elements of `v` that are less than `value` (`count_less`), 
elements of `v` that are equal to `value` (`count_equal`) and the length of `v` (`n`). 
It also extracts the highest value of elements below `value` (`last_less`) and the lowest 
value of elements above `value` (`less_greater`). Using them, different methods and 
definitions of the `quantilerank` are obtained by changing the `method` keyword argument: 

`:inc` (default) - It calculates a value in the range 0 to 1 inclusive. 
If `value ∈ v`, it returns `count_less / (n - 1)`, if not, apply interpolation based on 
def. 7 in Hyndman and Fan (1996). 
(equivalent to Excel `PERCENTRANK` and `PERCENTRANK.INC`)

`:exc` - It calculates a value in the range 0 to 1 exclusive. 
If `value ∈ v`, it returns `(count_less + 1) / (n + 1)`,  if not, apply interpolation 
based on def. 6 in Hyndman and Fan (1996). 
(equivalent to Excel `PERCENTRANK.EXC`)

`:compete` - Based on the `competerank` function from StatsBase.jl, if `value ∈ v`, it returns 
`count_less / (n - 1)`, if not, returns `(count_less - 1) / (n - 1)`.
Also, there is no interpolation.
(equivalent to MariaDB `PERCENT_RANK`, dplyr `percent_rank`)

`:tied` - Based on the def. in Roscoe, J. T. (1975), it returns 
`(count_less + count_equal/2) / n` and there is no interpolation. 
(equivalent to `mean` argument of Scipy `percentileofscore`)

`:strict` - It returns `count_less / n` and there is no interpolation.
(equivalent to `strict` method of Scipy `percentileofscore`)

`:weak` - It returns `(count_less + count_equal) / n` and there is no interpolation.
(equivalent to `weak` method of Scipy `percentileofscore`)

!!! note
    An `ArgumentError` is thrown if `v` contains `NaN` or `missing` values
    or if `v` is empty or contains only one element.

# References
[Percentile Rank on Wikipedia](https://en.wikipedia.org/wiki/Percentile_rank) covers 
definitions and examples.

Roscoe, J. T. (1975). "[Fundamental Research Statistics for the Behavioral Sciences 
(2nd ed.)](http://www.bryanburnham.net/wp-content/uploads/2014/07/Fundamental-Statistics-
for-the-Behavioral-Sciences-v2.0.pdf#page=57)".
    ISBN 0-03-091934-7.

Hyndman, R.J and Fan, Y. (1996) "[Sample Quantiles in Statistical Packages]
(https://www.amherst.edu/media/view/129116/original/Sample+Quantiles.pdf)",
    *The American Statistician*, Vol. 50, No. 4, pp. 361-365

[Quantile on Wikipedia](https://en.m.wikipedia.org/wiki/Quantile) details the different 
quantile definitions.

# Examples
```julia
julia> using StatsBase

julia> v1 = [1, 1, 1, 2, 3, 4, 8, 11, 12, 13];

julia> v2 = [1, 2, 3, 5, 6, missing, 8];

julia> v3 = [1, 2, 3, 4, 4, 5, 6, 7, 8, 9];

julia> quantilerank(v1, 2)
0.3333333333333333

julia> quantilerank(v1, 2, method=:exc), quantilerank(v1, 2, method=:tied)
(0.36363636363636365, 0.35)

# use `skipmissing` in vectors with missing entries.
julia> quantilerank(skipmissing(v2), 4)
0.5

# use `Ref` to treat vector `v3` as a scalar during broadcasting.
```julia
julia> quantilerank.(Ref(v3), [4, 8])
2-element Vector{Float64}:
 0.3333333333333333
 0.8888888888888888
```
"""
function quantilerank(v, value; method::Symbol=:inc)
    # checks
    (value isa Number && isnan(value)) || ismissing(value) &&
        throw(ArgumentError("`value` cannot be NaN or missing"))
    any(x -> ismissing(x) || (x isa Number && isnan(x)), v) &&
        throw(ArgumentError("`v` cannot contain missing or NaN entries"))

    count_less = count_equal   = n = 0
    last_less  = first_greater = value
    for x in v
        if x == value
            count_equal += 1
        elseif x < value
            count_less += 1
            if last_less == value || last_less < x
                last_less = x
            end
        else
            if first_greater == value || first_greater > x
                first_greater = x
            end
        end
        n += 1
    end

    n == 0 && throw(ArgumentError("`v` is empty. Insert a non-empty vector"))
    n == 1 && throw(ArgumentError("`v` has only 1 value. Use a vector with more elements"))

    if method == :inc
        if last_less == value
            return 0.0
        elseif count_equal > 0
            return count_less / (n - 1)
        elseif first_greater == value
            return 1.0
        else
            lower = (count_less - 1) / (n - 1)
            upper = count_less / (n - 1)
            ratio = (value - last_less) / (first_greater - last_less)
            return lower + ratio * (upper - lower)
        end
    elseif method == :exc
        if count_less == 0 && count_equal == 0
            return 0.0
        elseif count_less == 0
            return 1.0 / (n + 1)
        elseif count_equal > 0
            return (count_less + 1) / (n + 1)
        elseif first_greater == value
            return 1.0
        else
            lower = count_less / (n + 1)
            upper = (count_less + 1) / (n + 1)
            ratio = (value - last_less) / (first_greater - last_less)
            return lower + ratio * (upper - lower)
        end
    elseif method == :compete
        if value > maximum(v)
            return 1.0
        elseif value ≤ minimum(v) 
            return 0.0
        else
            value ∈ v && (count_less += 1)
            return (count_less - 1) / (n - 1)
        end 
    elseif method == :tied
        return (count_less + count_equal/2) / n
    elseif method == :strict
        return count_less / n
    elseif method == :weak
        return (count_less + count_equal) / n
    else
        throw(ArgumentError("method=:$method is not valid. Use :inc, :exc, :compete, :tied, :strict or :weak."))
    end
end

"""
    percentrank(v, value; method=:inc)

Return the `q`th percentile of a collection `value`, i.e. [`quantilerank`](@ref) * 100.

Read the [`quantilerank`](@ref) docstring for more details of all available methods. 
"""
percentrank(v, value; method::Symbol=:inc) = quantilerank(v, value, method=method) * 100