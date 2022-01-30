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
    quantilerank(itr, value; method=:inc)

Compute the quantile(s)-position in [0, 1] interval of a `value` relative to a collection 
`itr`, e.g. a quantile rank of x means that (x*100)% of the elements in `itr` are less than 
(or less than or equal) the given `value`, according to the respective method used. 

Different definitions can be chosen via the `method` keyword argument.
Let `count_less` be the number of elements of `itr` that are less than `value`, 
`count_equal` the number of elements of `itr` that are equal to `value`, `n` the length of `itr`, 
`greatest_smaller` the highest value below `value` and `smallest_greater` the lowest value above `value`. 
Then `method` supports the following definitions:

`:inc` (default) - It calculates a value in the range 0 to 1 inclusive. 
If `value ∈ itr`, it returns `count_less / (n - 1)`, if not, apply interpolation based on 
def. 7 in Hyndman and Fan (1996). 
(equivalent to Excel `PERCENTRANK` and `PERCENTRANK.INC`)

`:exc` - It calculates a value in the range 0 to 1 exclusive. 
If `value ∈ itr`, it returns `(count_less + 1) / (n + 1)`,  if not, apply interpolation 
based on def. 6 in Hyndman and Fan (1996). 
(equivalent to Excel `PERCENTRANK.EXC`)

`:compete` - If `value ∈ itr`, it returns `count_less / (n - 1)`, if not, 
returns `(count_less - 1) / (n - 1)` and there is no interpolation.
(equivalent to MariaDB `PERCENT_RANK`, dplyr `percent_rank`)

`:tied` - Based on the def. in Roscoe, J. T. (1975), it returns 
`(count_less + count_equal/2) / n` and there is no interpolation. 
(equivalent to `"mean"` kind of SciPy `percentileofscore`)

`:strict` - It returns `count_less / n` and there is no interpolation.
(equivalent to `"strict"` kind of SciPy `percentileofscore`)

`:weak` - It returns `(count_less + count_equal) / n` and there is no interpolation.
(equivalent to `"weak"` kind of SciPy `percentileofscore`)

!!! note
    An `ArgumentError` is thrown if `itr` contains `NaN` or `missing` values
    or if `itr` contains fewer than two elements.

# References
Roscoe, J. T. (1975). "[Fundamental Research Statistics for the Behavioral Sciences 
(2nd ed.)](http://www.bryanburnham.net/wp-content/uploads/2014/07/Fundamental-Statistics-
for-the-Behavioral-Sciences-v2.0.pdf#page=57)".
    ISBN 0-03-091934-7.

Hyndman, R.J and Fan, Y. (1996) "[Sample Quantiles in Statistical Packages]
(https://www.amherst.edu/media/view/129116/original/Sample+Quantiles.pdf)",
    *The American Statistician*, Vol. 50, No. 4, pp. 361-365

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

# use `skipmissing` for vectors with missing entries.
julia> quantilerank(skipmissing(v2), 4)
0.5

# use broadcasting with `Ref` to compute quantile rank for multiple values
julia> quantilerank.(Ref(v3), [4, 8])
2-element Vector{Float64}:
 0.3333333333333333
 0.8888888888888888
```
"""
function quantilerank(itr, value; method::Symbol=:inc)
    ((value isa Number && isnan(value)) || ismissing(value)) &&
        throw(ArgumentError("`value` cannot be NaN or missing"))
    any(x -> ismissing(x) || (x isa Number && isnan(x)), itr) &&
        throw(ArgumentError("`itr` cannot contain missing or NaN entries"))

    count_less = count_equal = n = 0
    greatest_smaller = smallest_greater = value
    for x in itr
        if x == value
            count_equal += 1
        elseif x < value
            count_less += 1
            if greatest_smaller == value || greatest_smaller < x
                greatest_smaller = x
            end
        else
            if smallest_greater == value || smallest_greater > x
                smallest_greater = x
            end
        end
        n += 1
    end

    n == 0 && throw(ArgumentError("`itr` is empty. Pass a collection with at least two elements"))
    n == 1 && throw(ArgumentError("`itr` has only 1 value. Pass a collection with at least two elements"))

    if method == :inc
        if greatest_smaller == value
            return 0.0
        elseif count_equal > 0
            return count_less / (n - 1)
        elseif smallest_greater == value
            return 1.0
        else
            lower = (count_less - 1) / (n - 1)
            upper = count_less / (n - 1)
            ratio = (value - greatest_smaller) / (smallest_greater - greatest_smaller)
            return lower + ratio * (upper - lower)
        end
    elseif method == :exc
        if count_less == 0 && count_equal == 0
            return 0.0
        elseif count_less == 0
            return 1.0 / (n + 1)
        elseif count_equal > 0
            return (count_less + 1) / (n + 1)
        elseif smallest_greater == value
            return 1.0
        else
            lower = count_less / (n + 1)
            upper = (count_less + 1) / (n + 1)
            ratio = (value - greatest_smaller) / (smallest_greater - greatest_smaller)
            return lower + ratio * (upper - lower)
        end
    elseif method == :compete
        if value > maximum(itr)
            return 1.0
        elseif value ≤ minimum(itr) 
            return 0.0
        else
            value ∈ itr && (count_less += 1)
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
    percentilerank(itr, value; method=:inc)

Return the `q`th percentile of a collection `value`, i.e. [`quantilerank`](@ref) * 100.

Read the [`quantilerank`](@ref) docstring for more details of all available methods. 
"""
percentilerank(itr, value; method::Symbol=:inc) = quantilerank(itr, value, method=method) * 100
