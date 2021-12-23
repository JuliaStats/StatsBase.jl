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
    quantilerank(v::AbstractVector, value; method=:inc, sorted=false)
    quantilerank(v::AbstractVector, values::AbstractVector; method=:inc, sorted=false)

Compute the quantile(s)-position [0-1] of a `value` or `values` relative to a collection `v`. 
The keyword argument `sorted` indicates whether `v` can be assumed to be sorted.

For example, a quantilerank of x% means that x% of the elements in `v` are below the given `value`.

The keyword argument `method` (default `:inc`) correspond to different calculation methodologies of the quantilerank, which are as follows:
- `:inc`     - Def. 7 in Hyndman and Fan (1996) of quantile. (Excel `PERCENTRANK` and `PERCENTRANK.INC`)
- `:exc`     - Def. 6 in Hyndman and Fan (1996) of quantile. (Excel `PERCENTRANK.EXC`)
- `:compete` - Based on the `competerank` of StatsBase (MariaDB `PERCENT_RANK`, dplyr `percent_rank`)
- `:mean`    - Def. in Roscoe, J. T. (1975) (equivalent to `mean` argument of Scipy `percentileofscore`)
- `:strict`  - (`strict` argument of Scipy `percentileofscore`)
- `:weak`    - (`weak` argument of Scipy `percentileofscore`)

!!! note
    An `ArgumentError` is thrown if `v` contains `NaN` or `missing` values.

# References
- [Percentile Rank on Wikipedia](https://en.wikipedia.org/wiki/Percentile_rank) covers definitions and examples.

- Roscoe, J. T. (1975). "[Fundamental Research Statistics for the Behavioral Sciences (2nd ed.)](http://www.bryanburnham.net/wp-content/uploads/2014/07/Fundamental-Statistics-for-the-Behavioral-Sciences-v2.0.pdf#page=57)".
    ISBN 0-03-091934-7.

- Hyndman, R.J and Fan, Y. (1996) "[Sample Quantiles in Statistical Packages](https://www.amherst.edu/media/view/129116/original/Sample+Quantiles.pdf)",
    *The American Statistician*, Vol. 50, No. 4, pp. 361-365

- [Quantile on Wikipedia](https://en.m.wikipedia.org/wiki/Quantile) details the different quantile definitions.


# Examples
```julia
julia> using StatsBase

julia> v1 = [1, 1, 1, 2, 3, 4, 8, 11, 12, 13];

julia> v2 = [1, 2, 3, 4, 4, 5, 6, 7, 8, 9];

julia> quantilerank(v1, 2), quantilerank(v1, 8)
(0.3333333333333333, 0.5555555555555556)

julia> quantilerank(v2, [4, 8])
2-element Vector{Float64}:
 0.3333333333333333
 0.8888888888888888
 
# also works with DataType
julia> using Dates

julia> d1 = Date("2021-01-01");

julia> daterange = d1:Day(1):(d1+Day(99)) |> collect;

julia> quantilerank(daterange, d1 + Day(20))
0.20202020202020202
```
"""
function quantilerank(v::AbstractVector, value; method::Symbol=:inc, sorted::Bool=false)
    # checks
    if value isa Number && eltype(v) <: Number
        isnan(value)    && throw(ArgumentError("value evaluated cannot contain NaN entries."))
        any(isnan.(v))  && throw(ArgumentError("vector cannot contain NaN entries."))
    end
    eltype(v) == Missing && throw(ArgumentError("vector cannot have only missings entries."))
    Missing <: eltype(v) && throw(ArgumentError("vector cannot contain missing entries. Use `skipmissing` function."))

    n = length(v)

    n == 0 && return 1.0

    if method == :inc
        if value ∈ v
            return count(<(value), v) / (n - 1)
        elseif value > maximum(v)
            return 1.0
        else
            !sorted && (v = sort(v))
            idx = searchsortedfirst(v, value)
            smrank, lgrank = 0, 0
            @inbounds for x in v
                x < value  && (smrank += 1)
                x < v[idx] && (lgrank += 1)
            end
            smallrank = (smrank - 1)  / (n - 1)
            largerank =  lgrank / (n - 1)
            step      = (value - v[idx-1])  / (v[idx] - v[idx-1])
            return smallrank + step*(largerank - smallrank)
        end
    elseif method == :exc
        if value ∈ v
            return (count(<(value), v) + 1) / (n + 1)
        elseif value > maximum(v)
            return 1.0
        else
            !sorted && (v = sort(v))
            idx  = searchsortedfirst(v, value)
            step = (value - v[idx-1]) / (v[idx] - v[idx-1])
            return ((idx - 1) + step) / (n + 1)
        end
    elseif method == :compete
        !sorted && (v = sort(v))
        return (searchsortedfirst(v, value) - 1) / (n - 1)
    elseif method == :mean
        smallrank, largerank = 0, 0
        for x in v
            x < value && (smallrank += 1)
            x ≤ value && (largerank += 1)
        end
        return (smallrank + largerank) / 2n
    elseif method == :strict
        return count(<(value), v) / n
    elseif method == :weak
        return count(≤(value), v) / n
    else
        throw(ArgumentError("method=:$method is not valid. Use :inc, :exc, :compete, :mean, :strict or :weak."))
    end
end

quantilerank(itr, value; kwargs...) = quantilerank(collect(itr), value; kwargs...)

function quantilerank(v::AbstractVector, values::AbstractVector; kwargs...)
    qtl = Vector{Float64}(undef, length(values))
    for i in eachindex(values)
        qtl[i] = quantilerank(v, values[i]; kwargs...)
    end
    return qtl
end

"""
    percentrank(v::AbstractVector, value)
    percentrank(v::AbstractVector, values::AbstractVector)

Return the `q`th percentile of a collection `value`, i.e. `quantilerank(v, value) * 100`.
"""
percentrank(v::AbstractVector, value; kwargs...) = quantilerank(v, value; kwargs...) * 100

percentrank(v::AbstractVector, values::AbstractVector; kwargs...) = quantilerank(v, values; kwargs...) .* 100

function percentrank(itr, value; kwargs...)
    v = collect(itr)
    eltype(v) == Missing && throw(ArgumentError("vector cannot have only missings entries."))
    Missing <: eltype(v) && throw(ArgumentError("vector cannot contain missing entries. Use `skipmissing` function."))
    return percentrank(v, value; kwargs...)
end