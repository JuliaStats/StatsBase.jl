# a variety of rankings
#
# Please refer to http://en.wikipedia.org/wiki/Ranking#Strategies_for_assigning_rankings
# to see the definitions of a variety of ranking strategies
#
# The implementations here follow this wikipedia page.
#

function _check_randparams(rks, x, p)
    n = length(rks)
    nx = length(x)
    np = length(p)
    nx == np == n || throw(
        DimensionMismatch("lengths of x $nx and p $np do not match that of ranks $n"))
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
function _ordinalrank!(rks::AbstractArray, x::AbstractArray, p::AbstractArray{<:Integer})
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
function _competerank!(rks::AbstractArray, x::AbstractArray, p::AbstractArray{<:Integer})
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
function _denserank!(rks::AbstractArray, x::AbstractArray, p::AbstractArray{<:Integer})
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
function _tiedrank!(rks::AbstractArray, x::AbstractArray, p::AbstractArray{<:Integer})
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
