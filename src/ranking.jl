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



# Ordinal ranking ("1234 ranking") -- use the literal order resulted from sort
function ordinalrank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    if n > 0
        i = 1
        while i <= n
            rks[p[i]] = i
            i += 1
        end
    end

    return rks
end


"""
    ordinalrank(x; lt = isless, rev::Bool = false)

Return the [ordinal ranking](https://en.wikipedia.org/wiki/Ranking#Ordinal_ranking_.28.221234.22_ranking.29)
("1234" ranking) of an array. The `lt` keyword allows providing a custom "less
than" function; use `rev=true` to reverse the sorting order.
All items in `x` are given distinct, successive ranks based on their
position in `sort(x; lt = lt, rev = rev)`.
Missing values are assigned rank `missing`.
"""
ordinalrank(x::AbstractArray; lt = isless, rev::Bool = false) =
    ordinalrank!(Array{Int}(undef, size(x)), x, sortperm(x; lt = lt, rev = rev))


# Competition ranking ("1224" ranking) -- resolve tied ranks using min
function competerank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    if n > 0
        p1 = p[1]
        v = x[p1]
        rks[p1] = k = 1

        i = 2
        while i <= n
            pi = p[i]
            xi = x[pi]
            if xi == v
                rks[pi] = k
            else
                rks[pi] = k = i
                v = xi
            end
            i += 1
        end
    end

    return rks
end


"""
    competerank(x; lt = isless, rev::Bool = false)

Return the [standard competition ranking](http://en.wikipedia.org/wiki/Ranking#Standard_competition_ranking_.28.221224.22_ranking.29)
("1224" ranking) of an array. The `lt` keyword allows providing a custom "less
than" function; use `rev=true` to reverse the sorting order.
Items that compare equal are given the same rank, then a gap is left
in the rankings the size of the number of tied items - 1.
Missing values are assigned rank `missing`.
"""
competerank(x::AbstractArray; lt = isless, rev::Bool = false) =
    competerank!(Array{Int}(undef, size(x)), x, sortperm(x; lt = lt, rev = rev))


# Dense ranking ("1223" ranking) -- resolve tied ranks using min
function denserank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    if n > 0
        p1 = p[1]
        v = x[p1]
        rks[p1] = k = 1

        i = 2
        while i <= n
            pi = p[i]
            xi = x[pi]
            if xi == v
                rks[pi] = k
            else
                rks[pi] = (k += 1)
                v = xi
            end
            i += 1
        end
    end

    return rks
end


"""
    denserank(x)

Return the [dense ranking](http://en.wikipedia.org/wiki/Ranking#Dense_ranking_.28.221223.22_ranking.29)
("1223" ranking) of an array. The `lt` keyword allows providing a custom "less
than" function; use `rev=true` to reverse the sorting order. Items that
compare equal receive the same ranking, and the next subsequent rank is
assigned with no gap.
Missing values are assigned rank `missing`.
"""
denserank(x::AbstractArray; lt = isless, rev::Bool = false) =
    denserank!(Array{Int}(undef, size(x)), x, sortperm(x; lt = lt, rev = rev))


# Tied ranking ("1 2.5 2.5 4" ranking) -- resolve tied ranks using average
function tiedrank!(rks::AbstractArray, x::AbstractArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    if n > 0
        v = x[p[1]]

        s = 1  # starting index of current range
        e = 2  # pass-by-end index of current range
        while e <= n
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
            e += 1
        end

        # the last range (e == n+1)
        ar = (s + n) / 2
        for i = s : n
            rks[p[i]] = ar
        end
    end

    return rks
end

# order (aka. rank), resolving ties using the mean rank
"""
    tiedrank(x)

Return the [tied ranking](http://en.wikipedia.org/wiki/Ranking#Fractional_ranking_.28.221_2.5_2.5_4.22_ranking.29),
also called fractional or "1 2.5 2.5 4" ranking,
of an array. The `lt` keyword allows providing a custom "less
than" function; use `rev=true` to reverse the sorting order.
Items that compare equal receive the mean of the
rankings they would have been assigned under ordinal ranking.
Missing values are assigned rank `missing`.
"""
tiedrank(x::AbstractArray; lt = isless, rev::Bool = false) =
    tiedrank!(Array{Float64}(undef, size(x)), x, sortperm(x; lt = lt, rev = rev))

for (f, f!, S) in zip([:ordinalrank, :competerank, :denserank, :tiedrank],
                      [:ordinalrank!, :competerank!, :denserank!, :tiedrank!],
                      [Int, Int, Int, Float64])
    @eval begin
        function $f(x::AbstractArray{>: Missing}; lt = isless, rev::Bool = false)
            inds = findall(!ismissing, x)
            isempty(inds) && return missings($S, size(x))
            xv = disallowmissing(view(x, inds))
            sp = sortperm(xv; lt = lt, rev = rev)
            rks = missings($S, length(x))
            $(f!)(view(rks, inds), xv, sp)
            rks
        end
    end
end