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
function ordinalrank!(rks::RealArray, x::RealArray, p::IntegerArray)
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
    ordinalrank(x)

Return the ordinal ranking ("1234" ranking) of a real-valued array.
All items in `x` are given distinct, successive ranks based on their
position in `sort(x)`.
"""
ordinalrank(x::RealArray) = ordinalrank!(Array{Int}(size(x)), x, sortperm(x))


# Competition ranking ("1224" ranking) -- resolve tied ranks using min
function competerank!(rks::RealArray, x::RealArray, p::IntegerArray)
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
    competerank(x)

Return the standard competition ranking ("1224" ranking) of a real-valued
array. Items that compare equal are given the same rank, then a gap is left
in the rankings the size of the number of tied items - 1.
"""
competerank(x::RealArray) = competerank!(Array{Int}(size(x)), x, sortperm(x))


# Dense ranking ("1223" ranking) -- resolve tied ranks using min
function denserank!(rks::RealArray, x::RealArray, p::IntegerArray)
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

Return the dense ranking ("1223" ranking) of a real-valued array. Items that
compare equal receive the same ranking, and the next subsequent rank is
assigned with no gap.
"""
denserank(x::RealArray) = denserank!(Array{Int}(size(x)), x, sortperm(x))


# Tied ranking ("1 2.5 2.5 4" ranking) -- resolve tied ranks using average
function tiedrank!(rks::RealArray, x::RealArray, p::IntegerArray)
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

Return the tied ranking, also called fractional or "1 2.5 2.5 4" ranking,
of a real-valued array. Items that compare equal receive the mean of the
rankings they would have been assigned under ordinal ranking.
"""
tiedrank(x::RealArray) = tiedrank!(Array{Float64}(size(x)), x, sortperm(x))

