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

ordinalrank(x::RealArray) = ordinalrank!(Array(Int, size(x)), x, sortperm(x))


# Competition ranking ("1224" ranking) -- resolve tied ranks using min
function competerank!(rks::RealArray, x::RealArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)
    
    if n > 0        
        v = x[1]
        rks[p[1]] = k = 1

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

competerank(x::RealArray) = competerank!(Array(Int, size(x)), x, sortperm(x))


# Dense ranking ("1223" ranking) -- resolve tied ranks using min
function denserank!(rks::RealArray, x::RealArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)
    
    if n > 0        
        v = x[1]
        rks[p[1]] = k = 1

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

denserank(x::RealArray) = denserank!(Array(Int, size(x)), x, sortperm(x))


# Tied ranking ("1 2.5 2.5 4" ranking) -- resolve tied ranks using average
function tiedrank!(rks::RealArray, x::RealArray, p::IntegerArray)
    n = _check_randparams(rks, x, p)

    i = 1

    while i <= n
        j = i
        while j + 1 <= n && x[p[i]] == x[p[j + 1]]
            j += 1
        end

        if j > i
            m = sum(i:j) / (j - i + 1)
            for k = i:j
                rks[p[k]] = m
            end
        else
            rks[p[i]] = i
        end

        i = j + 1
    end 

    return rks
end

# order (aka. rank), resolving ties using the mean rank
tiedrank(x::RealArray) = tiedrank!(Array(Float64, size(x)), x, sortperm(x))



