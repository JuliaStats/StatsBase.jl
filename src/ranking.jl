# a variety of rankings


# order (aka. rank), resolving ties using the mean rank
function tiedrank{T<:Real}(v::AbstractArray{T})
    n     = length(v)
    place = sortperm(v)
    ord   = Array(Float64, n)

    i = 1
    while i <= n
        j = i
        while j + 1 <= n && v[place[i]] == v[place[j + 1]]
            j += 1
        end

        if j > i
            m = sum(i:j) / (j - i + 1)
            for k = i:j
                ord[place[k]] = m
            end
        else
            ord[place[i]] = i
        end

        i = j + 1
    end

    return ord
end

tiedrank{T<:Real}(X::AbstractMatrix{T}) = tiedrank(reshape(X, length(X)))


