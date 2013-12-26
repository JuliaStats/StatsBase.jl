# Functions for tabulation

function table{T}(a::AbstractArray{T})
    counts = Dict{T, Int}()
    for i = 1:length(a)
        tmp = a[i]
        counts[tmp] = get(counts, tmp, 0) + 1
    end
    return counts
end

