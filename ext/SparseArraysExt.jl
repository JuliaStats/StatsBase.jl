module SparseArraysExt

using StatsBase
using SparseArrays

function StatsBase.indicatormat(::Type{SparseMatrixCSC}, x::AbstractArray{<:Integer}, k::Integer)
    n = length(x)
    return sparse(x, 1:n, true, k, n)
end

function StatsBase.indicatormat(::Type{SparseMatrixCSC}, x::AbstractArray{T}, c::AbstractArray{T}=sort!(unique(x))) where T
    d = indexmap(c)
    m = length(c)
    n = length(x)

    rinds = Vector{Int}(undef, n)
    for i = 1 : n
        rinds[i] = d[x[i]]
    end
    return sparse(rinds, 1:n, true, m, n)
end

end # module
