module StatsBaseSparseArraysExt

using SparseArrays
using StatsBase
import StatsBase: _indicatormat_sparse

_indicatormat_sparse(x::AbstractArray{<:Integer}, k::Integer) = (n = length(x); sparse(x, 1:n, true, k, n))

function _indicatormat_sparse(x::AbstractArray{T}, c::AbstractArray{T}) where T
    d = indexmap(c)
    m = length(c)
    n = length(x)

    rinds = Vector{Int}(undef, n)
    @inbounds for i = 1 : n
        rinds[i] = d[x[i]]
    end
    return sparse(rinds, 1:n, true, m, n)
end

end # module