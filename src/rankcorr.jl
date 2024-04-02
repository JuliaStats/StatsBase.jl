# Rank-based correlations
#
# - Spearman's correlation
# - Kendall's correlation
#

#######################################
#
#   Spearman correlation
#
#######################################

"""
    corspearman(x, y=x; skipmissing::Symbol=:none)

Compute Spearman's rank correlation coefficient. If `x` and `y` are vectors, the
output is a float, otherwise it's a matrix corresponding to the pairwise correlations
of the columns of `x` and `y`.

Uses multiple threads when either `x` or `y` is a matrix and `skipmissing` is `:pairwise`.

# Keyword argument

- `skipmissing::Symbol=:none`: If `:none` (the default), `missing` entries in `x` or `y`
    give rise to `missing` entries in the return. If `:pairwise` when calculating an element
    of the return, both `i`th entries of the input vectors are skipped if either is missing.
    If `:listwise` the `i`th rows of both `x` and `y` are skipped if `missing` appears in
    either; note that this might skip a high proportion of entries. Only allowed when `x` or
    `y` is a matrix.
"""
function corspearman(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)
    check_rankcor_args(x, y, skipmissing, false)
    if x === y
        return corspearman(x)
    else
        return corspearman_kernel!(x, y, skipmissing)
    end
end

function corspearman(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return corspearman(x, reshape(y, (length(y), 1)); skipmissing=skipmissing)
end

function corspearman(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corspearman(reshape(x, (length(x), 1)), y; skipmissing=skipmissing)
end

function corspearman(x::AbstractMatrix, y::AbstractMatrix=x;
    skipmissing::Symbol=:none)
    check_rankcor_args(x, y, skipmissing, true)
    return pairwise(corspearman, _eachcol(x), _eachcol(y); skipmissing=skipmissing)
end

function corspearman(x::AbstractVector{T}) where {T}
    return T === Missing ? missing : 1.0
end

function _pairwise!(::Val{:listwise}, f::typeof(corspearman), dest::AbstractMatrix, x, y,
    symmetric::Bool)
    return _pairwise!(Val(:none), f, dest, handle_listwise(x, y)..., symmetric)
end

function _pairwise!(::Val{:none}, f::typeof(corspearman),
    dest::AbstractMatrix, x, y, symmetric::Bool)

    symmetric = x === y
    if symmetric && promoted_type(x) === Missing
        offdiag = (length(x[1]) < 2) ? NaN : missing
        @inbounds for i in axes(dest, 1), j in axes(dest, 2)
            dest[i, j] = i == j ? missing : offdiag
        end
        return dest
    end

    ranksx = ranks_matrix(x)
    if symmetric
        dest .= _cor(ranksx, ranksx)
    else
        ranksy = ranks_matrix(y)
        dest .= _cor(ranksx, ranksy)
    end

    #=When elements x[i] and y[i] are identical (according to `===`) then dest[i,i] should
    be 1.0 even in the presence of missing and NaN values. But the return from `_cor` does
    not respect that requirement. So amend.=#
    autocor = (eltype(dest) === Missing && skipmissing == :none) ? missing : 1.0
    @inbounds for i in 1:min(size(dest, 1), size(dest, 2))
        x[i] === y[i] && (dest[i, i] = autocor)
    end
    return dest
end

function _pairwise!(::Val{:pairwise}, f::typeof(corspearman),
    dest::AbstractMatrix{V}, x, y, symmetric::Bool) where {V}

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))
    symmetric = x === y

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest′ = reshape(dest, size(dest, 2), size(dest, 1))
        _pairwise!(Val(:pairwise), f, dest′, y, x, symmetric)
        dest .= transpose(dest′)
        return dest
    end

    tempx = sortperm_matrix(x)
    tempy = symmetric ? tempx : sortperm_matrix(y)

    int64 = Int64[]
    fl64 = Float64[]
    nmtx = promoted_nmtype(x)[]
    nmty = promoted_nmtype(y)[]
    Threads.@threads for subset in equal_sum_subsets(nr, Threads.nthreads())

        for j in subset

            inds = task_local_vector(:inds, int64, m)
            spnmx = task_local_vector(:spnmx, int64, m)
            spnmy = task_local_vector(:spnmy, int64, m)
            nmx = task_local_vector(:nmx, nmtx, m)
            nmy = task_local_vector(:nmy, nmty, m)
            ranksx = task_local_vector(:ranksx, fl64, m)
            ranksy = task_local_vector(:ranksy, fl64, m)

            for i = 1:(symmetric ? j : nc)
                # For performance, diagonal is special-cased
                if x[j] === y[i] && i == j && V !== Union{}
                    if missing isa V && eltype(x[j]) == Missing
                        dest[j, i] = missing
                    else
                        dest[j, i] = 1.0
                    end
                else
                    dest[j, i] = corspearman_kernel!(x[j], y[i], :pairwise,
                        view(tempx, :, j), view(tempy, :, i), inds, spnmx, spnmy, nmx,
                        nmy, ranksx, ranksy)
                end
                symmetric && (dest[i, j] = dest[j, i])
            end
        end
    end
    return dest
end

"""
    corspearman_kernel!(x::AbstractVector{T}, y::AbstractVector{U},
    skipmissing::Symbol, sortpermx=sortperm(x), sortpermy=sortperm(y),
    inds=zeros(Int64, length(x)), spnmx=zeros(Int64, length(x)),
    spnmy=zeros(Int64, length(x)), nmx=similar(x, nonmissingtype(T)),
    nmy=similar(y, nonmissingtype(U)), ranksx=similar(x, Float64),
    ranksy=similar(y, Float64)) where {T,U}

Compute Spearman's rank correlation coefficient between `x` and `y` with no allocations if
all arguments are provided.

All vector arguments must have the same axes.
- `sortpermx`: The sort permutation of `x`.
- `sortpermy`: The sort permutation of `y`.
- `inds` ... `ranksy`: Pre-allocated "scratch" arguments.

## Example
```julia-repl
julia> using StatsBase, BenchmarkTools, Random

julia> Random.seed!(1);

julia> x = ifelse.(rand(1000) .< 0.05,missing,randn(1000));y = ifelse.(rand(1000) .< 0.05, missing,randn(1000));

julia> sortpermx=sortperm(x);sortpermy=sortperm(y);inds=zeros(Int64,1000);spnmx=zeros(Int64,1000);spnmy=zeros(Int64,1000);

julia> nmx=zeros(1000);nmy=zeros(1000);ranksx=similar(x,Float64);ranksy=similar(y,Float64);

julia> @btime corspearman_kernel!(\$x,\$y,:pairwise,\$sortpermx,\$sortpermy,\$inds,\$spnmx,\$spnmy,\$nmx,\$nmy,\$ranksx,\$ranksy)
4.671 μs (0 allocations: 0 bytes)
-0.016058512110609713
```
"""
function corspearman_kernel!(x::AbstractVector{T}, y::AbstractVector{U},
    skipmissing::Symbol, sortpermx=sortperm(x), sortpermy=sortperm(y),
    inds=zeros(Int64, length(x)), spnmx=zeros(Int64, length(x)),
    spnmy=zeros(Int64, length(x)), nmx=similar(x, nonmissingtype(T)),
    nmy=similar(y, nonmissingtype(U)), ranksx=similar(x, Float64),
    ranksy=similar(y, Float64)) where {T,U}

    (axes(x) == axes(sortpermx) == axes(y) == axes(sortpermy) == axes(inds) ==
     axes(spnmx) == axes(spnmy) == axes(nmx) == axes(nmy) == axes(ranksx) ==
     axes(ranksy)) || throw(ArgumentError("Axes of inputs must match"))

    if skipmissing == :pairwise

        lb = first(axes(x, 1))
        k = lb
        #= Process (x,y) to obtain (nmx,nmy) by filtering out elements at position k if
        either x[k] or y[k] is missing. inds provides the mapping of elements of x (or y) to
        elements of nmx (or nmy) i.e. x[k] maps to nmx[inds[k]]. inds is then used to obtain
        spnmx and spnmy more efficiently than calling sortperm(nmx) and sortperm(nmy).
         =#
        @inbounds for i in axes(x, 1)
            if !(ismissing(x[i]) || ismissing(y[i]))
                inds[i] = k
                nmx[k] = x[i]
                nmy[k] = y[i]
                k += 1
            else
                inds[i] = lb - 1
            end
        end

        nnm = k - 1
        if nnm <= 1
            return NaN
        end
        nmx = view(nmx, lb:nnm)
        nmy = view(nmy, lb:nnm)

        if any(_isnan, nmx) || any(_isnan, nmy)
            return NaN
        end

        k = lb
        @inbounds for i in axes(x, 1)
            if (inds[sortpermx[i]]) != lb - 1
                spnmx[k] = inds[sortpermx[i]]
                k += 1
            end
        end
        spnmx = view(spnmx, lb:nnm)

        k = lb
        @inbounds for i in axes(y, 1)
            if (inds[sortpermy[i]]) != lb - 1
                spnmy[k] = inds[sortpermy[i]]
                k += 1
            end
        end
        spnmy = view(spnmy, lb:nnm)

        if nnm <= 1
            return NaN
        end

        _tiedrank!(view(ranksx, 1:nnm), nmx, spnmx)
        _tiedrank!(view(ranksy, 1:nnm), nmy, spnmy)

        return cor(view(ranksx, 1:nnm), view(ranksy, 1:nnm))

    else
        if length(x) <= 1
            return NaN
        elseif skipmissing == :none && (missing isa T || missing isa U) &&
               (any(ismissing, x) || any(ismissing, y))
            return missing
        elseif any(_isnan, x) || any(_isnan, y)
            return NaN
        end

        _tiedrank!(ranksx, x, sortpermx)
        _tiedrank!(ranksy, y, sortpermy)
        return cor(ranksx, ranksy)
    end
end

# Auxiliary functions for Spearman's rank correlation

"""
    _cor(x, y)
Work-around various "unhelpful" features of cor:
a) Ensures that on-diagonal elements of the return are as they should be in the symmetric
case i.e. 1.0 unless eltype(x) is Missing in which case on-diagonal elements should be
missing.
b) Ensure that _cor(a,b) is NaN when a and b are vectors of equal length less than 2
c) Works around some edge-case bugs in cor's handling of `missing` where the function throws
if `x` and `y` are matrices but nevertheless looping around the columns of `x` and `y`
works. https://github.com/JuliaStats/Statistics.jl/issues/63

# Example
```julia-repl
julia> x = y = fill(missing,2,2)
2×2 Matrix{Missing}:
 missing  missing
 missing  missing

julia> Statistics.cor(x,y)
ERROR: MethodError: no method matching copy(::Missing)

julia> StatsBase._cor(x,y)
2×2 Matrix{Union{Missing, Float64}}:
 1.0        missing
  missing  1.0

julia>

```
"""
function _cor(ranksx::AbstractMatrix{T}, ranksy::AbstractMatrix{U}) where {T,U}
    symmetric = ranksx === ranksy

    if size(ranksx, 1) < 2
        if symmetric && T === Missing
            return ifelse.(axes(ranksx, 2) .== axes(ranksx, 2)', missing, NaN)
        elseif symmetric
            return ifelse.(axes(ranksx, 2) .== axes(ranksx, 2)', 1.0, NaN)
        else
            return fill(NaN, size(ranksx, 2), size(ranksy, 2))
        end
    end
    try
        if symmetric
            return cor(ranksx)
        else
            return cor(ranksx, ranksy)
        end
    catch
        #=This catch block is hit when e.g.
        ranksx === ranksy = [missing missing;missing missing]
        =#
        nr, nc = size(ranksx, 2), size(ranksy, 2)
        if ranksx === ranksy && T === Missing
            return fill(missing, nr, nc)
        elseif missing isa T || missing isa U
            C = ones(Union{Missing,Float64}, nr, nc)
        else
            C = ones(Float64, nr, nc)
        end

        for j = (symmetric ? 2 : 1):nr
            for i = 1:(symmetric ? j - 1 : nc)
                C[j, i] = cor(view(ranksx, :, j), view(ranksy, :, i))
                symmetric && (C[i, j] = C[j, i])
            end
        end
        return C
    end
end

"""
    sortperm_matrix(x)

Given `x`, a vector of vectors, return a matrix who's ith column is the sort permutation of
the ith element of x.
"""
function sortperm_matrix(x)
    m = length(x) == 0 ? 0 : length(first(x))
    nc = length(x)
    int64 = Int64[]
    out = Array{Int}(undef, m, nc)

    Threads.@threads for i in 1:nc
        ints = task_local_vector(:ints, int64, m)
        sortperm!(ints, x[i])
        out[:, i] .= ints
    end
    return out
end

"""
    ranks_matrix(x)

Given `x`, a vector of vectors, return a matrix such that the (Pearson) correlaton between
columns of the return is the Spearman rank correlation between the elements of x.
"""
function ranks_matrix(x)

    m = length(x) == 0 ? 0 : length(first(x))
    nc = length(x)
    int64 = Int64[]

    if promoted_type(x) === Missing
        return fill(missing, m, nc)
    end

    out = Array{Union{Missing,Int,Float64}}(undef, m, nc)

    Threads.@threads for i in 1:nc
        ints = task_local_vector(:ints, int64, m)

        if any(ismissing, x[i])
            out[:, i] .= missing
        elseif any(_isnan, x[i])
            out[:, i] .= NaN
        else
            sortperm!(ints, x[i])
            _tiedrank!(view(out, :, i), x[i], ints)
        end
    end
    return out
end

#######################################
#
#   Kendall correlation
#
#######################################

"""
    corkendall(x, y=x; skipmissing::Symbol=:none)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must be either vectors or
matrices, and entries may be `missing`.

Uses multiple threads when either `x` or `y` is a matrix.

# Keyword argument

- `skipmissing::Symbol=:none`: If `:none` (the default), `missing` entries in `x` or `y`
    give rise to `missing` entries in the return. If `:pairwise` when calculating an element
    of the return, both `i`th entries of the input vectors are skipped if either is missing.
    If `:listwise` the `i`th rows of both `x` and `y` are skipped if `missing` appears in
    either; note that this might skip a high proportion of entries. Only allowed when `x` or
    `y` is a matrix.
"""
function corkendall(x::AbstractMatrix, y::AbstractMatrix=x;
    skipmissing::Symbol=:none)
    check_rankcor_args(x, y, skipmissing, true)
    return pairwise(corkendall, _eachcol(x), _eachcol(y); skipmissing=skipmissing)
end

function corkendall(x::AbstractVector, y::AbstractVector;
    skipmissing::Symbol=:none)
    check_rankcor_args(x, y, skipmissing, false)
    if x === y
        return corkendall(x)
    else
        x = copy(x)
        permx = sortperm(x)
        permute!(x, permx)
        return corkendall_kernel!(x, y, permx, skipmissing)
    end
end

_eachcol(x) = VERSION >= v"1.1" ? eachcol(x) : [view(x,:,j) for j in axes(x,2)]

#= corkendall returns a vector in this case, inconsistent with with Statistics.cor and
StatsBase.corspearman, but consistent with StatsBase.corkendall.
 =#
function corkendall(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return vec(corkendall(x, reshape(y, (length(y), 1)); skipmissing=skipmissing))
end

function corkendall(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corkendall(reshape(x, (length(x), 1)), y; skipmissing=skipmissing)
end

function corkendall(x::AbstractVector{T}) where {T}
    return T === Missing ? missing : 1.0
end

function _pairwise!(::Val{:none}, f::typeof(corkendall), dest::AbstractMatrix, x, y,
    symmetric::Bool)
    return corkendall_loop!(:none, f, dest, x, y, symmetric)
end

function _pairwise!(::Val{:pairwise}, f::typeof(corkendall), dest::AbstractMatrix, x, y,
    symmetric::Bool)
    return corkendall_loop!(:pairwise, f, dest, x, y, symmetric)
end

function _pairwise!(::Val{:listwise}, f::typeof(corkendall), dest::AbstractMatrix, x, y,
    symmetric::Bool)
    return corkendall_loop!(:none, f, dest, handle_listwise(x, y)..., symmetric)
end

function corkendall_loop!(skipmissing::Symbol, f::typeof(corkendall), dest::AbstractMatrix{V},
    x, y, symmetric::Bool) where {V}

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest′ = reshape(dest, size(dest, 2), size(dest, 1))
        corkendall_loop!(skipmissing, f, dest′, y, x, symmetric)
        dest .= transpose(dest′)
        return dest
    end

    intvec = Int[]
    t = promoted_type(x)[]
    u = promoted_type(y)[]
    t′ = promoted_nmtype(x)[]
    u′ = promoted_nmtype(y)[]

    symmetric = x === y

    Threads.@threads for subset in equal_sum_subsets(nr, Threads.nthreads())

        for j in subset

            sortedxj = task_local_vector(:sortedxj, t, m)
            scratch_py = task_local_vector(:scratch_py, u, m)
            yi = task_local_vector(:yi, u, m)
            permx = task_local_vector(:permx, intvec, m)
            # Ensuring missing is not an element type of scratch_sy, scratch_fx, scratch_fy
            # gives improved performance.
            scratch_sy = task_local_vector(:scratch_sy, u′, m)
            scratch_fx = task_local_vector(:scratch_fx, t′, m)
            scratch_fy = task_local_vector(:scratch_fy, t′, m)

            sortperm!(permx, x[j])
            @inbounds for k in eachindex(sortedxj)
                sortedxj[k] = x[j][permx[k]]
            end

            for i = 1:(symmetric ? j : nc)
                # For performance, diagonal is special-cased
                if x[j] === y[i] && i == j && V !== Union{}
                    if missing isa V && eltype(x[j]) == Missing
                        dest[j, i] = missing
                    else
                        dest[j, i] = 1.0
                    end
                else
                    yi .= y[i]
                    dest[j, i] = corkendall_kernel!(sortedxj, yi, permx, skipmissing;
                        scratch_py=scratch_py, scratch_sy=scratch_sy, scratch_fx=scratch_fx,
                        scratch_fy=scratch_fy)
                end
                symmetric && (dest[i, j] = dest[j, i])
            end
        end
    end
    return dest
end

# Auxiliary functions for Kendall's rank correlation

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.

function check_rankcor_args(x, y, skipmissing, allowlistwise::Bool)
    #  Base.require_one_based_indexing(x, y) #TODO find how to reject non-1-based in a way that works in Julia 1.0.5
    size(x, 1) == size(y, 1) ||
        throw(DimensionMismatch("x and y have inconsistent dimensions"))
    if allowlistwise
        skipmissing in (:none, :pairwise, :listwise) ||
            throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise, but got :$skipmissing"))
    else
        skipmissing in (:none, :pairwise) ||
            throw(ArgumentError("skipmissing must be either :none or :pairwise, but got :$skipmissing"))
    end
end

"""
    corkendall_kernel!(sortedx::AbstractVector, y::AbstractVector,
    permx::AbstractVector{<:Integer}, skipmissing::Symbol;
    scratch_py::AbstractVector=similar(y),
    scratch_sy::AbstractVector=similar(y),
    scratch_fx::AbstractVector=similar(sortedx),
    scratch_fy::AbstractVector=similar(y))

Kendall correlation between two vectors but omitting the initial sorting of the first
argument. Calculating Kendall correlation between `x` and `y` is thus a two stage process:
a) sort `x` to get `sortedx`; b) call this function on `sortedx` and `y`, with
subsequent arguments:
- `permx`: The permutation that sorts `x` to yield `sortedx`.
- `scratch_py`: A vector used to permute `y` without allocation.
- `scratch_sy`: A vector used to sort `y` without allocation.
- `scratch_fx, scratch_fy`: Vectors used to filter `missing`s from `x` and `y` without
   allocation.
"""
function corkendall_kernel!(sortedx::AbstractVector{T}, y::AbstractVector{U},
    permx::AbstractVector{<:Integer}, skipmissing::Symbol;
    scratch_py::AbstractVector{V}=similar(y),
    scratch_sy::AbstractVector=similar(y),
    scratch_fx::AbstractVector=similar(sortedx),
    scratch_fy::AbstractVector=similar(y)) where {T,U,V}

    length(sortedx) >= 2 || return NaN

    if skipmissing == :none
        if missing isa T && any(ismissing, sortedx)
            return missing
        elseif missing isa U && any(ismissing, y)
            return missing
        end
    end

    @inbounds for i in eachindex(y)
        scratch_py[i] = y[permx[i]]
    end

    if missing isa T || missing isa V
        sortedx, permutedy = handle_pairwise(sortedx, scratch_py; scratch_fx=scratch_fx, scratch_fy=scratch_fy)
    else
        permutedy = scratch_py
    end

    (any(_isnan, sortedx) || any(_isnan, permutedy)) && return NaN

    n = length(sortedx)

    # Use widen to avoid overflows on both 32bit and 64bit
    npairs = div(widen(n) * (n - 1), 2)
    ntiesx = ndoubleties = nswaps = widen(0)
    k = 0

    @inbounds for i = 2:n
        if sortedx[i-1] == sortedx[i]
            k += 1
        elseif k > 0
            #=
            Sort the corresponding chunk of permutedy, so rows of hcat(sortedx,permutedy)
            are sorted first on sortedx, then (where sortedx values are tied) on permutedy.
            Hence double ties can be counted by calling countties.
            =#
            sort!(view(permutedy, (i-k-1):(i-1)))
            ntiesx += div(widen(k) * (k + 1), 2) # Must use wide integers here
            ndoubleties += countties(permutedy, i - k - 1, i - 1)
            k = 0
        end
    end
    if k > 0
        sort!(view(permutedy, (n-k):n))
        ntiesx += div(widen(k) * (k + 1), 2)
        ndoubleties += countties(permutedy, n - k, n)
    end

    nswaps = merge_sort!(permutedy, 1, n, scratch_sy)
    ntiesy = countties(permutedy, 1, n)

    # Calls to float below prevent possible overflow errors when
    # length(sortedx) exceeds 77_936 (32 bit) or 5_107_605_667 (64 bit)
    return (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
           sqrt(float(npairs - ntiesx) * float(npairs - ntiesy))
end

"""
    countties(x::AbstractVector{<:Real}, lo::Integer, hi::Integer)

Return the number of ties within `x[lo:hi]`. Assumes `x` is sorted.
"""
function countties(x::AbstractVector, lo::Integer, hi::Integer)
    # Use of widen below prevents possible overflow errors when
    # length(x) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    thistiecount = result = widen(0)
    checkbounds(x, lo:hi)
    @inbounds for i = (lo+1):hi
        if x[i] == x[i-1]
            thistiecount += 1
        elseif thistiecount > 0
            result += div(thistiecount * (thistiecount + 1), 2)
            thistiecount = widen(0)
        end
    end

    if thistiecount > 0
        result += div(thistiecount * (thistiecount + 1), 2)
    end
    return result
end

# Tests appear to show that a value of 64 is optimal,
# but note that the equivalent constant in base/sort.jl is 20.
const SMALL_THRESHOLD = 64

# merge_sort! copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
"""
    merge_sort!(v::AbstractVector, lo::Integer, hi::Integer,
    t::AbstractVector=similar(v, 0))

Mutates `v` by sorting elements `x[lo:hi]` using the merge sort algorithm.
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort
distance.
"""
function merge_sort!(v::AbstractVector, lo::Integer, hi::Integer,
    t::AbstractVector=similar(v, 0))
    # Use of widen below prevents possible overflow errors when
    # length(v) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    nswaps = widen(0)
    @inbounds if lo < hi
        hi - lo <= SMALL_THRESHOLD && return insertion_sort!(v, lo, hi)

        m = midpoint(lo, hi)
        (length(t) < m - lo + 1) && resize!(t, m - lo + 1)

        nswaps = merge_sort!(v, lo, m, t)
        nswaps += merge_sort!(v, m + 1, hi, t)

        i, j = 1, lo
        while j <= m
            t[i] = v[j]
            i += 1
            j += 1
        end

        i, k = 1, lo
        while k < j <= hi
            if v[j] < t[i]
                v[k] = v[j]
                j += 1
                nswaps += m - lo + 1 - (i - 1)
            else
                v[k] = t[i]
                i += 1
            end
            k += 1
        end
        while k < j
            v[k] = t[i]
            k += 1
            i += 1
        end
    end
    return nswaps
end

# insertion_sort! and midpoint copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
midpoint(lo::T, hi::T) where {T<:Integer} = lo + ((hi - lo) >>> 0x01)
midpoint(lo::Integer, hi::Integer) = midpoint(promote(lo, hi)...)

"""
    insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)

Mutates `v` by sorting elements `x[lo:hi]` using the insertion sort algorithm.
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort
distance.
"""
function insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)
    if lo == hi
        return widen(0)
    end
    nswaps = widen(0)
    @inbounds for i = lo+1:hi
        j = i
        x = v[i]
        while j > lo
            if x < v[j-1]
                nswaps += 1
                v[j] = v[j-1]
                j -= 1
                continue
            end
            break
        end
        v[j] = x
    end
    return nswaps
end

# Auxiliary functions for both Kendall's and Spearman's rank correlations

# _isnan required so that corkendall and corspearman have correct handling of NaNs and
# can also accept arguments with element type for which isnan is not defined but isless is
# is defined, so that rank correlation makes sense.
_isnan(x::T) where {T<:Number} = isnan(x)
_isnan(x) = false

