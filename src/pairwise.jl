function _pairwise!(::Val{:none}, f, dest::AbstractMatrix{V}, x, y,
    symmetric::Bool) where {V}

    nr, nc = size(dest)

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest2 = reshape(dest, size(dest, 2), size(dest, 1))
        _pairwise!(Val(:none), f, dest2, y, x, symmetric)
        dest .= transpose(dest2)
        return dest
    end

    #cor and friends are special-cased.
    diag_is_1 = (f in (corkendall, corspearman, cor))
    (diag_is_1 || f == cov) && (symmetric = x === y)
    #cov(x) is faster than cov(x, x)
    (f == cov) && (f = ((x, y) -> x === y ? cov(x) : cov(x, y)))

    Threads.@threads for subset in equal_sum_subsets(nr, Threads.nthreads())
        for j in subset
            for i = 1:(symmetric ? j : nc)
                # For performance, diagonal is special-cased
                if diag_is_1 && i == j && y[i] === x[j] && V !== Union{}
                    dest[j, i] = V === Missing ? missing : 1.0
                else
                    dest[j, i] = f(x[j], y[i])
                end
            end
        end
    end
    symmetric && LinearAlgebra.copytri!(dest, 'L')
    return dest
end

#Input validation for both pairwise and pairwise!
function check_vectors(x, y, skipmissing::Symbol, symmetric::Bool)

    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    if !(skipmissing in (:none, :pairwise, :listwise))
        throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise"))
    end

    #When skipmissing is :none, elements of x/y can have unequal length.
    skipmissing == :none && return

    m = length(x)
    n = length(y)
    if !(all(xi -> xi isa AbstractVector, x) && all(yi -> yi isa AbstractVector, y))
        throw(ArgumentError("All entries in x and y must be vectors " *
                            "when skipmissing=:$skipmissing"))
    end
    if m > 1
        indsx = keys(first(x))
        for i in 2:m
            keys(x[i]) == indsx ||
                throw(ArgumentError("All input vectors must have the same indices"))
        end
    end
    if n > 1
        indsy = keys(first(y))
        for j in 2:n
            keys(y[j]) == indsy ||
                throw(ArgumentError("All input vectors must have the same indices"))
        end
    end
    if m > 1 && n > 1
        indsx == indsy ||
            throw(ArgumentError("All input vectors must have the same indices"))
    end
end

function _pairwise!(::Val{:pairwise}, f, dest::AbstractMatrix{V}, x, y, symmetric::Bool) where {V}

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest2 = reshape(dest, size(dest, 2), size(dest, 1))
        _pairwise!(Val(:pairwise), f, dest2, y, x, symmetric)
        dest .= transpose(dest2)
        return dest
    end

    #cor and friends are special-cased.
    diag_is_1 = (f in (corkendall, corspearman, cor))
    (diag_is_1 || f == cov) && (symmetric = x === y)
    #cov(x) is faster than cov(x, x)
    (f == cov) && (f = ((x, y) -> x === y ? cov(x) : cov(x, y)))

    nmtx = promoted_nmtype(x)[]
    nmty = promoted_nmtype(y)[]

    Threads.@threads for subset in equal_sum_subsets(nr, Threads.nthreads())
        scratch_fx = task_local_vector(:scratch_fx, nmtx, m)
        scratch_fy = task_local_vector(:scratch_fy, nmty, m)
        for j in subset
            for i = 1:(symmetric ? j : nc)
                if diag_is_1 && i == j && y[i] === x[j] && V !== Union{} && V !== Missing
                    dest[j, i] = 1.0
                else
                    dest[j, i] = f(handle_pairwise(x[j], y[i]; scratch_fx=scratch_fx, scratch_fy=scratch_fy)...)
                end
            end
        end
    end
    symmetric && LinearAlgebra.copytri!(dest, 'L')
    return dest
end

function _pairwise!(::Val{:listwise}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    return _pairwise!(Val(:none), f, dest, handle_listwise(x, y)..., symmetric)
end

function _pairwise!(f, dest::AbstractMatrix, x, y; symmetric::Bool=false,
    skipmissing::Symbol=:none)

    x′ = x isa Union{AbstractArray,Tuple,NamedTuple} ? x : collect(x)
    y′ = y isa Union{AbstractArray,Tuple,NamedTuple} ? y : collect(y)
    m = length(x′)
    n = length(y′)

    size(dest) != (m, n) &&
        throw(DimensionMismatch("dest has dimensions $(size(dest)) but expected ($m, $n)"))
    Base.has_offset_axes(dest) && throw("dest indices must start at 1")

    return _pairwise!(Val(skipmissing), f, dest, x′, y′, symmetric)
end

if VERSION >= v"1.6.0-DEV"
    # Function has moved in Julia 1.7
    if isdefined(Base, :typejoin_union_tuple)
        using Base: typejoin_union_tuple
    else
        using Base.Broadcast: typejoin_union_tuple
    end
else
    typejoin_union_tuple(::Type) = Any
end

# Identical to `Base.promote_typejoin` except that it uses `promote_type`
# instead of `typejoin` to combine members of `Union` types
function promote_type_union(::Type{T}) where {T}
    if T === Union{}
        return Union{}
    elseif T isa UnionAll
        return Any # TODO: compute more precise bounds
    elseif T isa Union
        return promote_type(promote_type_union(T.a), promote_type_union(T.b))
    elseif T <: Tuple
        return typejoin_union_tuple(T)
    else
        return T
    end
end

function _pairwise(::Val{skipmissing}, f, x, y, symmetric::Bool) where {skipmissing}
    x′ = x isa Union{AbstractArray,Tuple,NamedTuple} ? x : collect(x)
    y′ = y isa Union{AbstractArray,Tuple,NamedTuple} ? y : collect(y)
    m = length(x′)
    n = length(y′)

    T = Core.Compiler.return_type(f, Tuple{eltype(x′),eltype(y′)})
    Tsm = Core.Compiler.return_type((x, y) -> f(handle_pairwise(x, y)...),
        Tuple{eltype(x′),eltype(y′)})

    if skipmissing === :none
        dest = Matrix{T}(undef, m, n)
    elseif skipmissing in (:pairwise, :listwise)
        dest = Matrix{Tsm}(undef, m, n)
    else
        throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise"))
    end

    # Preserve inferred element type
    isempty(dest) && return dest

    _pairwise!(f, dest, x′, y′, symmetric=symmetric, skipmissing=skipmissing)

    if isconcretetype(eltype(dest))
        return dest
    else
        # Final eltype depends on actual contents (consistent with `map` and `broadcast`
        # but using `promote_type` rather than `promote_typejoin`)
        U = mapreduce(typeof, promote_type, dest)
        # V is inferred (contrary to U), but it only gives an upper bound for U
        V = promote_type_union(Union{T,Tsm})
        return convert(Matrix{U}, dest)::Matrix{<:V}
    end
end

"""
    pairwise!(f, dest::AbstractMatrix, x[, y];
              symmetric::Bool=false, skipmissing::Symbol=:none)

Store in matrix `dest` the result of applying `f` to all possible pairs
of entries in iterators `x` and `y`, and return it. Rows correspond to
entries in `x` and columns to entries in `y`, and `dest` must therefore
be of size `length(x) × length(y)`.
If `y` is omitted then `x` is crossed with itself.

As a special case, if `f` is `cor`, `corspearman` or `corkendall`, diagonal cells for
which entries from `x` and `y` are identical (according to `===`) are set to one even in the
presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted and ignored (taken as `true`)
  if `f` is `cov`, `cor`, `corkendall` or `corspearman`.
- `skipmissing::Symbol=:none`: If `:none` (the default), missing values
  in inputs are passed to `f` without any modification.
  Use `:pairwise` to skip entries with a `missing` value in either
  of the two vectors passed to `f` for a given pair of vectors in `x` and `y`.
  Use `:listwise` to skip entries with a `missing` value in any of the
  vectors in `x` or `y`; note that this might drop a large part of entries.
  Only allowed when entries in `x` and `y` are vectors.

# Examples
```jldoctest
julia> using StatsBase, Statistics

julia> dest = zeros(3, 3);

julia> x = [1 3 7
            2 5 6
            3 8 4
            4 6 2];

julia> pairwise!(cor, dest, eachcol(x));

julia> dest
3×3 Matrix{Float64}:
  1.0        0.744208  -0.989778
  0.744208   1.0       -0.68605
 -0.989778  -0.68605    1.0

julia> y = [1 3 missing
            2 5 6
            3 missing 2
            4 6 2];

julia> pairwise!(cor, dest, eachcol(y), skipmissing=:pairwise);

julia> dest
3×3 Matrix{Float64}:
  1.0        0.928571  -0.866025
  0.928571   1.0       -1.0
 -0.866025  -1.0        1.0
```
"""
function pairwise!(f, dest::AbstractMatrix, x, y=x;
    symmetric::Bool=false, skipmissing::Symbol=:none)
    check_vectors(x, y, skipmissing, symmetric)
    return _pairwise!(f, dest, x, y, symmetric=symmetric, skipmissing=skipmissing)
end

"""
    pairwise(f, x[, y];
             symmetric::Bool=false, skipmissing::Symbol=:none)

Return a matrix holding the result of applying `f` to all possible pairs
of entries in iterators `x` and `y`. Rows correspond to
entries in `x` and columns to entries in `y`. If `y` is omitted then a
square matrix crossing `x` with itself is returned.

As a special case, if `f` is `cor`, `corspearman` or `corkendall`, diagonal cells for
which entries from `x` and `y` are identical (according to `===`) are set to one even in the
presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted and ignored (taken as `true`)
  if `f` is `cov`, `cor`, `corkendall` or `corspearman`.
- `skipmissing::Symbol=:none`: If `:none` (the default), missing values
  in inputs are passed to `f` without any modification.
  Use `:pairwise` to skip entries with a `missing` value in either
  of the two vectors passed to `f` for a given pair of vectors in `x` and `y`.
  Use `:listwise` to skip entries with a `missing` value in any of the
  vectors in `x` or `y`; note that this might drop a large part of entries.
  Only allowed when entries in `x` and `y` are vectors.

# Examples
```jldoctest
julia> using StatsBase, Statistics

julia> x = [1 3 7
            2 5 6
            3 8 4
            4 6 2];

julia> pairwise(cor, eachcol(x))
3×3 Matrix{Float64}:
  1.0        0.744208  -0.989778
  0.744208   1.0       -0.68605
 -0.989778  -0.68605    1.0

julia> y = [1 3 missing
            2 5 6
            3 missing 2
            4 6 2];

julia> pairwise(cor, eachcol(y), skipmissing=:pairwise)
3×3 Matrix{Float64}:
  1.0        0.928571  -0.866025
  0.928571   1.0       -1.0
 -0.866025  -1.0        1.0
```
"""
function pairwise(f, x, y=x; symmetric::Bool=false, skipmissing::Symbol=:none)
    check_vectors(x, y, skipmissing, symmetric)
    return _pairwise(Val(skipmissing), f, x, y, symmetric)
end

# Auxiliary functions for pairwise

promoted_type(x) = mapreduce(eltype, promote_type, x, init=Union{})
promoted_nmtype(x) = mapreduce(nonmissingtype ∘ eltype, promote_type, x, init=Union{})

"""
    handle_listwise(x, y)

Remove missings in a listwise manner. Assumes `x` and `y` are vectors of iterables which
have been validated via `check_vectors`.

## Example
```julia-repl
julia> a = [6,7,8,9,10,missing];b = [4,5,6,7,missing,8];c = [2,3,4,missing,5,6];d = [1,2,3,4,5,6];

julia> StatsBase.handle_listwise([a,b],[c,d])
([[6, 7, 8], [4, 5, 6]], [[2, 3, 4], [1, 2, 3]])
```
"""
function handle_listwise(x, y)
    if !(missing isa promoted_type(x) || missing isa promoted_type(y))
        return x, y
    end

    nminds = .!ismissing.(first(x))
    @inbounds for xi in Iterators.drop(x, 1)
        nminds .&= .!ismissing.(xi)
    end
    if x !== y
        @inbounds for yj in y
            nminds .&= .!ismissing.(yj)
        end
    end

    # Computing integer indices once for all vectors is faster
    nminds′ = findall(nminds)
    # TODO: check whether wrapping views in a custom array type which asserts
    # that entries cannot be `missing` (similar to `skipmissing`)
    # could offer better performance

    x′ = [disallowmissing(view(xi, nminds′)) for xi in x]
    if x === y
        return x′, x′
    else
        y′ = [disallowmissing(view(yi, nminds′)) for yi in y]
        return x′, y′
    end
end

"""
    handle_pairwise(x::AbstractVector, y::AbstractVector;
    scratch_fx::AbstractVector=similar(x, nonmissingtype(eltype(x))),
    scratch_fy::AbstractVector=similar(y, nonmissingtype(eltype(y))))

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which elements `x[i]` and
`y[i]` are excluded if  `ismissing(x[i])||ismissing(y[i])`.
"""
function handle_pairwise(x::AbstractVector, y::AbstractVector;
    scratch_fx::AbstractVector=similar(x, nonmissingtype(eltype(x))),
    scratch_fy::AbstractVector=similar(y, nonmissingtype(eltype(y))))

    axes(x) == axes(y) || throw(DimensionMismatch("x and y have inconsistent dimensions"))
    lb = first(axes(x, 1))
    j = lb - 1
    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            scratch_fx[j] = x[i]
            scratch_fy[j] = y[i]
        end
    end

    return view(scratch_fx, lb:j), view(scratch_fy, lb:j)
end

#=Condition a) makes equal_sum_subsets useful for load-balancing the multi-threaded section
of _pairwise! in the non-symmetric case, and condition b) for the symmetric case.=#
"""
    equal_sum_subsets(n::Int, num_subsets::Int)::Vector{Vector{Int}}

Divide the integers 1:n into a number of subsets such that a) each subset has
(approximately) the same number of elements; and b) the sum of the elements in each subset
is nearly equal. If `n` is a multiple of `2 * num_subsets` both conditions hold exactly.

## Example
```julia-repl
julia> StatsBase.equal_sum_subsets(30,5)
5-element Vector{Vector{Int64}}:
 [30, 21, 20, 11, 10, 1]
 [29, 22, 19, 12, 9, 2]
 [28, 23, 18, 13, 8, 3]
 [27, 24, 17, 14, 7, 4]
 [26, 25, 16, 15, 6, 5]
```
"""
function equal_sum_subsets(n::Int, num_subsets::Int)::Vector{Vector{Int}}
    subsets = [Int[] for _ in 1:min(n, num_subsets)]
    writeto, scanup = 1, true
    for i = n:-1:1
        push!(subsets[writeto], i)
        if scanup && writeto == num_subsets
            scanup = false
        elseif (!scanup) && writeto == 1
            scanup = true
        else
            writeto += scanup ? 1 : -1
        end
    end
    return subsets
end

"""
    task_local_vector(key::Symbol, similarto::AbstractArray{V},
    length::Int)::Vector{V} where {V}

Retrieve from task local storage a vector of length `length` and matching the element
type of `similarto`, with initialisation on first call during a task.
"""
function task_local_vector(key::Symbol, similarto::AbstractArray{V},
    length::Int)::Vector{V} where {V}
    haskey(task_local_storage(), key) || task_local_storage(key, similar(similarto, length))
    return task_local_storage(key)
end