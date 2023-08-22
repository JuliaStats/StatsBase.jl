function _pairwise!(::Val{:none}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    @inbounds for (i, xi) in enumerate(x), (j, yj) in enumerate(y)
        symmetric && i > j && continue

        # For performance, diagonal is special-cased
        if f === cor && eltype(dest) !== Union{} && i == j && xi === yj
            # TODO: float() will not be needed after JuliaLang/Statistics.jl#61
            dest[i, j] = float(cor(xi))
        else
            dest[i, j] = f(xi, yj)
        end
    end
    if symmetric
        m, n = size(dest)
        @inbounds for j in 1:n, i in (j+1):m
            dest[i, j] = dest[j, i]
        end
    end
    return dest
end

function check_vectors(x, y, skipmissing::Symbol)
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

function _pairwise!(::Val{:pairwise}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    check_vectors(x, y, :pairwise)
    @inbounds for (j, yj) in enumerate(y)
        ynminds = .!ismissing.(yj)
        @inbounds for (i, xi) in enumerate(x)
            symmetric && i > j && continue

            if xi === yj
                ynm = view(yj, ynminds)
                # For performance, diagonal is special-cased
                if f === cor && eltype(dest) !== Union{} && i == j
                    # TODO: float() will not be needed after JuliaLang/Statistics.jl#61
                    dest[i, j] = float(cor(xi))
                else
                    dest[i, j] = f(ynm, ynm)
                end
            else
                nminds = .!ismissing.(xi) .& ynminds
                xnm = view(xi, nminds)
                ynm = view(yj, nminds)
                dest[i, j] = f(xnm, ynm)
            end
        end
    end
    if symmetric
        m, n = size(dest)
        @inbounds for j in 1:n, i in (j+1):m
            dest[i, j] = dest[j, i]
        end
    end
    return dest
end

function _pairwise!(::Val{:listwise}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    check_vectors(x, y, :listwise)
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
    return _pairwise!(Val(:none), f, dest,
                      [view(xi, nminds′) for xi in x],
                      [view(yi, nminds′) for yi in y],
                      symmetric)
end

function _pairwise!(f, dest::AbstractMatrix, x, y;
                    symmetric::Bool=false, skipmissing::Symbol=:none)
    if !(skipmissing in (:none, :pairwise, :listwise))
        throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise"))
    end

    x′ = x isa Union{AbstractArray, Tuple, NamedTuple} ? x : collect(x)
    y′ = y isa Union{AbstractArray, Tuple, NamedTuple} ? y : collect(y)
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
function promote_type_union(::Type{T}) where T
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
    x′ = x isa Union{AbstractArray, Tuple, NamedTuple} ? x : collect(x)
    y′ = y isa Union{AbstractArray, Tuple, NamedTuple} ? y : collect(y)
    m = length(x′)
    n = length(y′)

    T = Core.Compiler.return_type(f, Tuple{eltype(x′), eltype(y′)})
    Tsm = Core.Compiler.return_type((x, y) -> f(disallowmissing(x), disallowmissing(y)),
                                     Tuple{eltype(x′), eltype(y′)})

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
        V = promote_type_union(Union{T, Tsm})
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

As a special case, if `f` is `cor`, diagonal cells for which entries
from `x` and `y` are identical (according to `===`) are set to one even
in the presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted.
  Defaults to `true` when `f` is `cor` or `cov`.
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
    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    return _pairwise!(f, dest, x, y, symmetric=symmetric, skipmissing=skipmissing)
end

"""
    pairwise(f, x[, y];
             symmetric::Bool=false, skipmissing::Symbol=:none)

Return a matrix holding the result of applying `f` to all possible pairs
of entries in iterators `x` and `y`. Rows correspond to
entries in `x` and columns to entries in `y`. If `y` is omitted then a
square matrix crossing `x` with itself is returned.

As a special case, if `f` is `cor`, diagonal cells for which entries
from `x` and `y` are identical (according to `===`) are set to one even
in the presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted.
  Defaults to `true` when `f` is `cor` or `cov`.
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
    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    return _pairwise(Val(skipmissing), f, x, y, symmetric)
end

# cov(x) is faster than cov(x, x)
_cov(x, y) = x === y ? cov(x) : cov(x, y)

pairwise!(::typeof(cov), dest::AbstractMatrix, x, y;
          symmetric::Bool=false, skipmissing::Symbol=:none) =
    pairwise!(_cov, dest, x, y, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cov), x, y; symmetric::Bool=false, skipmissing::Symbol=:none) =
    pairwise(_cov, x, y, symmetric=symmetric, skipmissing=skipmissing)

pairwise!(::typeof(cov), dest::AbstractMatrix, x;
         symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise!(_cov, dest, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cov), x; symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise(_cov, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise!(::typeof(cor), dest::AbstractMatrix, x;
          symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise!(cor, dest, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cor), x; symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise(cor, x, x, symmetric=symmetric, skipmissing=skipmissing)
