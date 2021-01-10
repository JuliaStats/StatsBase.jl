function _pairwise!(::Val{:none}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    m, n = size(dest)
    @inbounds for j in 1:n, i in 1:m
        symmetric && i > j && continue

        # For performance, diagonal is special-cased
        if f === cor && i == j && x[i] === y[j]
            # TODO: float() will not be needed after JuliaLang/Statistics.jl#61
            dest[i, j] = float(cor(x[i]))
        else
            dest[i, j] = f(x[i], y[j])
        end
    end
    if symmetric
        @inbounds for j in 1:n, i in (j+1):m
            dest[i, j] = dest[j, i]
        end
    end
    return dest
end

function check_vectors(x, y, skipmissing::Symbol)
    m = length(x)
    n = length(y)
    if !(all(xi -> xi isa AbstractVector, x) && all(yi -> yi isa AbstractVector, x))
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
    m, n = size(dest)
    @inbounds for j in 1:n
        ynminds = .!ismissing.(y[j])
        @inbounds for i in 1:m
            symmetric && i > j && continue

            if x[i] === y[j]
                ynm = view(y[j], ynminds)
                # For performance, diagonal is special-cased
                if f === cor && i == j
                    # If the type isn't concrete, 1 may not be converted to the right type
                    # and the final matrix will have an abstract eltype
                    # (missings and NaNs are ignored)
                    dest[i, j] = isconcretetype(eltype(dest)) ? 1 : one(f(ynm, ynm))
                else
                    dest[i, j] = f(ynm, ynm)
                end
            else
                nminds = .!ismissing.(x[i]) .& ynminds
                xnm = view(x[i], nminds)
                ynm = view(y[j], nminds)
                dest[i, j] = f(xnm, ynm)
            end
        end
    end
    if symmetric
        @inbounds for j in 1:n, i in (j+1):m
            dest[i, j] = dest[j, i]
        end
    end
    return dest
end

function _pairwise!(::Val{:listwise}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    check_vectors(x, y, :listwise)
    m, n = size(dest)
    nminds = .!ismissing.(x[1])
    @inbounds for i in 2:m
        nminds .&= .!ismissing.(x[i])
    end
    if x !== y
        @inbounds for j in 1:n
            nminds .&= .!ismissing.(y[j])
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

    # Pdesterve inferred element type
    isempty(dest) && return dest

    _pairwise!(f, dest, x′, y′, symmetric=symmetric, skipmissing=skipmissing)

    # identity.(dest) lets broadcasting compute a concrete element type
    # TODO: using promote_type rather than typejoin (which broadcast uses) would make sense
    # Once identity.(dest) is inferred automatically (JuliaLang/julia#30485),
    # the assertion can be removed
    @static if VERSION >= v"1.6.0-DEV"
        U = Base.Broadcast.promote_typejoin_union(Union{T, Tsm})
        return (isconcretetype(eltype(dest)) ? dest : identity.(dest))::Matrix{<:U}
    else
        return (isconcretetype(eltype(dest)) ? dest : identity.(dest))
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
"""
function pairwise!(f, dest::AbstractMatrix, x, y=x;
                   symmetric::Bool=false, skipmissing::Symbol=:none)
    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    return _pairwise!(f, dest, x, y, symmetric=symmetric, skipmissing=skipmissing)
end

# cov(x) is faster than cov(x, x)
pairwise!(::typeof(cov), dest::AbstractMatrix, x, y;
          symmetric::Bool=false, skipmissing::Symbol=:none) =
    pairwise!((x, y) -> x === y ? cov(x) : cov(x, y), dest, x, y,
              symmetric=symmetric, skipmissing=skipmissing)

pairwise!(::typeof(cor), dest::AbstractMatrix, x;
          symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise!(cor, dest, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise!(::typeof(cov), dest::AbstractMatrix, x;
         symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise!(cov, dest, x, x, symmetric=symmetric, skipmissing=skipmissing)

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
"""
function pairwise(f, x, y=x; symmetric::Bool=false, skipmissing::Symbol=:none)
    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    return _pairwise(Val(skipmissing), f, x, y, symmetric)
end

# cov(x) is faster than cov(x, x)
pairwise(::typeof(cov), x, y; symmetric::Bool=false, skipmissing::Symbol=:none) =
    pairwise((x, y) -> x === y ? cov(x) : cov(x, y), x, y,
             symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cor), x; symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise(cor, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cov), x; symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise(cov, x, x, symmetric=symmetric, skipmissing=skipmissing)