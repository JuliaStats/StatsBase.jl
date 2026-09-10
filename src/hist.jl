using Base.Cartesian

import Base: show, ==, push!, append!, float
import LinearAlgebra: norm, normalize, normalize!


## Fast getindex function for multiple arrays, returns a tuple of array elements
@inline Base.@propagate_inbounds @generated function _multi_getindex(i::Union{Integer, CartesianIndex}, c::AbstractArray...)
    N = length(c)
    result_expr = Expr(:tuple)
    for j in 1:N
        push!(result_expr.args, :(c[$j][i]))
    end
    result_expr
end


# Need a generated function to promote edge types, because a simple
# promote_type(map(eltype, h.edges)...) isn't type stable (tested
# with Julia v0.5).
@generated function _promote_edge_types(edges::NTuple{N,AbstractVector}) where N
    promote_type(map(eltype, edges.parameters)...)
end


"""
    UniformEdges{F} <: AbstractVector{F}

Bin edges of equal width, as produced by `fit(Histogram, v; nbins)`. The edges are stored
as a `Vector{F}` and `step(edges)` gives the bin width. Because each edge is the value of a
decimal number rounded to `F`, differences of neighbouring edges may deviate from `step` by
an ulp, but the bins are equal-width by construction.
"""
struct UniformEdges{F<:AbstractFloat} <: AbstractVector{F}
    edges::Vector{F}
    step::F
end

Base.size(e::UniformEdges) = size(e.edges)
Base.IndexStyle(::Type{<:UniformEdges}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(e::UniformEdges, i::Int) = e.edges[i]
Base.step(e::UniformEdges) = e.step
# Print like a range: first edge, width, last edge
Base.show(io::IO, e::UniformEdges) = print(io, first(e), ':', step(e), ':', last(e))
Base.show(io::IO, ::MIME"text/plain", e::UniformEdges) = show(io, e)

## nice-valued edges for histograms
function histrange(v::AbstractArray{T}, n::Integer, closed::Symbol=:left) where T
    F = float(T)
    nv = length(v)
    if nv == 0 && n < 0
        throw(ArgumentError("number of bins must be ≥ 0 for an empty array, got $n"))
    elseif nv > 0 && n < 1
        throw(ArgumentError("number of bins must be ≥ 1 for a non-empty array, got $n"))
    elseif nv == 0
        return UniformEdges([zero(F)], one(F))
    end

    lo, hi = extrema(v)
    histrange(F(lo), F(hi), n, closed)
end

# Return `UniformEdges{F}` of strictly increasing bin edges covering `[lo, hi]` with
# approximately `n` bins of equal width, where the width is a "nice" decimal number: 1, 2 or 5
# times a power of ten.
# The edges are the decimal numbers `k * 10^e` for consecutive multiples `k` of the width,
# each rounded to the nearest `F`. Rounding each edge individually is what makes this
# work for every floating point type: no arithmetic progression has to be represented in F,
# and an observation that is the rounding of the same decimal as an edge compares equal to
# that edge. The endpoints are then adjusted by comparing `lo` and `hi` against the F
# edges themselves, i.e. against exactly what `binindex` compares against, so no
# observation can fall outside the edges (#1009).
function histrange(lo::F, hi::F, n::Integer, closed::Symbol=:left) where F<:AbstractFloat
    isfinite(lo) && isfinite(hi) ||
        throw(ArgumentError("histogram edges cannot be computed for non-finite data"))
    # Choose the decimal exponent e and the multiplier m ∈ {1, 2, 5, 10} of the bin width
    # from the raw width in F. Only the integers e and m are used from here on, so the
    # rounding of these few operations does not affect the exactness of the edges.
    # If all values are identical, use a single bin of unit width around the value.
    bw = hi == lo ? one(F) : hi / n - lo / n
    # The nice width below is at least bw/1.1. For consecutive decimal edges to round to
    # distinct values of F, the width must exceed twice the floating-point spacing of the
    # edges, and the spacing at the last edge can be twice that at hi when the edge lies in
    # the next binade. Requiring bw ≥ 3 eps thus guarantees strictly increasing edges.
    bw = max(bw, 3 * eps(lo), 3 * eps(hi))
    lbw = log10(bw)
    e0 = floor(Int, lbw)
    r = exp10(lbw - e0)  # bw / 10^e0 in [1, 10), without forming 10^e0 in F
    m0 = r <= 1.1 ? 1 : r <= 2.2 ? 2 : r <= 5.5 ? 5 : 10
    # write 10 * 10^e0 as 1 * 10^(e0 + 1) so that the multiples k below stay as small as possible
    m, e = m0 == 10 ? (1, e0 + 1) : (m0, e0)
    # Multiples of the width near lo and hi. The width is at least three times the spacing of
    # the data, so these quotients are at most 2^precision / 3 and their rounding error is
    # below one; the adjustment below corrects the remaining off-by-one.
    edge(k) = _decimal(F, k, e)
    stepF = edge(m)
    K = _multiple_type(F)
    kfirst = m * floor(K, lo / stepF)
    klast = m * ceil(K, hi / stepF)
    # Adjust the endpoints against the F edges so that lo and hi are inside the edges, with
    # the first edge as large and the last edge as small as possible. Note that e.g. the
    # Float32 observation 0.7f0 is below the decimal 0.7 but equals the Float32 edge 0.7f0
    if closed == :right #(,]
        while edge(kfirst + m) < lo
            kfirst += m
        end
        while lo <= edge(kfirst)
            kfirst -= m
        end
        while edge(klast - m) >= hi
            klast -= m
        end
        while edge(klast) < hi
            klast += m
        end
    else #[,)
        while edge(kfirst + m) <= lo
            kfirst += m
        end
        while lo < edge(kfirst)
            kfirst -= m
        end
        while edge(klast - m) > hi
            klast -= m
        end
        while edge(klast) <= hi
            klast += m
        end
    end
    return UniformEdges(F[edge(k) for k in kfirst:m:klast], stepF)
end

# Integer type for the multiples of the width, which are bounded by 5/3 * maxintfloat(F)
_multiple_type(::Type{F}) where F<:AbstractFloat = 2 * maxintfloat(F) <= typemax(Int) ? Int : BigInt

# The decimal number k * 10^e rounded to F. When k and 10^|e| are both exactly representable
# in F, the single multiplication or division is correctly rounded by IEEE arithmetic; this
# covers all but extreme magnitudes. Beyond that, the power of ten is taken from a table of
# correctly rounded powers and at most three roundings occur (k, the power, the product or
# quotient), so the result is within two ulps of the decimal. Exact conversion for all
# exponents is the job of a decimal parser and is not attempted here. Note that `F(k // 10^d)`
# is not an option: Base performs that division in F itself, so the denominator is already
# rounded for F = Float32 from 10^11 on (JuliaLang/julia#49749).
function _decimal(::Type{F}, k::Integer, e::Integer) where F<:AbstractFloat
    if abs(k) <= maxintfloat(F) && abs(e) <= _maxexp10(F)
        p = F(10)^abs(e)
        return e >= 0 ? F(k) * p : F(k) / p
    elseif e >= 0
        return F(k) * _pow10(F, e)
    else
        p = _pow10(F, -e)
        isfinite(p) && return F(k) / p
        # 10^(-e) overflows F although k * 10^e may not: divide by an exact power first
        a = _maxexp10(F)
        return F(k) / F(10)^a / _pow10(F, -e - a)
    end
end

# 10^d for d ≥ 0 rounded to the nearest F, from tables for the IEEE types (Inf beyond floatmax)
const _POW10_FLOAT64 = [Float64(big(10)^d) for d in 0:308]
const _POW10_FLOAT32 = [Float32(big(10)^d) for d in 0:38]
const _POW10_FLOAT16 = [Float16(big(10)^d) for d in 0:4]
_pow10(::Type{Float64}, d::Integer) = d < length(_POW10_FLOAT64) ? @inbounds(_POW10_FLOAT64[d + 1]) : Inf
_pow10(::Type{Float32}, d::Integer) = d < length(_POW10_FLOAT32) ? @inbounds(_POW10_FLOAT32[d + 1]) : Inf32
_pow10(::Type{Float16}, d::Integer) = d < length(_POW10_FLOAT16) ? @inbounds(_POW10_FLOAT16[d + 1]) : Inf16
_pow10(::Type{F}, d::Integer) where F<:AbstractFloat = F(big(10)^d)

# The largest d such that 10^d = 2^d * 5^d is exactly representable in F, i.e. 5^d < 2^precision(F)
_maxexp10(::Type{Float16}) = 4
_maxexp10(::Type{Float32}) = 10
_maxexp10(::Type{Float64}) = 22
function _maxexp10(::Type{F}) where F<:AbstractFloat
    # 5^d computed in F is exact while below 2^precision and at least 2^precision otherwise,
    # so the comparison decides exactly; start from the closed form and correct if needed
    lim = ldexp(one(F), precision(F))
    d = floor(Int, precision(F) / log2(5))
    while F(5)^d >= lim
        d -= 1
    end
    while F(5)^(d + 1) < lim && isfinite(F(10)^(d + 1))
        d += 1
    end
    return d
end

histrange(vs::NTuple{N,AbstractVector},nbins::NTuple{N,Integer},closed::Symbol) where {N} =
    map((v,n) -> histrange(v,n,closed),vs,nbins)
histrange(vs::NTuple{N,AbstractVector},nbins::Integer,closed::Symbol) where {N} =
    map(v -> histrange(v,nbins,closed),vs)



## histograms ##
function sturges(n)  # Sturges' formula
    n==0 && return one(n)
    ceil(Integer, log2(n))+1
end

abstract type AbstractHistogram{T<:Real,N,E} end

# N-dimensional histogram object
"""
    Histogram <: AbstractHistogram

The `Histogram` type represents data that has been tabulated into intervals
(known as *bins*) along the real line, or in higher dimensions, over a real space.
Histograms can be fitted to data using the `fit` method.

# Fields
* edges: An iterator that contains the boundaries of the bins in each dimension.
* weights: An array that contains the weight of each bin.
* closed: A symbol with value `:right` or `:left` indicating on which side bins
  (half-open intervals or higher-dimensional analogues thereof) are closed.
  See below for an example.
* isdensity: There are two interpretations of a `Histogram`. If `isdensity=false` the weight of a bin corresponds to the amount of a quantity in the bin.
  If `isdensity=true` then it corresponds to the density (amount / volume) of the quantity in the bin. See below for an example.

# Examples
## Example illustrating `closed`
```jldoctest
julia> using StatsBase

julia> fit(Histogram, [2.],  1:3, closed=:left)
Histogram{Int64, 1, Tuple{UnitRange{Int64}}}
edges:
  1:3
weights: [0, 1]
closed: left
isdensity: false

julia> fit(Histogram, [2.],  1:3, closed=:right)
Histogram{Int64, 1, Tuple{UnitRange{Int64}}}
edges:
  1:3
weights: [1, 0]
closed: right
isdensity: false
```
## Example illustrating `isdensity`
```jldoctest
julia> using StatsBase, LinearAlgebra

julia> bins = [0,1,7]; # a small and a large bin

julia> obs = [0.5, 1.5, 1.5, 2.5]; # one observation in the small bin and three in the large

julia> h = fit(Histogram, obs, bins)
Histogram{Int64, 1, Tuple{Vector{Int64}}}
edges:
  [0, 1, 7]
weights: [1, 3]
closed: left
isdensity: false

julia> # observe isdensity = false and the weights field records the number of observations in each bin

julia> normalize(h, mode=:density)
Histogram{Float64, 1, Tuple{Vector{Int64}}}
edges:
  [0, 1, 7]
weights: [1.0, 0.5]
closed: left
isdensity: true

julia> # observe isdensity = true and weights tells us the number of observation per binsize in each bin
```
"""
mutable struct Histogram{T<:Real,N,E} <: AbstractHistogram{T,N,E}
    edges::E
    weights::Array{T,N}
    closed::Symbol
    isdensity::Bool
    function Histogram{T,N,E}(edges::NTuple{N,AbstractArray}, weights::Array{T,N},
                              closed::Symbol, isdensity::Bool=false) where {T,N,E}
        closed == :right || closed == :left || error("closed must :left or :right")
        isdensity && !(T <: AbstractFloat) && error("Density histogram must have float-type weights")
        _edges_nbins(edges) == size(weights) || error("Histogram edge vectors must be 1 longer than corresponding weight dimensions")
        new{T,N,E}(edges,weights,closed,isdensity)
    end
end

Histogram(edges::NTuple{N,AbstractVector}, weights::AbstractArray{T,N},
          closed::Symbol=:left, isdensity::Bool=false) where {T,N} =
    Histogram{T,N,typeof(edges)}(edges,weights,closed,isdensity)

Histogram(edges::NTuple{N,AbstractVector}, ::Type{T}, closed::Symbol=:left,
          isdensity::Bool=false) where {T,N} =
    Histogram(edges,zeros(T,_edges_nbins(edges)...),closed,isdensity)

Histogram(edges::NTuple{N,AbstractVector}, closed::Symbol=:left,
          isdensity::Bool=false) where {N} =
    Histogram(edges,Int,closed,isdensity)

function show(io::IO, h::AbstractHistogram)
    println(io, typeof(h))
    println(io,"edges:")
    for e in h.edges
        println(IOContext(io, :limit => true),"  ",e)
    end
    println(io,"weights: ",h.weights)
    println(io,"closed: ",h.closed)
    print(io,"isdensity: ",h.isdensity)
end

(==)(h1::Histogram,h2::Histogram) = (==)(h1.edges,h2.edges) && (==)(h1.weights,h2.weights) && (==)(h1.closed,h2.closed) && (==)(h1.isdensity,h2.isdensity)


binindex(h::AbstractHistogram{T,1}, x::Real) where {T} = binindex(h, (x,))[1]

binindex(h::Histogram{T,N}, xs::NTuple{N,Real}) where {T,N} =
    map((edge, x) -> _edge_binindex(edge, h.closed, x), h.edges, xs)

# Compare with `<` rather than the default `isless`: `isless(-0.0, 0.0)` is true, so -0.0
# would be binned differently from 0.0, whereas `-0.0 < 0.0` is false and the two are
# treated as equal. NaN compares false with everything under `<`, so it ends up outside
# the edges (bin index 0 or `length(edge)`) and is dropped by `push!` as before. `<` is also
# cheaper than `isless` and keeps the arithmetic fast path for ranges.
@inline function _edge_binindex(edge::AbstractVector, closed::Symbol, x::Real)
    if closed === :right
        return searchsortedfirst(edge, x, lt = <) - 1
    else
        return searchsortedlast(edge, x, lt = <)
    end
end
# For equal-width bins, estimate the index from the width and correct it against the stored
# edges. The estimate is off by at most one, so this is a few operations instead of a binary
# search, and the result is exactly what the search would give.
@inline function _edge_binindex(edge::UniformEdges, closed::Symbol, x::Real)
    v = edge.edges
    n = length(v)
    @inbounds begin
        lo = v[1]
        hi = v[n]
        if closed === :right
            # number of edges strictly below x
            lo < x || return 0
            hi < x && return n
            i = clamp(floor(Int, (x - lo) / edge.step) + 1, 1, n - 1)
            while i > 0 && !(v[i] < x)
                i -= 1
            end
            while i < n && v[i + 1] < x
                i += 1
            end
        else
            # number of edges at or below x
            lo <= x || return 0
            hi <= x && return n
            i = clamp(floor(Int, (x - lo) / edge.step) + 1, 1, n - 1)
            while i > 0 && !(v[i] <= x)
                i -= 1
            end
            while i < n && v[i + 1] <= x
                i += 1
            end
        end
        return i
    end
end


binvolume(h::AbstractHistogram{T,1}, binidx::Integer) where {T} = binvolume(h, (binidx,))
binvolume(::Type{V}, h::AbstractHistogram{T,1}, binidx::Integer) where {V,T} = binvolume(V, h, (binidx,))

binvolume(h::Histogram{T,N}, binidx::NTuple{N,Integer}) where {T,N} =
    binvolume(_promote_edge_types(h.edges), h, binidx)

binvolume(::Type{V}, h::Histogram{T,N}, binidx::NTuple{N,Integer}) where {V,T,N} =
    prod(map((edge, i) -> _edge_binvolume(V, edge, i), h.edges, binidx))

@inline _edge_binvolume(::Type{V}, edge::AbstractVector, i::Integer) where {V} = V(edge[i+1]) - V(edge[i])
@inline _edge_binvolume(::Type{V}, edge::AbstractRange, i::Integer) where {V} = V(step(edge))
@inline _edge_binvolume(::Type{V}, edge::UniformEdges, i::Integer) where {V} = V(step(edge))
@inline _edge_binvolume(edge::AbstractVector, i::Integer) = _edge_binvolume(eltype(edge), edge, i)


@inline _edges_nbins(edges::NTuple{N,AbstractVector}) where {N} = map(_edge_nbins, edges)

@inline _edge_nbins(edge::AbstractVector) = length(edge) - 1


# 1-dimensional

Histogram(edge::AbstractVector, weights::AbstractVector{T}, closed::Symbol=:left, isdensity::Bool=false) where {T} =
    Histogram((edge,), weights, closed, isdensity)

Histogram(edge::AbstractVector, ::Type{T}, closed::Symbol=:left, isdensity::Bool=false) where {T} =
    Histogram((edge,), T, closed, isdensity)

Histogram(edge::AbstractVector, closed::Symbol=:left, isdensity::Bool=false) =
    Histogram((edge,), closed, isdensity)


push!(h::AbstractHistogram{T,1}, x::Real, w::Real) where {T} = push!(h, (x,), w)
push!(h::AbstractHistogram{T,1}, x::Real) where {T} = push!(h,x,one(T))
append!(h::AbstractHistogram{T,1}, v::AbstractVector) where {T} = append!(h, (v,))
append!(h::AbstractHistogram{T,1}, v::AbstractVector, wv::Union{AbstractVector,AbstractWeights}) where {T} = append!(h, (v,), wv)

fit(::Type{Histogram{T}},v::AbstractVector, edg::AbstractVector; closed::Symbol=:left) where {T} =
    fit(Histogram{T},(v,), (edg,), closed=closed)
fit(::Type{Histogram{T}},v::AbstractVector; closed::Symbol=:left, nbins=sturges(length(v))) where {T} =
    fit(Histogram{T},(v,); closed=closed, nbins=nbins)
fit(::Type{Histogram{T}},v::AbstractVector, wv::AbstractWeights, edg::AbstractVector; closed::Symbol=:left) where {T} =
    fit(Histogram{T},(v,), wv, (edg,), closed=closed)
fit(::Type{Histogram{T}},v::AbstractVector, wv::AbstractWeights; closed::Symbol=:left, nbins=sturges(length(v))) where {T} =
    fit(Histogram{T}, (v,), wv; closed=closed, nbins=nbins)

fit(::Type{Histogram}, v::AbstractVector, wv::AbstractWeights{W}, args...; kwargs...) where {W} = fit(Histogram{W}, v, wv, args...; kwargs...)

# N-dimensional

function push!(h::Histogram{T,N},xs::NTuple{N,Real},w::Real) where {T,N}
    h.isdensity && error("Density histogram must have float-type weights")
    idx = binindex(h, xs)
    if checkbounds(Bool, h.weights, idx...)
        h.weights[idx...] += w
    end
    h
end

function push!(h::Histogram{T,N},xs::NTuple{N,Real},w::Real) where {T<:AbstractFloat,N}
    idx = binindex(h, xs)
    if checkbounds(Bool, h.weights, idx...)
        h.weights[idx...] += h.isdensity ? w / binvolume(h, idx) : w
    end
    h
end

push!(h::AbstractHistogram{T,N},xs::NTuple{N,Real}) where {T,N} = push!(h,xs,one(T))


function append!(h::AbstractHistogram{T,N}, vs::NTuple{N,AbstractVector}) where {T,N}
    for i in eachindex(vs...)
        xs = _multi_getindex(i, vs...)
        push!(h, xs, one(T))
    end
    h
end
function append!(h::AbstractHistogram{T,N}, vs::NTuple{N,AbstractVector}, wv::AbstractVector) where {T,N}
    for i in eachindex(wv, vs...)
        xs = _multi_getindex(i, vs...)
        push!(h, xs, wv[i])
    end
    h
end

# Turn kwargs nbins into a type-stable tuple of integers:
function _nbins_tuple(vs::NTuple{N,AbstractVector}, nbins) where N
    template = map(length, vs)
    result = broadcast((t, x) -> typeof(t)(x), template, nbins)
    result::typeof(template)
end

fit(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}, edges::NTuple{N,AbstractVector}; closed::Symbol=:left) where {T,N} =
    append!(Histogram(edges, T, closed, false), vs)

fit(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}; closed::Symbol=:left, nbins=sturges(length(vs[1]))) where {T,N} =
    fit(Histogram{T}, vs, histrange(vs,_nbins_tuple(vs, nbins),closed); closed=closed)

fit(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}, wv::AbstractWeights{W}, edges::NTuple{N,AbstractVector}; closed::Symbol=:left) where {T,N,W} =
    append!(Histogram(edges, T, closed, false), vs, wv)

fit(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}, wv::AbstractWeights; closed::Symbol=:left, nbins=sturges(length(vs[1]))) where {T,N} =
    fit(Histogram{T}, vs, wv, histrange(vs,_nbins_tuple(vs, nbins),closed); closed=closed)

"""
    fit(Histogram, data[, weight][, edges]; closed=:left[, nbins])

Fit a histogram to `data`.

# Arguments

* `data`: either a vector (for a 1-dimensional histogram), or a tuple of
  vectors of equal length (for an *n*-dimensional histogram).

* `weight`: an optional `AbstractWeights` (of the same length as the
  data vectors), denoting the weight each observation contributes to the
  bin. If no weight vector is supplied, each observation has weight 1.

* `edges`: a vector (for example an `AbstractRange` object), or tuple of vectors, that gives
  the edges of the bins along each dimension. If no edges are provided, they are chosen
  so that approximately `nbins` bins of equal width are constructed along each dimension.

!!! note
    In most cases, the number of bins will be `nbins`. However, to ensure that the bins have
    equal width, more or fewer than `nbins` bins may be used. The automatically chosen bin
    width is a "nice" decimal number (1, 2 or 5 times a power of ten) and the edges are
    multiples of it rounded to the floating point type of the data, returned as
    [`StatsBase.UniformEdges`](@ref), a vector of edges which also records the bin width as `step`.
    All observations are guaranteed to fall inside the automatically chosen edges.
    For data of extreme magnitude (beyond about `1e±22` for `Float64`, `1e±10` for `Float32`),
    the edges may differ from the decimal number by up to two units in the last place.

!!! note
    Observations that fall outside the supplied `edges` are not counted. No error or warning
    is raised, so the edges should span the data unless dropping observations is intended.
    For the observations and the edges to compare as the decimal numbers they were written
    as, the edges should have the same floating point type as the data. For example,
    `0.7f0` is slightly smaller than the `Float64` value `0.7` and would be counted in the bin
    below `0.7` if the edges are `Float64`, but in the bin starting at `0.7f0` if the edges
    are `Float32`.

# Keyword arguments

* `closed`: if `:left` (the default), the bin intervals are left-closed [a,b);
  if `:right`, intervals are right-closed (a,b].

* `nbins`: if no `edges` argument is supplied, the approximate number of bins to use
  along each dimension (can be either a single integer, or a tuple of integers).
  If omitted, it is computed using Sturges's formula, i.e. `ceil(log2(length(n))) + 1`
  with `n` the number of data points.

# Examples

```julia
# Univariate
h = fit(Histogram, rand(100))
h = fit(Histogram, rand(100), 0:0.1:1.0)
h = fit(Histogram, rand(100), nbins=10)
h = fit(Histogram, rand(100), weights(rand(100)), 0:0.1:1.0)
h = fit(Histogram, [20], 0:20:100)
h = fit(Histogram, [20], 0:20:100, closed=:right)

# Multivariate
h = fit(Histogram, (rand(100),rand(100)))
h = fit(Histogram, (rand(100),rand(100)),nbins=10)
```
"""
fit(::Type{Histogram}, args...; kwargs...) = fit(Histogram{Int}, args...; kwargs...)
fit(::Type{Histogram}, vs::NTuple{N,AbstractVector}, wv::AbstractWeights{W}, args...; kwargs...) where {N,W} = fit(Histogram{W}, vs, wv, args...; kwargs...)


# Get a suitable high-precision type for the norm of a histogram.
norm_type(h::Histogram{T,N}) where {T,N} =
    promote_type(T, _promote_edge_types(h.edges))

norm_type(::Type{T}) where {T<:Integer} = promote_type(T, Int64)
norm_type(::Type{T}) where {T<:AbstractFloat} = promote_type(T, Float64)


"""
    norm(h::Histogram)

Calculate the norm of histogram `h` as the absolute value of its integral.
"""
@generated function norm(h::Histogram{T,N}) where {T,N}
    quote
        edges = h.edges
        weights = h.weights
        SumT = norm_type(h)
        v_0 = 1
        s_0 = zero(SumT)
        @nloops(
            $N, i, weights,
            d -> begin
                v_{$N-d+1} = v_{$N-d} * _edge_binvolume(SumT, edges[d], i_d)
                s_{$N-d+1} = zero(SumT)
            end,
            d -> begin
                s_{$N-d} += s_{$N-d+1}
            end,
            begin
                $(Symbol("s_$(N)")) += (@nref $N weights i) * $(Symbol("v_$N"))
            end
        )
        s_0
    end
end


float(h::Histogram{T,N}) where {T<:AbstractFloat,N} = h

float(h::Histogram{T,N}) where {T,N} = Histogram(h.edges, float(h.weights), h.closed, h.isdensity)



"""
    normalize!(h::Histogram{T,N}, aux_weights::Array{T,N}...;
               mode::Symbol=:pdf) where {T<:AbstractFloat,N}

Normalize the histogram `h` and optionally scale one or more auxiliary weight
arrays appropriately. See description of `normalize` for details. Returns `h`.
"""
@generated function normalize!(h::Histogram{T,N}, aux_weights::Array{T,N}...; mode::Symbol=:pdf) where {T<:AbstractFloat,N}
    quote
        edges = h.edges
        weights = h.weights

        for A in aux_weights
            (size(A) != size(weights)) && throw(DimensionMismatch("aux_weights must have same size as histogram weights"))
        end

        if mode == :none
            # nothing to do
        elseif mode == :pdf || mode == :density || mode == :probability
            if h.isdensity
                if mode == :pdf || mode == :probability
                    # histogram already represents a density, just divide weights by norm
                    s = 1/norm(h)
                    weights .*= s
                    for A in aux_weights
                        A .*= s
                    end
                else
                    # :density - histogram already represents a density, nothing to do
                end
            else
                if mode == :pdf || mode == :density
                    # Divide weights by bin volume, for :pdf also divide by sum of weights
                    SumT = norm_type(h)
                    vs_0 = (mode == :pdf) ? sum(SumT, weights) : one(SumT)
                    @nloops $N i weights d->(vs_{$N-d+1} = vs_{$N-d} * _edge_binvolume(SumT, edges[d], i_d)) begin
                        (@nref $N weights i) /= $(Symbol("vs_$N"))
                        for A in aux_weights
                            (@nref $N A i) /= $(Symbol("vs_$N"))
                        end
                    end
                    h.isdensity = true
                else
                    # :probability - divide weights by sum of weights
                    nf = inv(sum(weights))
                    weights .*= nf
                    for A in aux_weights
                        A .*= nf
                    end
                end
            end
        else
            throw(ArgumentError("Normalization mode must be :pdf, :density, :probability or :none"))
        end
        h
    end
end


"""
    normalize(h::Histogram{T,N}; mode::Symbol=:pdf) where {T,N}

Normalize the histogram `h`.

Valid values for `mode` are:

*  `:pdf`: Normalize by sum of weights and bin sizes. Resulting histogram
   has norm 1 and represents a PDF.
* `:density`: Normalize by bin sizes only. Resulting histogram represents
   count density of input and does not have norm 1. Will not modify the
   histogram if it already represents a density (`h.isdensity == 1`).
* `:probability`: Normalize by sum of weights only. Resulting histogram
   represents the fraction of probability mass for each bin and does not have
   norm 1.
*  `:none`: Leaves histogram unchanged. Useful to simplify code that has to
   conditionally apply different modes of normalization.

Successive application of both `:probability` and `:density` normalization (in
any order) is equivalent to `:pdf` normalization.
"""
normalize(h::Histogram{T,N}; mode::Symbol=:pdf) where {T,N} =
    normalize!(deepcopy(float(h)), mode = mode)


"""
    normalize(h::Histogram{T,N}, aux_weights::Array{T,N}...; mode::Symbol=:pdf) where {T,N}

Normalize the histogram `h` and rescales one or more auxiliary weight arrays
at the same time (`aux_weights` may, e.g., contain estimated statistical
uncertainties). The values of the auxiliary arrays are scaled by the same
factor as the corresponding histogram weight values. Returns a tuple of the
normalized histogram and scaled auxiliary weights.
"""
function normalize(h::Histogram{T,N}, aux_weights::Array{T,N}...; mode::Symbol=:pdf) where {T,N}
    h_fltcp = deepcopy(float(h))
    aux_weights_fltcp = map(x -> deepcopy(float(x)), aux_weights)
    normalize!(h_fltcp, aux_weights_fltcp..., mode = mode)
    (h_fltcp, aux_weights_fltcp...)
end


"""
    zero(h::Histogram)

Create a new histogram with the same binning, type and shape of weights
and the same properties (`closed` and `isdensity`) as `h`, with all weights
set to zero.
"""
Base.zero(h::Histogram{T,N,E}) where {T,N,E} =
    Histogram{T,N,E}(deepcopy(h.edges), zero(h.weights), h.closed, h.isdensity)


"""
    merge!(target::Histogram, others::Histogram...)

Update histogram `target` by merging it with the histograms `others`. See
`merge(histogram::Histogram, others::Histogram...)` for details.
"""
function Base.merge!(target::Histogram, others::Histogram...)
    for h in others
        target.edges != h.edges && throw(ArgumentError("can't merge histograms with different binning"))
        size(target.weights) != size(h.weights) && throw(ArgumentError("can't merge histograms with different dimensions"))
        target.closed != h.closed && throw(ArgumentError("can't merge histograms with different closed left/right settings"))
        target.isdensity != h.isdensity && throw(ArgumentError("can't merge histograms with different isdensity settings"))
    end
    for h in others
        target.weights .+= h.weights
    end
    target
end


"""
    merge(h::Histogram, others::Histogram...)

Construct a new histogram by merging `h` with `others`. All histograms must
have the same binning, shape of weights and properties (`closed` and
`isdensity`). The weights of all histograms are summed up for each bin, the
weights of the resulting histogram will have the same type as those of `h`.
"""
Base.merge(h::Histogram, others::Histogram...) = merge!(zero(h), h, others...)

"""
    StatsBase.midpoints(v)

Calculate the midpoints (pairwise mean of consecutive elements).
"""
midpoints(v::AbstractVector) = [middle(v[i - 1], v[i]) for i in 2:length(v)]

midpoints(r::AbstractRange) = r[1:(end - 1)] .+ (step(r) / 2)
