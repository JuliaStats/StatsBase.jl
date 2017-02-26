using Base.Cartesian

import Base: show, ==, push!, append!, norm, normalize, normalize!

# Mechanism for temporary deprecation of default for "closed" (because default
# value has changed). After deprecation is lifed, remove "_check_closed_arg"
# and all calls to it, and replace every ":default_left" with ":left". Also
# remove "closed=:left" in tests for all lines marked "FIXME: closed".
function _check_closed_arg(closed::Symbol, funcsym)
    if closed == :default_left
        Base.depwarn("Default for keyword argument \"closed\" has changed from :right to :left.", funcsym)
        :left
    else
        closed
    end
end


## nice-valued ranges for histograms
function histrange{T}(v::AbstractArray{T}, n::Integer, closed::Symbol=:default_left)
    closed = _check_closed_arg(closed,:histrange)
    F = float(T)
    nv = length(v)
    if nv == 0 && n < 0
        throw(ArgumentError("number of bins must be ≥ 0 for an empty array, got $n"))
    elseif nv > 0 && n < 1
        throw(ArgumentError("number of bins must be ≥ 1 for a non-empty array, got $n"))
    elseif nv == 0
        return zero(F):zero(F)
    end

    lo, hi = extrema(v)
    histrange(F(lo), F(hi), n, closed)
end

function histrange{F}(lo::F, hi::F, n::Integer, closed::Symbol=:default_left)
    closed = _check_closed_arg(closed,:histrange)
    if hi == lo
        start = F(hi)
        step = one(F)
        divisor = one(F)
        len = one(F)
    else
        bw = (F(hi) - F(lo)) / n
        lbw = log10(bw)
        if lbw >= 0
            step = exp10(floor(lbw))
            r = bw / step
            if r <= 1.1
                nothing
            elseif r <= 2.2
                step *= 2
            elseif r <= 5.5
                step *= 5
            else
                step *= 10
            end
            divisor = one(F)
            start = step*floor(lo/step)
            len = ceil((hi - start)/step)
        else
            divisor = exp10(-floor(lbw))
            r = bw * divisor
            if r <= 1.1
                nothing
            elseif r <= 2.2
                divisor /= 2
            elseif r <= 5.5
                divisor /= 5
            else
                divisor /= 10
            end
            step = one(F)
            start = floor(lo*divisor)
            len = ceil(hi*divisor - start)
        end
    end
    # fix up endpoints
    if closed == :right #(,]
        while lo <= start/divisor
            start -= step
        end
        while (start + (len-1)*step)/divisor < hi
            len += one(F)
        end
    else
        while lo < start/divisor
            start -= step
        end
        while (start + (len-1)*step)/divisor <= hi
            len += one(F)
        end
    end
    @static if VERSION < v"0.6.0-dev.2376" # Julia PR 18777
        FloatRange(start,step,len,divisor)
    else
        Base.floatrange(start,step,len,divisor)
    end
end

histrange{N}(vs::NTuple{N,AbstractVector},nbins::NTuple{N,Integer},closed::Symbol) =
    map((v,n) -> histrange(v,n,closed),vs,nbins)
histrange{N}(vs::NTuple{N,AbstractVector},nbins::Integer,closed::Symbol) =
    map(v -> histrange(v,nbins,closed),vs)



## histograms ##
function sturges(n)  # Sturges' formula
    n==0 && return one(n)
    ceil(Integer, log2(n))+1
end

@compat abstract type AbstractHistogram{T<:Real,N,E} end

# N-dimensional histogram object
type Histogram{T<:Real,N,E} <: AbstractHistogram{T,N,E}
    edges::E
    weights::Array{T,N}
    closed::Symbol
    function (::Type{Histogram{T,N,E}}){T,N,E}(edges::NTuple{N,AbstractArray},
                                               weights::Array{T,N}, closed::Symbol)
        closed == :right || closed == :left || error("closed must :left or :right")
        map(x -> length(x)-1,edges) == size(weights) || error("Histogram edge vectors must be 1 longer than corresponding weight dimensions")
        new{T,N,E}(edges,weights,closed)
    end
end

Histogram{T,N}(edges::NTuple{N,AbstractVector},weights::AbstractArray{T,N},closed::Symbol=:default_left) =
    Histogram{T,N,typeof(edges)}(edges,weights,_check_closed_arg(closed,:Histogram))

Histogram{T,N}(edges::NTuple{N,AbstractVector},::Type{T},closed::Symbol=:default_left) =
    Histogram(edges,zeros(T,map(x -> length(x)-1,edges)...),_check_closed_arg(closed,:Histogram))

Histogram{N}(edges::NTuple{N,AbstractVector},closed::Symbol=:default_left) =
    Histogram(edges,Int,_check_closed_arg(closed,:Histogram))

function show(io::IO, h::AbstractHistogram)
    println(io, typeof(h))
    println(io,"edges:")
    for e in h.edges
        println(io,"  ",e)
    end
    println(io,"weights: ",h.weights)
    print(io,"closed: ",h.closed)
end

(==)(h1::Histogram,h2::Histogram) = (==)(h1.edges,h2.edges) && (==)(h1.weights,h2.weights) && (==)(h1.closed,h2.closed)

# 1-dimensional
Histogram{T}(edge::AbstractVector, weights::AbstractVector{T}, closed::Symbol=:default_left) =
    Histogram((edge,), weights, _check_closed_arg(closed,:Histogram))

Histogram{T}(edge::AbstractVector, ::Type{T}, closed::Symbol=:default_left) =
    Histogram(edge, zeros(T,length(edge)-1), _check_closed_arg(closed,:Histogram))

Histogram(edge::AbstractVector,closed::Symbol=:default_left) =
    Histogram(edge, Int, _check_closed_arg(closed,:Histogram))

function push!{T,E}(h::Histogram{T,1,E}, x::Real,w::Real)
    i = if h.closed == :right
        searchsortedfirst(h.edges[1], x) - 1
    else
        searchsortedlast(h.edges[1], x)
    end
    if 1 <= i <= length(h.weights)
        @inbounds h.weights[i] += w
    end
    h
end
push!{T,E}(h::AbstractHistogram{T,1,E}, x::Real) = push!(h,x,one(T))

function append!{T}(h::AbstractHistogram{T,1}, v::AbstractVector)
    for x in v
        push!(h,x)
    end
    h
end
function append!{T}(h::AbstractHistogram{T,1}, v::AbstractVector,wv::WeightVec)
    for (x,w) in zip(v,wv.values)
        push!(h,x,w)
    end
    h
end

fit(::Type{Histogram},v::AbstractVector, edg::AbstractVector; closed::Symbol=:default_left) =
    append!(Histogram(edg,_check_closed_arg(closed,:fit)), v)
fit(::Type{Histogram},v::AbstractVector; closed::Symbol=:default_left, nbins=sturges(length(v))) = begin
    closed = _check_closed_arg(closed,:fit)
    fit(Histogram, v, histrange(v,nbins,closed); closed=closed)
end

fit{W}(::Type{Histogram},v::AbstractVector, wv::WeightVec{W}, edg::AbstractVector; closed::Symbol=:default_left) =
    append!(Histogram(edg,W,_check_closed_arg(closed,:fit)), v, wv)
fit(::Type{Histogram},v::AbstractVector, wv::WeightVec; closed::Symbol=:default_left, nbins=sturges(length(v))) = begin
    closed = _check_closed_arg(closed,:fit)
    fit(Histogram, v, wv, histrange(v,nbins,closed); closed=closed)
end

# N-dimensional
function push!{T,N}(h::Histogram{T,N},xs::NTuple{N,Real},w::Real)
    is = if h.closed == :right
        map((edge, x) -> searchsortedfirst(edge,x) - 1, h.edges, xs)
    else
        map(searchsortedlast, h.edges, xs)
    end
    try
        h.weights[is...] += w
    catch e
        isa(e,BoundsError) || rethrow(e)
    end
    h
end
push!{T,N}(h::AbstractHistogram{T,N},xs::NTuple{N,Real}) = push!(h,xs,one(T))

function append!{T,N}(h::AbstractHistogram{T,N}, vs::NTuple{N,AbstractVector})
    for xs in zip(vs...)
        push!(h,xs)
    end
    h
end
function append!{T,N}(h::AbstractHistogram{T,N}, vs::NTuple{N,AbstractVector},wv::WeightVec)
    for (xs,w) in zip(zip(vs...),wv.values)
        push!(h,xs,w)
    end
    h
end

fit{N}(::Type{Histogram}, vs::NTuple{N,AbstractVector}, edges::NTuple{N,AbstractVector}; closed::Symbol=:default_left) =
    append!(Histogram(edges,_check_closed_arg(closed,:fit)), vs)
fit{N}(::Type{Histogram}, vs::NTuple{N,AbstractVector}; closed::Symbol=:default_left, nbins=sturges(length(vs[1]))) = begin
    closed = _check_closed_arg(closed,:fit)
    fit(Histogram, vs, histrange(vs,nbins,closed); closed=closed)
end

fit{N,W}(::Type{Histogram}, vs::NTuple{N,AbstractVector}, wv::WeightVec{W}, edges::NTuple{N,AbstractVector}; closed::Symbol=:default_left) =
    append!(Histogram(edges,W,_check_closed_arg(closed,:fit)), vs, wv)
fit{N}(::Type{Histogram},vs::NTuple{N,AbstractVector}, wv::WeightVec; closed::Symbol=:default_left, nbins=sturges(length(vs[1]))) = begin
    closed = _check_closed_arg(closed,:fit)
    fit(Histogram, vs, wv, histrange(vs,nbins,closed); closed=closed)
end


# Get a suitable high-precision type for the norm of a histogram.
@generated function norm_type{T, N, E}(h::Histogram{T, N, E})
    args = [:( eltype(edges[$d]) ) for d = 1:N]
    quote
        edges = h.edges
        norm_type(promote_type(T, $(args...)))
    end
end

norm_type{T<:Integer}(::Type{T}) = promote_type(T, Int64)
norm_type{T<:AbstractFloat}(::Type{T}) = promote_type(T, Float64)


@generated function norm{T, N, E}(h::Histogram{T, N, E})
    quote
        edges = h.edges
        weights = h.weights
        S = norm_type(h)
        v_0 = 1
        s_0 = zero(S)
        @inbounds @nloops(
            $N, i, weights,
            d -> begin
                v_{$N-d+1} = v_{$N-d} * (edges[d][i_d + 1] - edges[d][i_d])
                s_{$N-d+1} = zero(S)
            end,
            d -> begin
                s_{$N-d} += s_{$N-d+1}
            end,
            begin
                $(Symbol("s_$(N)")) += (@nref $N weights i) * $(Symbol("v_$N"))
            end
        )
        norm(s_0)
    end
end


# Deep copy of h, ensuring result has floating point weights
function _deepcopy_fltweights{T, N, E}(h::Histogram{T, N, E})
    weights_flt = float(h.weights)
    weights_fltcp = is(weights_flt, h.weights) ? deepcopy(weights_flt) : weights_flt
    Histogram{eltype(weights_fltcp), N, E}(deepcopy(h.edges), weights_fltcp, h.closed)
end


# Normalization modes:
#
# * `:norm`: Normalize by norm(h). Resulting histogram has norm 1. Histograms
#    of same data but with different binning will have different weight values
#    at same coordinates after normalization. Operation is idempotent.
# * `:pdf`: Normalize by sum of weights and bin sizes. Resulting histogram
#    has norm 1 and represents a PDF. Histograms of same data but with
#    different binning will have same weight values at same coordinates
#    after normalization. Operation is *not* idempotent.
# * `:density`: Normalize by bin sizes only. Resulting histogram represents
#    count density of input and does not have norm 1. Histograms of same data
#    but with different binning will have same weight values at same
#    coordinates after normalization. Operation is *not* idempotent.

@generated function normalize!{T, N, E}(h::Histogram{T, N, E}; mode::Symbol = :norm)
    quote
        edges = h.edges
        weights = h.weights

        if mode == :norm
            weights ./= norm(h)
        elseif mode == :pdf || mode == :density
            SumT = promote_type($T, Float64)
            vs_0 = (mode == :pdf) ? sum(SumT(x) for x in weights) : one(SumT)
            @inbounds @nloops $N i weights d->(vs_{$N-d+1} = vs_{$N-d} * (edges[d][i_d + 1] - edges[d][i_d])) begin
                (@nref $N weights i) /= $(Symbol("vs_$N"))
            end
        else
            error("mode must be :norm, :pdf or :density")
        end
        h
    end
end


normalize(h::Histogram; mode::Symbol = :norm) =
    normalize!(_deepcopy_fltweights(h), mode = mode)
