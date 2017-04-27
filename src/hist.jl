using Base.Cartesian

import Base: show, ==, push!, append!, float, norm, normalize, normalize!

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


## Fast getindex function for multiple arrays, returns a tuple of array elements
@inline Base.@propagate_inbounds @generated function _multi_getindex(i::Integer, c::AbstractArray...)
    N = length(c)
    result_expr = Expr(:tuple)
    for j in 1:N
        push!(result_expr.args, :(c[$j][i]))
    end
    result_expr
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
    isdensity::Bool
    function (::Type{Histogram{T,N,E}}){T,N,E}(edges::NTuple{N,AbstractArray},
                                               weights::Array{T,N}, closed::Symbol, isdensity::Bool=false)
        closed == :right || closed == :left || error("closed must :left or :right")
        isdensity && !(T <: AbstractFloat) && error("Density histogram must have float-type weights")
        map(x -> length(x)-1,edges) == size(weights) || error("Histogram edge vectors must be 1 longer than corresponding weight dimensions")
        new{T,N,E}(edges,weights,closed,isdensity)
    end
end

Histogram{T,N}(edges::NTuple{N,AbstractVector},weights::AbstractArray{T,N},closed::Symbol=:default_left, isdensity::Bool=false) =
    Histogram{T,N,typeof(edges)}(edges,weights,_check_closed_arg(closed,:Histogram),isdensity)

Histogram{T,N}(edges::NTuple{N,AbstractVector},::Type{T},closed::Symbol=:default_left, isdensity::Bool=false) =
    Histogram(edges,zeros(T,map(x -> length(x)-1,edges)...),_check_closed_arg(closed,:Histogram),isdensity)

Histogram{N}(edges::NTuple{N,AbstractVector},closed::Symbol=:default_left, isdensity::Bool=false) =
    Histogram(edges,Int,_check_closed_arg(closed,:Histogram),isdensity)

function show(io::IO, h::AbstractHistogram)
    println(io, typeof(h))
    println(io,"edges:")
    for e in h.edges
        println(io,"  ",e)
    end
    println(io,"weights: ",h.weights)
    println(io,"closed: ",h.closed)
    print(io,"isdensity: ",h.isdensity)
end

(==)(h1::Histogram,h2::Histogram) = (==)(h1.edges,h2.edges) && (==)(h1.weights,h2.weights) && (==)(h1.closed,h2.closed) && (==)(h1.isdensity,h2.isdensity)


binindex{T,E}(h::AbstractHistogram{T,1,E}, x::Real) = binindex(h, (x,))[1]

binindex{T,N,E}(h::Histogram{T,N,E}, xs::NTuple{N,Real}) =
    map((edge, x) -> _edge_binindex(edge, h.closed, x), h.edges, xs)

@inline function _edge_binindex(edge::AbstractVector, closed::Symbol, x::Real)
    if closed == :right
        searchsortedfirst(edge, x) - 1
    else
        searchsortedlast(edge, x)
    end
end


binvolume{T,E}(h::AbstractHistogram{T,1,E}, binidx::Integer) = binvolume(h, (binidx,))
binvolume{V,T,E}(::Type{V}, h::AbstractHistogram{T,1,E}, binidx::Integer) = binvolume(V, h, (binidx,))

binvolume{T,N,E}(h::Histogram{T,N,E}, binidx::NTuple{N,Integer}) =
    binvolume(promote_type(map(eltype, h.edges)...), h, binidx)

binvolume{V,T,N,E}(::Type{V}, h::Histogram{T,N,E}, binidx::NTuple{N,Integer}) =
    prod(map((edge, i) -> _edge_binvolume(V, edge, i), h.edges, binidx))

@inline _edge_binvolume{V}(::Type{V}, edge::AbstractVector, i::Integer) = V(edge[i+1]) - V(edge[i])
@inline _edge_binvolume{V}(::Type{V}, edge::Range, i::Integer) = V(step(edge))
@inline _edge_binvolume(edge::AbstractVector, i::Integer) = _edge_binvolume(eltype(edge), edge, i)


# 1-dimensional

Histogram{T}(edge::AbstractVector, weights::AbstractVector{T}, closed::Symbol=:default_left, isdensity::Bool=false) =
    Histogram((edge,), weights, closed, isdensity)

Histogram{T}(edge::AbstractVector, ::Type{T}, closed::Symbol=:default_left, isdensity::Bool=false) =
    Histogram((edge,), T, closed, isdensity)

Histogram(edge::AbstractVector, closed::Symbol=:default_left, isdensity::Bool=false) =
    Histogram((edge,), closed, isdensity)


push!{T,E}(h::AbstractHistogram{T,1,E}, x::Real, w::Real) = push!(h, (x,), w)
push!{T,E}(h::AbstractHistogram{T,1,E}, x::Real) = push!(h,x,one(T))
append!{T}(h::AbstractHistogram{T,1}, v::AbstractVector) = append!(h, (v,))
append!{T}(h::AbstractHistogram{T,1}, v::AbstractVector, wv::Union{AbstractVector,WeightVec}) = append!(h, (v,), wv)


fit{T}(::Type{Histogram{T}},v::AbstractVector, edg::AbstractVector; closed::Symbol=:default_left) =
    fit(Histogram{T},(v,), (edg,), closed=closed)
fit{T}(::Type{Histogram{T}},v::AbstractVector; closed::Symbol=:default_left, nbins=sturges(length(v))) =
    fit(Histogram{T},(v,); closed=closed, nbins=nbins)
fit{T}(::Type{Histogram{T}},v::AbstractVector, wv::WeightVec, edg::AbstractVector; closed::Symbol=:default_left) =
    fit(Histogram{T},(v,), wv, (edg,), closed=closed)
fit{T}(::Type{Histogram{T}},v::AbstractVector, wv::WeightVec; closed::Symbol=:default_left, nbins=sturges(length(v))) =
    fit(Histogram{T}, (v,), wv; closed=closed, nbins=nbins)

fit{W}(::Type{Histogram}, v::AbstractVector, wv::WeightVec{W}, args...; kwargs...) = fit(Histogram{W}, v, wv, args...; kwargs...)


# N-dimensional

function push!{T,N}(h::Histogram{T,N},xs::NTuple{N,Real},w::Real)
    (h.isdensity == true) && error("Density histogram must have float-type weights")
    idx = binindex(h, xs)
    if checkbounds(Bool, h.weights, idx...)
        @inbounds h.weights[idx...] += w
    end
    h
end

function push!{T<:AbstractFloat,N}(h::Histogram{T,N},xs::NTuple{N,Real},w::Real)
    idx = binindex(h, xs)
    if checkbounds(Bool, h.weights, idx...)
        @inbounds h.weights[idx...] += h.isdensity ? w / binvolume(h, idx) : w
    end
    h
end

push!{T,N}(h::AbstractHistogram{T,N},xs::NTuple{N,Real}) = push!(h,xs,one(T))


function append!{T,N}(h::AbstractHistogram{T,N}, vs::NTuple{N,AbstractVector})
    @inbounds for i in eachindex(vs...)
        xs = _multi_getindex(i, vs...)
        push!(h, xs, one(T))
    end
    h
end
function append!{T,N}(h::AbstractHistogram{T,N}, vs::NTuple{N,AbstractVector}, wv::AbstractVector)
    @inbounds for i in eachindex(wv, vs...)
        xs = _multi_getindex(i, vs...)
        push!(h, xs, wv[i])
    end
    h
end
append!{T,N}(h::AbstractHistogram{T,N}, vs::NTuple{N,AbstractVector}, wv::WeightVec) = append!(h, vs, values(wv))


fit{T,N}(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}, edges::NTuple{N,AbstractVector}; closed::Symbol=:default_left) =
    append!(Histogram(edges, T, _check_closed_arg(closed,:fit), false), vs)

fit{T,N}(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}; closed::Symbol=:default_left, isdensity::Bool=false, nbins=sturges(length(vs[1]))) = begin
    closed = _check_closed_arg(closed,:fit)
    fit(Histogram{T}, vs, histrange(vs,nbins,closed); closed=closed)
end

fit{T,N,W}(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}, wv::WeightVec{W}, edges::NTuple{N,AbstractVector}; closed::Symbol=:default_left) =
    append!(Histogram(edges, T, _check_closed_arg(closed,:fit), false), vs, wv)

fit{T,N}(::Type{Histogram{T}}, vs::NTuple{N,AbstractVector}, wv::WeightVec; closed::Symbol=:default_left, isdensity::Bool=false, nbins=sturges(length(vs[1]))) = begin
    closed = _check_closed_arg(closed,:fit)
    fit(Histogram{T}, vs, wv, histrange(vs,nbins,closed); closed=closed)
end

fit(::Type{Histogram}, args...; kwargs...) = fit(Histogram{Int}, args...; kwargs...)
fit{N,W}(::Type{Histogram}, vs::NTuple{N,AbstractVector}, wv::WeightVec{W}, args...; kwargs...) = fit(Histogram{W}, vs, wv, args...; kwargs...)


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


"""
    norm(h::Histogram)

Calculate the norm of histogram `h` as the absolute value of its integral.
"""
@generated function norm{T, N, E}(h::Histogram{T, N, E})
    quote
        edges = h.edges
        weights = h.weights
        SumT = norm_type(h)
        v_0 = 1
        s_0 = zero(SumT)
        @inbounds @nloops(
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


float{T<:AbstractFloat, N, E}(h::Histogram{T, N, E}) = h

float{T, N, E}(h::Histogram{T, N, E}) = Histogram(h.edges, float(h.weights), h.closed, h.isdensity)



"""
    normalize!{T<:AbstractFloat, N, E}(h::Histogram{T, N, E}, aux_weights::Array{T,N}...; mode::Symbol = :pdf)

Normalize the histogram `h` and optionally scale one or more auxiliary weight
arrays appropriately. See description of `normalize` for details. Returns `h`.
"""
@generated function normalize!{T<:AbstractFloat, N, E}(h::Histogram{T, N, E}, aux_weights::Array{T,N}...; mode::Symbol = :pdf)
    quote
        edges = h.edges
        weights = h.weights

        for A in aux_weights
            (size(A) != size(weights)) && throw(DimensionMismatch("aux_weights must have same size as histogram weights"))
        end

        if mode == :none
            # nothing to do
        elseif mode == :pdf || mode == :density
            if h.isdensity
                if mode == :pdf
                    # histogram already represents a density, just divide weights by norm
                    s = 1/norm(h)
                    weights .*= s
                    for A in aux_weights
                        A .*= s
                    end
                else
                    # histogram already represents a density, nothing to do
                end
            else
                # Divide weights by bin volume, for :pdf also divide by sum of weights
                SumT = norm_type(h)
                vs_0 = (mode == :pdf) ? sum(SumT(x) for x in weights) : one(SumT)
                @inbounds @nloops $N i weights d->(vs_{$N-d+1} = vs_{$N-d} * _edge_binvolume(SumT, edges[d], i_d)) begin
                    (@nref $N weights i) /= $(Symbol("vs_$N"))
                    for A in aux_weights
                        (@nref $N A i) /= $(Symbol("vs_$N"))
                    end
                end
            end
            h.isdensity = true
        else mode != :pdf && mode != :density
            throw(ArgumentError("Normalization mode must be :pdf, :density or :none"))
        end
        h
    end
end


"""
    normalize{T, N, E}(h::Histogram{T, N, E}; mode::Symbol = :pdf)

Normalize the histogram `h`.

Valid values for `mode` are:

*  `:pdf`: Normalize by sum of weights and bin sizes. Resulting histogram
   has norm 1 and represents a PDF.
* `:density`: Normalize by bin sizes only. Resulting histogram represents
   count density of input and does not have norm 1. Will not modify the
   histogram if it already represents a density (`h.isdensity == 1`).
*  `:none`: Leaves histogram unchanged. Useful to simplify code that has to
   conditionally apply different modes of normalization.
"""
normalize{T, N, E}(h::Histogram{T, N, E}; mode::Symbol = :pdf) =
    normalize!(deepcopy(float(h)), mode = mode)


"""
    normalize{T, N, E}(h::Histogram{T, N, E}, aux_weights::Array{T,N}...; mode::Symbol = :pdf)

Normalize the histogram `h` and rescales one or more auxiliary weight arrays
at the same time (`aux_weights` may, e.g., contain estimated statistical
uncertainties). The values of the auxiliary arrays are scaled by the same
factor as the corresponding histogram weight values. Returns a tuple of the
normalized histogram and scaled auxiliary weights.
"""
function normalize{T, N, E}(h::Histogram{T, N, E}, aux_weights::Array{T,N}...; mode::Symbol = :pdf)
    h_fltcp = deepcopy(float(h))
    aux_weights_fltcp = map(x -> deepcopy(float(x)), aux_weights)
    normalize!(h_fltcp, aux_weights_fltcp..., mode = mode)
    (h_fltcp, aux_weights_fltcp...)
end
