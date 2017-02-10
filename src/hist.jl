import Base: show, ==, push!, append!

## nice-valued ranges for histograms
function histrange{T}(v::AbstractArray{T}, n::Integer, closed::Symbol=:right)
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

function histrange{F}(lo::F, hi::F, n::Integer, closed::Symbol=:right)
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

abstract AbstractHistogram{T<:Real,N,E}

# N-dimensional histogram object
type Histogram{T<:Real,N,E} <: AbstractHistogram{T,N,E}
    edges::E
    weights::Array{T,N}
    closed::Symbol
    function Histogram(edges::NTuple{N,AbstractArray},weights::Array{T,N},closed::Symbol)
        closed == :right || closed == :left || error("closed must :left or :right")
        map(x -> length(x)-1,edges) == size(weights) || error("Histogram edge vectors must be 1 longer than corresponding weight dimensions")
        new(edges,weights,closed)
    end
end

Histogram{T,N}(edges::NTuple{N,AbstractVector},weights::AbstractArray{T,N},closed::Symbol=:right) =
    Histogram{T,N,typeof(edges)}(edges,weights,closed)

Histogram{T,N}(edges::NTuple{N,AbstractVector},::Type{T},closed::Symbol=:right) =
    Histogram(edges,zeros(T,map(x -> length(x)-1,edges)...),closed)

Histogram{N}(edges::NTuple{N,AbstractVector},closed::Symbol=:right) =
    Histogram(edges,Int,closed)

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
Histogram{T}(edge::AbstractVector, weights::AbstractVector{T}, closed::Symbol=:right) =
    Histogram((edge,), weights, closed)

Histogram{T}(edge::AbstractVector, ::Type{T}, closed::Symbol=:right) =
    Histogram(edge, zeros(T,length(edge)-1), closed)

Histogram(edge::AbstractVector,closed::Symbol=:right) =
    Histogram(edge, Int, closed)

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

fit(::Type{Histogram},v::AbstractVector, edg::AbstractVector; closed::Symbol=:right) =
    append!(Histogram(edg,closed), v)
fit(::Type{Histogram},v::AbstractVector; closed::Symbol=:right, nbins=sturges(length(v))) =
    fit(Histogram, v, histrange(v,nbins,closed); closed=closed)

fit{W}(::Type{Histogram},v::AbstractVector, wv::WeightVec{W}, edg::AbstractVector; closed::Symbol=:right) =
    append!(Histogram(edg,W,closed), v, wv)
fit(::Type{Histogram},v::AbstractVector, wv::WeightVec; closed::Symbol=:right, nbins=sturges(length(v))) =
    fit(Histogram, v, wv, histrange(v,nbins,closed); closed=closed)

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

fit{N}(::Type{Histogram}, vs::NTuple{N,AbstractVector}, edges::NTuple{N,AbstractVector}; closed::Symbol=:right) =
    append!(Histogram(edges,closed), vs)
fit{N}(::Type{Histogram}, vs::NTuple{N,AbstractVector}; closed::Symbol=:right, nbins=sturges(length(vs[1]))) =
    fit(Histogram, vs, histrange(vs,nbins,closed); closed=closed)

fit{N,W}(::Type{Histogram}, vs::NTuple{N,AbstractVector}, wv::WeightVec{W}, edges::NTuple{N,AbstractVector}; closed::Symbol=:right) =
    append!(Histogram(edges,W,closed), vs, wv)
fit{N}(::Type{Histogram},vs::NTuple{N,AbstractVector}, wv::WeightVec; closed::Symbol=:right, nbins=sturges(length(vs[1]))) =
    fit(Histogram, vs, wv, histrange(vs,nbins,closed); closed=closed)
