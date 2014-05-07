import Base.push!
import Base.append!
import Base.isequal

## nice-valued ranges for histograms
function histrange{T<:FloatingPoint,N}(v::AbstractArray{T,N}, n::Integer, closed::Symbol)
    if length(v) == 0
        return 0.0:1.0:0.0
    end
    lo, hi = extrema(v)
    if hi == lo
        step = 1.0
    else
        bw = (hi - lo) / n
        e = 10.0^floor(log10(bw))
        r = bw / e
        if r <= 2
            step = 2*e
        elseif r <= 5
            step = 5*e
        else
            step = 10*e
        end
    end
    if closed == :right
        start = step*(ceil(lo/step)-1)
        nm1 = iceil((hi - start)/step)       
    else
        start = step*floor(lo/step)
        nm1 = ifloor((hi - start)/step)+1    
    end
    start:step:(start + nm1*step)
end

function histrange{T<:Integer,N}(v::AbstractArray{T,N}, n::Integer, closed::Symbol)
    if length(v) == 0
        return 0:1:0
    end
    lo, hi = extrema(v)
    if hi == lo
        step = 1
    else
        bw = (hi - lo) / n
        e = int(10^max(0,ifloor(log10(bw))))
        r = bw / e
        if r <= 1
            step = e
        elseif r <= 2
            step = 2*e
        elseif r <= 5
            step = 5*e
        else
            step = 10*e
        end
    end
    if closed == :right
        start = step*(iceil(lo/step)-1)
        nm1 = iceil((hi - start)/step)
    else
        start = step*ifloor(lo/step)
        nm1 = ifloor((hi - start)/step)+1
    end
    start:step:(start + nm1*step)
end

## midpoints of intervals
midpoints(r::Range) = r[1:length(r)-1] + 0.5*step(r)
midpoints(v::AbstractVector) = [0.5*(v[i] + v[i+1]) for i in 1:length(v)-1]

## hist ##
function sturges(n)  # Sturges' formula
    n==0 && return one(n)
    iceil(log2(n))+1
end

# N-dimensional histogram object
immutable Histogram{T<:Real,N,E}
    edges::E
    weights::Array{T,N}
    closed::Symbol
    function Histogram(edges::NTuple{N,AbstractArray},weights::Array{T,N},closed::Symbol)
        closed == :right || closed == :left || error("closed must :left or :right")
        map(x -> length(x)-1,edges) == size(weights) || error("Histogram edge vectors must be 1 longer than corresponding weight dimensions")
        new(edges,weights,closed)
    end
end
Histogram{T,N}(edges::NTuple{N,AbstractVector},weights::AbstractArray{T,N},closed::Symbol=:right) = Histogram{T,N,typeof(edges)}(edges,weights,closed)
Histogram{T,N}(edges::NTuple{N,AbstractVector},::Type{T},closed::Symbol=:right) = Histogram(edges,zeros(T,map(x -> length(x)-1,edges)...),closed)
Histogram{N}(edges::NTuple{N,AbstractVector},closed::Symbol=:right) = Histogram(edges,Int,closed)

isequal(h1::Histogram,h2::Histogram) = isequal(h1.edges,h2.edges) && isequal(h1.weights,h2.weights) && isequal(h1.closed,h2.closed)

# 1-dimensional
Histogram{T}(edge::AbstractVector,weights::AbstractVector{T},closed::Symbol=:right) = Histogram{T,1,(typeof(edge),)}((edge,),weights,closed)
Histogram{T}(edge::AbstractVector,::Type{T},closed::Symbol=:right) = Histogram(edge,zeros(T,length(edge)-1),closed)
Histogram(edge::AbstractVector,closed::Symbol=:right) = Histogram(edge,Int,closed)

function push!{T,E}(h::Histogram{T,1,E}, x::Real)
    i = if h.closed == :right 
        searchsortedfirst(h.edges[1], x) - 1 
    else
        searchsortedlast(h.edges[1], x)
    end
    if 1 <= i <= length(h.weights)
        @inbounds h.weights[i] += one(T)
    end
    h
end

function append!{T}(h::Histogram{T,1}, v::AbstractVector)
    for x in v
        push!(h,x)
    end
    h
end

hist(v::AbstractVector, edg::AbstractVector; closed::Symbol=:right) = append!(Histogram(edg,closed),v)
hist(v::AbstractVector, n::Integer; closed::Symbol=:right) = hist(v,histrange(v,n,closed);closed=closed)
hist(v::AbstractVector; closed::Symbol=:right) = hist(v,sturges(length(v));closed=closed)

# N-dimensional
function push!{T,N}(h::Histogram{T,N},xs::NTuple{N,Real})
    is = if h.closed == :right
        map((edge, x) -> searchsortedfirst(edge,x) - 1, h.edges, xs)
    else
        map(searchsortedlast, h.edges, xs)
    end
    try
        h.weights[is...] += one(T)
    catch e
        isa(e,BoundsError) || rethrow(e)
    end
    h
end
function append!{T,N}(h::Histogram{T,N}, vs::NTuple{N,AbstractVector})
    for xs in zip(vs...)
        push!(h,xs)
    end
    h
end

hist{N}(vs::NTuple{N,AbstractVector}, edges::NTuple{N,AbstractVector}; closed::Symbol=:right) = append!(Histogram(edges,closed),vs)
hist{N}(vs::NTuple{N,AbstractVector}, ns::NTuple{N,Integer}; closed::Symbol=:right) = hist(vs, map((v,n) -> histrange(v,n,closed),vs,ns);closed=closed)
hist{N}(vs::NTuple{N,AbstractVector}, n::Integer; closed::Symbol=:right) = hist(vs, map(v -> histrange(v,n,closed),vs);closed=closed)
hist{N}(vs::NTuple{N,AbstractVector}; closed::Symbol=:right) = hist(vs, sturges(length(vs[1]));closed=closed)
