# Empirical estimation of CDF and PDF

## Empirical CDF

struct ECDF{T <: Real, W <: Real, I}
    sorted_values::Vector{Tuple{T, W, W, W}}
end

weighttype(::Type{<:ECDF{<:Any, W}}) where W = W
weighttype(ecdf::ECDF) = weighttype(typeof(ecdf))
isinterpolating(::Type{<:ECDF{<:Any, <:Any, I}}) where I = I
isinterpolating(ecdf::ECDF) = isinterpolating(typeof(ecdf))

function Base.show(io::IO, e::ECDF)
    print(io, typeof(e))
    print(io, "(")
    print(io, length(e.sorted_values))
    println(io, " values)")
end

function (ecdf::ECDF)(x::Real)
    isnan(x) && return NaN
    pos = searchsortedlast(ecdf.sorted_values, x, by=first)
    (pos == 0) && return zero(weighttype(ecdf))
    @inbounds val, cdf_val, val_invdelta, cdf_delta = ecdf.sorted_values[pos]
    if !isinterpolating(ecdf) || (val == x) || (pos == length(ecdf.sorted_values))
        return cdf_val
    else
        return muladd(cdf_delta, min((x - val) * val_invdelta, one(weighttype(ecdf))), cdf_val)
    end
end

# broadcasts ecdf() over an array
# caching the previous calculated value
function Base.Broadcast.broadcasted(ecdf::ECDF, v::AbstractArray)
    res = similar(v, weighttype(ecdf))
    @inbounds for i in eachindex(v)
        res[i] = i == 1 || v[i] != v[i-1] ? ecdf(v[i]) : res[i-1]
    end
    return res
end

"""
    ecdf(X; weights::AbstractWeights)

Return an empirical cumulative distribution function (ECDF) based on a vector of samples
given in `X`. Optionally providing `weights` returns a weighted ECDF.

Note: this function that returns a callable composite type, which can then be applied to
evaluate CDF values on other samples.

`extrema`, `minimum`, and `maximum` are supported to for obtaining the range over which
function is inside the interval ``(0,1)``; the function is defined for the whole real line.
"""
function ecdf(X::RealVector;
              weights::Union{Nothing, RealVector}=nothing,
              interpolate::Bool=false)
    any(isnan, X) && throw(ArgumentError("ecdf can not include NaN values"))
    evenweights = isnothing(weights) || isempty(weights)
    evenweights || (length(X) == length(weights)) ||
        throw(ArgumentError("data and weight vectors must be the same size, " *
                            "got $(length(X)) and $(length(weights))"))
    T = eltype(X)
    W0 = evenweights ? Int : eltype(weights)
    W = isnothing(weights) ? Float64 : eltype(one(W0)/sum(weights))

    wsum = evenweights ? length(X) : sum(weights)
    ord = sortperm(X)

    sorted_vals = sizehint!(Vector{Tuple{T, W, W, W}}(), length(X))
    isempty(X) && return ECDF{T, W, interpolate}(sorted_vals)

    valprev = val = X[first(ord)]
    wsumprev = zero(W0)
    valw = zero(W0)

    push_valprev!() = push!(sorted_vals, (valprev, min(wsumprev/wsum, one(W)),
                                          inv(val - valprev), valw/wsum))

    @inbounds for i in ord
        valnew = X[i]
        if (val != valnew) || (i == last(ord))
            (wsumprev > 0) && push_valprev!()
            valprev = val
            val = valnew
            wsumprev += valw
            valw = zero(W0)
        end
        valw += evenweights ? one(W0) : weights[i]
    end
    #@assert valw + wsumprev == wsum # may fail due to fp-arithmetic
    (wsumprev > 0) && push_valprev!()
    # last value
    push!(sorted_vals, (val, one(W), zero(W), zero(W)))
    return ECDF{T,W,interpolate}(sorted_vals)
end

minimum(ecdf::ECDF) = first(ecdf.sorted_values)[1]

maximum(ecdf::ECDF) = last(ecdf.sorted_values)[1]

extrema(ecdf::ECDF) = (minimum(ecdf), maximum(ecdf))
