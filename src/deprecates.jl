using Base: @deprecate, @deprecate_binding, depwarn

if !isdefined(Base, :stderr)
    @deprecate stderr(obj::StatisticalModel) stderror(obj)
else
    function (io::typeof(stderr))(obj::StatisticalModel)
        Base.depwarn("stderr(obj::StatisticalModel) is deprecated, use stderror(obj) instead", :stderr)
        io === stderr ? stderror(obj) : throw(MethodErrror(io, (obj,)))
    end
end

@deprecate model_response(obj::StatisticalModel) response(obj)

@deprecate norepeats(a::AbstractArray) allunique(a)

@deprecate(mad!(v::AbstractArray{<:Real}, center;
                constant::Real = BigFloat("1.482602218505601860547076529360423431326703202590312896536266275245674447622701")),
           mad!(v, center=center, constant=constant))

### Deprecated January 2019
@deprecate scattermatm(x::DenseMatrix, mean, dims::Int) scattermat(x, mean=mean, dims=dims)
@deprecate scattermatm(x::DenseMatrix, mean, wv::AbstractWeights, dims::Int) scattermat(x, wv, mean=mean, dims=dims)
@deprecate scattermat(x::DenseMatrix, dims::Int) scattermat(x, dims=dims)
@deprecate scattermat(x::DenseMatrix, wv::AbstractWeights, dims::Int) scattermat(x, wv, dims=dims)
@deprecate scattermat_zm(x::DenseMatrix, dims::Int) scattermat_zm(x, dims=dims)
@deprecate scattermat_zm(x::DenseMatrix, wv::AbstractWeights, dims::Int) scattermat_zm(x::DenseMatrix, wv::AbstractWeights, dims=dims)
@deprecate mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights, dims::Int) mean!(R, A, w, dims=dims)
@deprecate mean(A::AbstractArray{T}, w::AbstractWeights{W}, dims::Int) where {T<:Number,W<:Real} mean(A, w, dims=dims)

@deprecate wquantile(v::AbstractVector{<:Real}, w::AbstractWeights{<:Real}, p::AbstractVector{<:Real}) quantile(v, w, p)
@deprecate wquantile(v::AbstractVector{<:Real}, w::AbstractWeights{<:Real}, p::Number) quantile(v, w, [p])[1]
@deprecate wquantile(v::AbstractVector{<:Real}, w::AbstractVector{<:Real}, p::AbstractVector{<:Real}) quantile(v, pweights(w), p)
@deprecate wquantile(v::AbstractVector{<:Real}, w::AbstractVector{<:Real}, p::Number) quantile(v, pweights(w), [p])[1]
@deprecate wmedian(v::AbstractVector{<:Real}, w::AbstractWeights{<:Real}) median(v, w)
@deprecate wmedian(v::AbstractVector{<:Real}, w::AbstractVector{<:Real}) median(v, weights(w))

@deprecate quantile(v::AbstractArray{<:Real}) quantile(v, [.0, .25, .5, .75, 1.0])

### Deprecated September 2019
@deprecate sum(A::AbstractArray, w::AbstractWeights, dims::Int) sum(A, w, dims=dims)
@deprecate values(wv::AbstractWeights) convert(Vector, wv)

### Deprecated November 2021
@deprecate stdm(x::AbstractArray{<:Real}, w::AbstractWeights, m::Real; corrected::Union{Bool, Nothing}=nothing) std(x, w, mean=m, corrected=corrected) false
@deprecate varm(x::AbstractArray{<:Real}, w::AbstractWeights, m::Real; corrected::Union{Bool, Nothing}=nothing) var(x, w, mean=m, corrected=corrected) false
@deprecate stdm(x::AbstractArray{<:Real}, w::AbstractWeights, m::AbstractArray{<:Real}, dim::Int; corrected::Union{Bool, Nothing}=nothing) std(x, w, dim, mean=m, corrected=corrected) false
@deprecate varm(x::AbstractArray{<:Real}, w::AbstractWeights, m::AbstractArray{<:Real}, dim::Int; corrected::Union{Bool, Nothing}=nothing) var(x, w, dim, mean=m, corrected=corrected) false
@deprecate varm!(R::AbstractArray, x::AbstractArray{<:Real}, w::AbstractWeights, m::AbstractArray{<:Real}, dim::Int; corrected::Union{Bool, Nothing}=nothing) var!(R, x, w, dim, mean=m, corrected=corrected) false

### This was never part of the public API
### Deprecated April 2024
function make_alias_table!(w::AbstractVector, wsum,
                           a::AbstractVector{Float64},
                           alias::AbstractVector{Int})
    Base.depwarn("make_alias_table! is both internal and deprecated, use AliasTables.jl instead", :make_alias_table!)
    # Arguments:
    #
    #   w [in]:         input weights
    #   wsum [in]:      pre-computed sum(w)
    #
    #   a [out]:        acceptance probabilities
    #   alias [out]:    alias table
    #
    # Note: a and w can be the same array, then that array will be
    #       overwritten inplace by acceptance probabilities
    #
    # Returns nothing
    #

    n = length(w)
    length(a) == length(alias) == n ||
        throw(DimensionMismatch("Inconsistent array lengths."))

    ac = n / wsum
    for i = 1:n
        @inbounds a[i] = w[i] * ac
    end

    larges = Vector{Int}(undef, n)
    smalls = Vector{Int}(undef, n)
    kl = 0  # actual number of larges
    ks = 0  # actual number of smalls

    for i = 1:n
        @inbounds ai = a[i]
        if ai > 1.0
            larges[kl+=1] = i  # push to larges
        elseif ai < 1.0
            smalls[ks+=1] = i  # push to smalls
        end
    end

    while kl > 0 && ks > 0
        s = smalls[ks]; ks -= 1  # pop from smalls
        l = larges[kl]; kl -= 1  # pop from larges
        @inbounds alias[s] = l
        @inbounds al = a[l] = (a[l] - 1.0) + a[s]
        if al > 1.0
            larges[kl+=1] = l  # push to larges
        else
            smalls[ks+=1] = l  # push to smalls
        end
    end

    # this loop should be redundant, except for rounding
    for i = 1:ks
        @inbounds a[smalls[i]] = 1.0
    end
    nothing
end
