import Base.@deprecate
import Base.depwarn

import Base.varm, Base.stdm
@deprecate varm(v::RealArray, m::Real, wv::AbstractWeights) varm(v, wv, m)
@deprecate varm(A::RealArray, M::RealArray, wv::AbstractWeights, dim::Int) varm(v, wv, m, dim)
@deprecate stdm(v::RealArray, m::Real, wv::AbstractWeights) stdm(v, wv, m)
@deprecate stdm(v::RealArray, m::RealArray, wv::AbstractWeights, dim::Int) stdm(v, wv, m, dim)

@deprecate _moment2(v::RealArray, m::Real, wv::AbstractWeights) _moment2(v, wv, m)
@deprecate _moment3(v::RealArray, m::Real, wv::AbstractWeights) _moment3(v, wv, m)
@deprecate _moment4(v::RealArray, m::Real, wv::AbstractWeights) _moment4(v, wv, m)
@deprecate _momentk(v::RealArray, k::Int, m::Real, wv::AbstractWeights) _momentk(v, k, wv, m)
@deprecate moment(v::RealArray, k::Int, m::Real, wv::AbstractWeights) moment(v, k, wv, m)

@deprecate AIC(obj::StatisticalModel) aic(obj)
@deprecate AICc(obj::StatisticalModel) aicc(obj)
@deprecate BIC(obj::StatisticalModel) bic(obj)

@deprecate R2(obj::StatisticalModel, variant::Symbol) r2(obj, variant)
@deprecate R²(obj::StatisticalModel, variant::Symbol) r²(obj, variant)
@deprecate adjR2(obj::StatisticalModel, variant::Symbol) adjr2(obj, variant)
@deprecate adjR²(obj::StatisticalModel, variant::Symbol) adjr²(obj, variant)

function findat!{T}(r::IntegerArray, a::AbstractArray{T}, b::AbstractArray{T})
    Base.depwarn("findat! is deprecated, use indexin instead", :findat!)
    length(r) == length(b) || raise_dimerror()
    d = indexmap(a)
    @inbounds for i = 1 : length(b)
        r[i] = get(d, b[i], 0)
    end
    return r
end


"""
    findat(a, b)

For each element in `b`, find its first index in `a`. If the value does
not occur in `a`, the corresponding index is 0.
"""
findat(a::AbstractArray, b::AbstractArray) = findat!(Array{Int}(size(b)), a, b)

@deprecate df(obj::StatisticalModel) dof(obj)
@deprecate df_residual(obj::StatisticalModel) dof_residual(obj)

@weights WeightVec

"""
    WeightVec(vs, wsum=sum(vs))

Construct a `WeightVec` with weight values `vs` and sum of weights `wsum`.
"""
function WeightVec{S<:Real, V<:RealVector}(vs::V, s::S=sum(vs))
    new_types = "AnalyticWeights, FrequencyWeights or ProbabilityWeights"
    Base.depwarn("WeightVec is deprecated, use $new_types instead", :WeightVec)
    WeightVec{S, eltype(vs), V}(vs, s)
end

"""
    weights(vs)

Construct a `WeightVec` from a given array.
"""
function weights(vs::RealArray)
    Base.depwarn("weights is deprecated, use aweights, fweights or pweights instead", :weights)
    v = vec(vs)
    s = sum(v)
    WeightVec{typeof(s), eltype(v), typeof(v)}(v, s)
end

"""
    varcorrection(w::WeightVec, corrected=false)

Returns ``\\frac{1}{\sum w}`` when corrected is false and throws an `ArgumentError`
when correction is true.
"""
function varcorrections(w::WeightVec, corrected::Bool=false)
    corrected && throw(ArgumentError("WeightVec does not support bias correction."))
    1 / w.sum
end
