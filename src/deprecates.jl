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

_correction_dep_msg(fname) =
    string(fname, " will default to `corrected=true` in the future.")

# The following methods are for wrapping the deprecated `correction=false` behaviour.
# When we default to `correction=true` these methods should be removed in favour of
# adding `corrected::Bool=true` in the appropriate methods.

function Base.varm(v::RealArray, wv::AbstractWeights, m::Real)
    Base.depwarn(_correction_dep_msg("`varm`"), :varm)
    varm(v, wv, m, false)
end

function Base.varm(A::RealArray, wv::AbstractWeights, M::RealArray, dim::Int)
    Base.depwarn(_correction_dep_msg("`varm`"), :varm)
    varm(A, wv, M, dim, false)
end

function Base.var(v::RealArray, wv::AbstractWeights; mean=nothing)
    Base.depwarn(_correction_dep_msg("`var`"), :var)
    var(v, wv, false; mean=mean)
end

function Base.var(A::RealArray, wv::AbstractWeights, dim::Int; mean=nothing)
    Base.depwarn(_correction_dep_msg("`var`"), :var)
    var(A, wv, dim, false; mean=mean)
end

function Base.stdm(v::RealArray, wv::AbstractWeights, m::Real)
    Base.depwarn(_correction_dep_msg("`stdm`"), :stdm)
    stdm(v, wv, m, false)
end

function Base.stdm(A::RealArray, wv::AbstractWeights, M::RealArray, dim::Int)
    Base.depwarn(_correction_dep_msg("`stdm`"), :stdm)
    stdm(A, wv, M, dim, false)
end

function Base.std(v::RealArray, wv::AbstractWeights; mean=nothing)
    Base.depwarn(_correction_dep_msg("`std`"), :std)
    std(v, wv, false; mean=mean)
end

function Base.std(A::RealArray, wv::AbstractWeights, dim::Int; mean=nothing)
    Base.depwarn(_correction_dep_msg("`std`"), :std)
    std(A, wv, dim, false; mean=mean)
end

function mean_and_var(A::RealArray, wv::AbstractWeights)
    m = mean(A, wv)
    v = varm(A, wv, m, true)
    m, v
end

function mean_and_var(A::RealArray, wv::AbstractWeights, dim::Int)
    m = mean(A, wv, dim)
    v = varm(A, wv, m, dim, true)
    m, v
end

function mean_and_std(A::RealArray, wv::AbstractWeights)
    m = mean(A, wv)
    s = stdm(A, wv, m, true)
    m, s
end

function mean_and_std(A::RealArray, dim::Int)
    m = mean(A, dim)
    s = stdm(A, m, dim, true)
    m, s
end

function mean_and_std(A::RealArray, wv::AbstractWeights, dim::Int)
    m = mean(A, wv, dim)
    s = stdm(A, wv, m, dim, true)
    m, s
end

function Base.cov(x::DenseMatrix, wv::AbstractWeights)
    Base.depwarn(_correction_dep_msg("`cov`"), :cov)
    Base.covm(x, Base.mean(x, wv, 1), wv, 1, false)
end

function Base.cov(x::DenseMatrix, wv::AbstractWeights, vardim::Int)
    Base.depwarn(_correction_dep_msg("`cov`"), :cov)
    Base.covm(x, Base.mean(x, wv, vardim), wv, vardim, false)
end

function Base.covm(x::DenseMatrix, mean, wv::AbstractWeights)
    Base.depwarn(_correction_dep_msg("`covm`"), :covm)
    scale!(scattermatm(x, mean, wv, 1), varcorrection(wv, false))
end

function Base.covm(x::DenseMatrix, mean, wv::AbstractWeights, vardim::Int)
    Base.depwarn(_correction_dep_msg("`covm`"), :covm)
    scale!(scattermatm(x, mean, wv, vardim), varcorrection(wv, false))
end

function mean_and_cov(x::DenseMatrix, wv::AbstractWeights)
    m = mean(x, wv, 1)
    return m, Base.cov(x, wv, 1)
end

function mean_and_cov(x::DenseMatrix, wv::AbstractWeights, vardim::Int)
    m = mean(x, wv, vardim)
    return m, Base.cov(x, wv, vardim)
end
