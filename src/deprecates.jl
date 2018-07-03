using Base: @deprecate, @deprecate_binding, depwarn

@deprecate varm(v::RealArray, m::Real, wv::AbstractWeights) varm(v, wv, m)
@deprecate varm(A::RealArray, M::RealArray, wv::AbstractWeights, dim::Int) varm(v, wv, m, dim)
@deprecate stdm(v::RealArray, m::Real, wv::AbstractWeights) stdm(v, wv, m)
@deprecate stdm(v::RealArray, m::RealArray, wv::AbstractWeights, dim::Int) stdm(v, wv, m, dim)

@deprecate trimmean(x::RealArray, p::Real) mean(trim(x, p/2))

@deprecate _moment2(v::RealArray, m::Real, wv::AbstractWeights) _moment2(v, wv, m)
@deprecate _moment3(v::RealArray, m::Real, wv::AbstractWeights) _moment3(v, wv, m)
@deprecate _moment4(v::RealArray, m::Real, wv::AbstractWeights) _moment4(v, wv, m)
@deprecate _momentk(v::RealArray, k::Int, m::Real, wv::AbstractWeights) _momentk(v, k, wv, m)
@deprecate moment(v::RealArray, k::Int, m::Real, wv::AbstractWeights) moment(v, k, wv, m)

@deprecate AIC(obj::StatisticalModel) aic(obj)
@deprecate AICc(obj::StatisticalModel) aicc(obj)
@deprecate BIC(obj::StatisticalModel) bic(obj)

if !isdefined(Base, :stderr)
    @deprecate stderr(obj::StatisticalModel) stderror(obj)
else
    function (io::typeof(stderr))(obj::StatisticalModel)
        Base.depwarn("stderr(obj::StatisticalModel) is deprecated, use stderror(obj) instead", :stderr)
        io === stderr ? stderror(obj) : throw(MethodErrror(io, (obj,)))
    end
end

@deprecate R2(obj::StatisticalModel, variant::Symbol) r2(obj, variant)
@deprecate R²(obj::StatisticalModel, variant::Symbol) r²(obj, variant)
@deprecate adjR2(obj::StatisticalModel, variant::Symbol) adjr2(obj, variant)
@deprecate adjR²(obj::StatisticalModel, variant::Symbol) adjr²(obj, variant)
@deprecate model_response(obj::StatisticalModel) response(obj)

@deprecate norepeats(a::AbstractArray) allunique(a)

function findat!(r::IntegerArray, a::AbstractArray{T}, b::AbstractArray{T}) where T
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
findat(a::AbstractArray, b::AbstractArray) = findat!(Array{Int}(undef, size(b)), a, b)

@deprecate df(obj::StatisticalModel) dof(obj)
@deprecate df_residual(obj::StatisticalModel) dof_residual(obj)

@deprecate_binding WeightVec Weights

struct RandIntSampler  # for generating Int samples in [0, K-1]
    a::Int
    Ku::UInt
    U::UInt

    function RandIntSampler(K::Int)
        Base.depwarn("RandIntSampler is deprecated, use Base.Random.RangeGenerator instead",
                     :RandIntSampler)
        Ku = UInt(K)
        new(1, Ku, div(typemax(UInt), Ku) * Ku)
    end
    function RandIntSampler(a::Int, b::Int)
        Base.depwarn("RandIntSampler is deprecated, use Base.Random.RangeGenerator instead",
                     :RandIntSampler)
        Ku = UInt(b-a+1)
        new(a, Ku, div(typemax(UInt), Ku) * Ku)
    end
end

function rand(rng::AbstractRNG, s::RandIntSampler)
    x = rand(rng, UInt)
    while x >= s.U
        x = rand(rng, UInt)
    end
    s.a + Int(rem(x, s.Ku))
end
rand(s::RandIntSampler) = rand(Random.GLOBAL_RNG, s)

@deprecate randi(rng::AbstractRNG, K::Int) rand(rng, 1:K)
@deprecate randi(K::Int) rand(1:K)
@deprecate randi(rng::AbstractRNG, a::Int, b::Int) rand(rng, a:b)
@deprecate randi(a::Int, b::Int) rand(a:b)

@deprecate(mad!(v::AbstractArray{<:Real}, center;
                constant::Real = BigFloat("1.482602218505601860547076529360423431326703202590312896536266275245674447622701")),
           mad!(v, center=center, constant=constant))
