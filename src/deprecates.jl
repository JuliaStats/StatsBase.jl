import Base.@deprecate
import Base.depwarn

import Base.varm, Base.stdm
@deprecate varm(v::RealArray, m::Real, wv::WeightVec) varm(v, wv, m)
@deprecate varm(A::RealArray, M::RealArray, wv::WeightVec, dim::Int) varm(v, wv, m, dim)
@deprecate stdm(v::RealArray, m::Real, wv::WeightVec) stdm(v, wv, m)
@deprecate stdm(v::RealArray, m::RealArray, wv::WeightVec, dim::Int) stdm(v, wv, m, dim)

@deprecate _moment2(v::RealArray, m::Real, wv::WeightVec) _moment2(v, wv, m)
@deprecate _moment3(v::RealArray, m::Real, wv::WeightVec) _moment3(v, wv, m)
@deprecate _moment4(v::RealArray, m::Real, wv::WeightVec) _moment4(v, wv, m)
@deprecate _momentk(v::RealArray, k::Int, m::Real, wv::WeightVec) _momentk(v, k, wv, m)
@deprecate moment(v::RealArray, k::Int, m::Real, wv::WeightVec) moment(v, k, wv, m)

@deprecate AIC(obj::StatisticalModel) aic(obj)
@deprecate AICc(obj::StatisticalModel) aicc(obj)
@deprecate BIC(obj::StatisticalModel) bic(obj)

@deprecate R2(obj::StatisticalModel, variant::Symbol) r2(obj, variant)
@deprecate R²(obj::StatisticalModel, variant::Symbol) r²(obj, variant)
@deprecate adjR2(obj::StatisticalModel, variant::Symbol) adjr2(obj, variant)
@deprecate adjR²(obj::StatisticalModel, variant::Symbol) adjr²(obj, variant)

