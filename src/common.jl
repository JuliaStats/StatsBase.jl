# common utilities


## convenient type alias
#
#  These types signficantly reduces the need of using
#  type parameters in functions (which are often just
#  for the purpose of restricting the arrays to real)
#

typealias RealArray{T<:Real,N} AbstractArray{T,N}
typealias RealVector{T<:Real} AbstractArray{T,1}
typealias RealMatrix{T<:Real} AbstractArray{T,2}

# weight vector: 
#
#    a wrapper that indicator a vector representing a sequence of weights
#
immutable WeightVec{Vec<:RealVector,W}
	values::Vec
	sum::W
end

WeightVec{Vec<:RealVector,W<:Real}(vs::Vec,wsum::W) = WeightVec{Vec,W}(vs, wsum)
WeightVec(vs::RealVector) = WeightVec(vs, sum(vs))

weights(vs::RealVector) = WeightVec(vs)
weights(vs::RealArray) = WeightVec(vec(vs))

length(wv::WeightVec) = length(wv.values)
values(wv::WeightVec) = wv.values
sum(wv::WeightVec) = wv.sum

