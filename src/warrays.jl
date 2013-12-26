# Arrays associated with weights

immutable WeightedVector{VVec<:AbstractVector,WVec<:AbstractVector,W<:Real}
	values::VVec
	weights::WVec
	wsum::W
end

function WeightedVector{VV<:AbstractVector,WV<:AbstractVector,W<:Real}(v::VV, w::WV, wsum::W)
	length(v) == length(w) || error("Inconsistent argument lengths.")
	WeightedVector{VV,WV,W}(v, w, wsum)
end

WeightedVector{VV<:AbstractVector,WV<:AbstractVector}(v::VV, w::WV) = WeightedVector(v, w, sum(w))

values(wv::WeightedVector) = wv.values
weights(wv::WeightedVector) = wv.weights
weightsum(wv::WeightedVector) = wv.wsum

