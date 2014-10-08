# Counts of discrete values

#################################################
#
#  counts on given levels
#  
#################################################

typealias IntRange1{T<:Integer} Range1{T}

#### functions for counting a single list of integers (1D)

function addcounts!(r::AbstractArray, x::IntegerArray, levels::IntRange1)
	# add counts of integers from x to r

	k = length(levels)
	length(r) == k || raise_dimerror()

	m0 = levels[1]
	m1 = levels[end]
	b = m0 - 1
	
	@inbounds for i in 1 : length(x)
		xi = x[i]
		if m0 <= xi <= m1
			r[xi - b] += 1  
		end
	end
	return r
end

function addcounts!(r::AbstractArray, x::IntegerArray, levels::IntRange1, wv::WeightVec)
	k = length(levels)
	length(r) == k || raise_dimerror()

	m0 = levels[1]
	m1 = levels[end]
	b = m0 - 1
	w = values(wv)
	
	@inbounds for i in 1 : length(x)
		xi = x[i]
		if m0 <= xi <= m1
			r[xi - b] += w[i]  
		end
	end
	return r
end

counts(x::IntegerArray, levels::IntRange1) = addcounts!(zeros(Int, length(levels)), x, levels)
counts(x::IntegerArray, levels::IntRange1, wv::WeightVec) = addcounts!(zeros(eltype(wv), length(levels)), x, levels, wv)
counts(x::IntegerArray, k::Integer) = counts(x, 1:k)
counts(x::IntegerArray, k::Integer, wv::WeightVec) = counts(x, 1:k, wv)

proportions(x::IntegerArray, levels::IntRange1) = counts(x, levels) .* inv(length(x))
proportions(x::IntegerArray, levels::IntRange1, wv::WeightVec) = counts(x, levels, wv) .* inv(sum(wv))
proportions(x::IntegerArray, k::Integer) = proportions(x, 1:k)
proportions(x::IntegerArray, k::Integer, wv::WeightVec) = proportions(x, 1:k, wv)

#### functions for counting a single list of integers (2D)

function addcounts!(r::AbstractArray, x::IntegerArray, y::IntegerArray, levels::(IntRange1, IntRange1))
	# add counts of integers from x to r

	n = length(x)
	length(y) == n || raise_dimerror()

	xlevels, ylevels = levels

	kx = length(xlevels)
	ky = length(ylevels)
	size(r) == (kx, ky) || raise_dimerror()

	mx0 = xlevels[1]
	mx1 = xlevels[end]
	my0 = ylevels[1]
	my1 = ylevels[end]

	bx = mx0 - 1
	by = my0 - 1
	
	for i = 1:n
		xi = x[i]
		yi = y[i]
		if (mx0 <= xi <= mx1) && (my0 <= yi <= my1)
			r[xi - bx, yi - by] += 1  
		end
	end
	return r
end

function addcounts!(r::AbstractArray, x::IntegerArray, y::IntegerArray, levels::(IntRange1, IntRange1), wv::WeightVec)
	# add counts of integers from x to r

	n = length(x)
	length(y) == length(wv) == n || raise_dimerror()

	xlevels, ylevels = levels

	kx = length(xlevels)
	ky = length(ylevels)
	size(r) == (kx, ky) || raise_dimerror()

	mx0 = xlevels[1]
	mx1 = xlevels[end]
	my0 = ylevels[1]
	my1 = ylevels[end]

	bx = mx0 - 1
	by = my0 - 1
	w = values(wv)
	
	for i = 1:n
		xi = x[i]
		yi = y[i]
		if (mx0 <= xi <= mx1) && (my0 <= yi <= my1)
			r[xi - bx, yi - by] += w[i] 
		end
	end
	return r
end

# facet functions

function counts(x::IntegerArray, y::IntegerArray, levels::(IntRange1, IntRange1))
	addcounts!(zeros(Int, length(levels[1]), length(levels[2])), x, y, levels)
end

function counts(x::IntegerArray, y::IntegerArray, levels::(IntRange1, IntRange1), wv::WeightVec)
	addcounts!(zeros(eltype(wv), length(levels[1]), length(levels[2])), x, y, levels, wv)
end

counts(x::IntegerArray, y::IntegerArray, levels::IntRange1) = counts(x, y, (levels, levels))
counts(x::IntegerArray, y::IntegerArray, levels::IntRange1, wv::WeightVec) = counts(x, y, (levels, levels), wv)

counts(x::IntegerArray, y::IntegerArray, ks::(Integer, Integer)) = counts(x, y, (1:ks[1], 1:ks[2]))
counts(x::IntegerArray, y::IntegerArray, ks::(Integer, Integer), wv::WeightVec) = counts(x, y, (1:ks[1], 1:ks[2]), wv)
counts(x::IntegerArray, y::IntegerArray, k::Integer) = counts(x, y, (1:k, 1:k))
counts(x::IntegerArray, y::IntegerArray, k::Integer, wv::WeightVec) = counts(x, y, (1:k, 1:k), wv)

proportions(x::IntegerArray, y::IntegerArray, levels::(IntRange1, IntRange1)) = counts(x, y, levels) .* inv(length(x))
proportions(x::IntegerArray, y::IntegerArray, levels::(IntRange1, IntRange1), wv::WeightVec) = counts(x, y, levels, wv) .* inv(sum(wv))

proportions(x::IntegerArray, y::IntegerArray, ks::(Integer, Integer)) = proportions(x, y, (1:ks[1], 1:ks[2]))
proportions(x::IntegerArray, y::IntegerArray, ks::(Integer, Integer), wv::WeightVec) = proportions(x, y, (1:ks[1], 1:ks[2]), wv)
proportions(x::IntegerArray, y::IntegerArray, k::Integer) = proportions(x, y, (1:k, 1:k))
proportions(x::IntegerArray, y::IntegerArray, k::Integer, wv::WeightVec) = proportions(x, y, (1:k, 1:k), wv)


#################################################
#
#  countmap on unknown levels
#
#  These methods are based on dictionaries, and 
#  can be used on any kind of hashable values.
#  
#################################################

## auxiliary functions

function _normalize_countmap{T}(cm::Dict{T}, s::Real)
	r = Dict{T,Float64}[]
	for (k, c) in cm
		r[k] = c / s 
	end
	return r
end

## 1D

function addcounts!{T}(cm::Dict{T}, x::AbstractArray{T})
	for v in x
		cm[v] = get(cm, v, 0) + 1
	end
	return cm
end

function addcounts!{T,W}(cm::Dict{T}, x::AbstractArray{T}, wv::WeightVec{W})
	n = length(x)
	length(wv) == n || raise_dimerror()
	w = values(wv)
	z = zero(W)

	for i = 1 : n
		@inbounds xi = x[i]
		@inbounds wi = w[i]
		cm[xi] = get(cm, xi, z) + wi
	end
	return cm
end

countmap{T}(x::AbstractArray{T}) = addcounts!(Dict{T,Int}[], x)
countmap{T,W}(x::AbstractArray{T}, wv::WeightVec{W}) = addcounts!(Dict{T,W}[], x, wv)

proportionmap(x::AbstractArray) = _normalize_countmap(countmap(x), length(x))
proportionmap(x::AbstractArray, wv::WeightVec) = _normalize_countmap(countmap(x, wv), sum(wv))

