# Counts of discrete values

#################################################
#
#  counts on integer ranges
#
#  These methods are the most efficient if one
#  knows that all values fall in a small range
#  of integers.
#  
#################################################

#### functions for counting a single list of integers (1D)

function addcounts!{T<:Integer}(r::AbstractArray, x::AbstractArray{T}, rgn::Range1{T})
	# add counts of integers from x to r

	k = length(rgn)
	length(r) == k || error("Inconsistent argument lengths.")

	m0 = rgn[1]
	m1 = rgn[end]
	b = m0 - one(T)
	
	for i in 1 : length(x)
		@inbounds xi = x[i]
		if m0 <= xi <= m1
			@inbounds r[xi - b] += 1  
		end
	end
	return r
end

function addwcounts!{T<:Integer,W<:Real}(r::AbstractArray, x::AbstractArray{T}, w::AbstractArray{W}, rgn::Range1{T})
	k = length(rgn)
	length(r) == k || error("Inconsistent argument lengths.")

	m0 = rgn[1]
	m1 = rgn[end]
	b = m0 - one(T)
	
	for i in 1 : length(x)
		@inbounds xi = x[i]
		if m0 <= xi <= m1
			@inbounds r[xi - b] += w[i]  
		end
	end
	return r
end

counts{T<:Integer}(x::AbstractArray{T}, rgn::Range1{T}) = addcounts!(zeros(Int, length(rgn)), x, rgn)
wcounts{T<:Integer,W<:Real}(x::AbstractArray{T}, w::AbstractArray{W}, rgn::Range1{T}) = addwcounts!(zeros(W, length(rgn)), x, w, rgn)

