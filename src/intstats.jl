# Statistics related to a collection of integers
#
# Note: Original versions of these functions were provided by MLBase.jl
#       After they are migrated here, there has been some clean-up of
#       both the interface & implementation.
#  

#### functions for counting a single list of integers (1D)

function addcounts!{T<:Real}(r::AbstractArray, x::AbstractArray{T})
	# add counts of integers from x to r

	k = convert(T, length(r))
	for i in 1 : length(x)
		@inbounds xi = x[i]
		if 1 <= xi <= k
			@inbounds r[xi] += 1  
		end
	end
	return r
end

function addwcounts!{T<:Real,W<:Real}(r::AbstractArray, x::AbstractArray{T}, w::AbstractArray{W})
	k = convert(T, length(r))
	for i in 1 : length(x)
		@inbounds xi = x[i]
		if 1 <= xi <= k
			@inbounds r[xi] += w[i]
		end
	end
	return r
end

# This function returns a vector r, s.t. r[i] == sum(x .== i)
counts{T<:Real}(x::AbstractArray{T}, k::Integer) = addcounts!(zeros(Int, k), x)
counts{T<:Real}(x::AbstractArray{T}) = counts(x, int(maximum(x)))

wcounts{T<:Real,W<:Real}(x::AbstractArray{T}, w::AbstractArray{W}, k::Integer) = addwcounts!(zeros(W, k), x, w)
wcounts{T<:Real,W<:Real}(x::AbstractArray{T}, w::AbstractArray{W}) = wcounts(x, w, int(maximum(x)))


### functions for counting over two lists of indexes
### (this is useful in many applications, e.g. computing confusion matrix)

function addcounts!{T<:Real}(r::AbstractMatrix, x1::AbstractArray{T}, x2::AbstractArray{T})
	n = length(x1)
	length(x2) == n || throw(ArgumentError("Inconsistent array lengths."))

	k1 = size(r, 1)
	k2 = size(r, 2)

	for i in 1 : n
		@inbounds x1i = x1[i]
		@inbounds x2i = x2[i]
		if 1 <= x1i <= k1 && 1 <= x2i <= k2
			r[x1i, x2i] += 1
		end
	end
	r
end

function addwcounts!{T<:Real,W<:Real}(r::AbstractMatrix, x1::AbstractArray{T}, x2::AbstractArray{T}, w::AbstractArray{W})
	n = length(x1)
	length(x2) == n || throw(ArgumentError("Inconsistent array lengths."))

	k1 = size(r, 1)
	k2 = size(r, 2)

	for i in 1 : n
		@inbounds x1i = x1[i]
		@inbounds x2i = x2[i]
		if 1 <= x1i <= k1 && 1 <= x2i <= k2
			r[x1i, x2i] += w[i]
		end
	end
	r
end

counts{T<:Real}(x1::AbstractArray{T}, x2::AbstractArray{T}, rsiz::(Int, Int)) = addcounts!(zeros(Int, rsiz), x1, x2)
counts{T<:Real}(x1::AbstractArray{T}, x2::AbstractArray{T}, rd::Int) = addcounts!(zeros(Int, rd, rd), x1, x2)
counts{T<:Real}(x1::AbstractArray{T}, x2::AbstractArray{T}) = counts(x1, x2, (int(maximum(x1)), int(maximum(x2))))

function wcounts{T<:Real,W<:Real}(x1::AbstractArray{T}, x2::AbstractArray{T}, w::AbstractArray{W}, rsiz::(Int, Int))
	addwcounts!(zeros(W, rsiz), x1, x2, w)
end

function wcounts{T<:Real,W<:Real}(x1::AbstractArray{T}, x2::AbstractArray{T}, w::AbstractArray{W}, rd::Int)
	addwcounts!(zeros(W, rd, rd), x1, x2, w)
end

function wcounts{T<:Real,W<:Real}(x1::AbstractArray{T}, x2::AbstractArray{T}, w::AbstractArray{W})
	wcounts(x1, x2, w, (int(maximum(x1)), int(maximum(x2))))
end

