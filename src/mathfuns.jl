# Statistics-related math functions

# scalar functions

xlogx(x::FloatingPoint) = x > zero(x) ? x * log(x) : zero(x)
xlogx(x::Real) = xlogx(float(x))

xlogy{T<:FloatingPoint}(x::T, y::T) = x > zero(T) ? x * log(y) : zero(x)
xlogy{T<:Real}(x::T, y::T) = xlogy(float(x), float(y))
xlogy(x::Real, y::Real) = xlogy(promote(x, y)...)

logistic(x::FloatingPoint) = one(x) / (one(x) + exp(-x))
logistic(x::Real) = logistic(float(x))

logit(x::FloatingPoint) = log(x / (one(x) - x))
logit(x::Real) = logit(float(x))

softplus(x::FloatingPoint) = x <= 0 ? log1p(exp(x)) : x + log1p(exp(-x))
softplus(x::Real) = softplus(float(x))

invsoftplus(x::FloatingPoint) = log(exp(x) - one(x))
invsoftplus(x::Real) = invsoftplus(float(x))

@vectorize_1arg Real xlogx
@vectorize_2arg Real xlogy
@vectorize_1arg Real logistic
@vectorize_1arg Real logit
@vectorize_1arg Real softplus
@vectorize_1arg Real invsoftplus

## logsumexp

logsumexp{T<:FloatingPoint}(x::T, y::T) = x > y ? x + log1p(exp(y - x)) : y + log1p(exp(x - y))
logsumexp{T<:Real}(x::T, y::T) = logsumexp(float(x), float(y))
logsumexp(x::Real, y::Real) = logsumexp(promote(x, y)...)

function logsumexp{T<:Real}(x::AbstractArray{T})
    isempty(x) && return -Inf
    u = maximum(x)
    s = 0.
    for i = 1:length(x)
    	@inbounds s += exp(x[i] - u)
    end
    log(s) + u
end

## softmax

function softmax!{T<:Real}(r::AbstractArray, x::AbstractArray{T})
	n = length(x)
	length(r) == n || throw(DimensionMismatch("Inconsistent array lengths."))
	u = maximum(x)
	s = 0.
	for i = 1:n
		@inbounds s += (r[i] = exp(x[i] - u))
	end
	for i = 1:n
		r[i] /= s
	end
	r
end

softmax!{T<:FloatingPoint}(x::AbstractArray{T}) = softmax!(x, x)
softmax{T<:Real}(x::AbstractArray{T}) = softmax!(Array(Float64, size(x)), x)

