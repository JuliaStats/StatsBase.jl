## WeightVec

using ArrayViews
using StatsBase
using Base.Test

@test isa(weights([1, 2, 3]), WeightVec{Int})
@test isa(weights([1., 2., 3.]), WeightVec{Float64})

@test isempty(weights(Float64[]))

w = [1., 2., 3.]
wv = weights(w)
@test eltype(wv) == Float64
@test length(wv) == 3
@test values(wv) === w
@test sum(wv) == 6.0
@test !isempty(wv) 

## wsum 

x = [6., 8., 9.]
w = [2., 3., 4.]

@test wsum(Float64[], Float64[]) === 0.0
@test wsum(x, w) == 72.0

## wsum along dimensions

@test wsum(x, w, 1) == [72.0]

x = rand(6, 8)
w1 = rand(6)
w2 = rand(8)

@test size(wsum(x, w1, 1)) == (1, 8)
@test size(wsum(x, w2, 2)) == (6, 1)

@test_approx_eq wsum(x, w1, 1) sum(x .* w1, 1)
@test_approx_eq wsum(x, w2, 2) sum(x .* w2', 2)

x = rand(6, 5, 4)
w1 = rand(6)
w2 = rand(5)
w3 = rand(4)

@test size(wsum(x, w1, 1)) == (1, 5, 4)
@test size(wsum(x, w2, 2)) == (6, 1, 4)
@test size(wsum(x, w3, 3)) == (6, 5, 1)

@test_approx_eq wsum(x, w1, 1) sum(x .* w1, 1)
@test_approx_eq wsum(x, w2, 2) sum(x .* w2', 2)
@test_approx_eq wsum(x, w3, 3) sum(x .* reshape(w3, 1, 1, 4), 3)

v = view(x, 2:4, :, :)

@test_approx_eq wsum(v, w1[1:3], 1) sum(v .* w1[1:3], 1)
@test_approx_eq wsum(v, w2, 2) sum(v .* w2', 2)
@test_approx_eq wsum(v, w3, 3) sum(v .* reshape(w3, 1, 1, 4), 3)

## wsum for Arrays with non-BlasReal elements

x = rand(1:100, 6, 8)
w1 = rand(6)
w2 = rand(8)

@test_approx_eq wsum(x, w1, 1) sum(x .* w1, 1)
@test_approx_eq wsum(x, w2, 2) sum(x .* w2', 2)

## wsum!

x = rand(6)
w = rand(6)

r = ones(1)
@test wsum!(r, x, w, 1; init=true) === r
@test_approx_eq r [dot(x, w)]

r = ones(1)
@test wsum!(r, x, w, 1; init=false) === r
@test_approx_eq r [dot(x, w) + 1.0]

x = rand(6, 8)
w1 = rand(6)
w2 = rand(8)

r = ones(1, 8)
@test wsum!(r, x, w1, 1; init=true) === r
@test_approx_eq r sum(x .* w1, 1)

r = ones(1, 8)
@test wsum!(r, x, w1, 1; init=false) === r
@test_approx_eq r sum(x .* w1, 1) .+ 1.0

r = ones(6)
@test wsum!(r, x, w2, 2; init=true) === r
@test_approx_eq r sum(x .* w2', 2)

r = ones(6)
@test wsum!(r, x, w2, 2; init=false) === r
@test_approx_eq r sum(x .* w2', 2) .+ 1.0

x = rand(8, 6, 5)
w1 = rand(8)
w2 = rand(6)
w3 = rand(5)

r = ones(1, 6, 5)
@test wsum!(r, x, w1, 1; init=true) === r
@test_approx_eq r sum(x .* w1, 1)

r = ones(1, 6, 5)
@test wsum!(r, x, w1, 1; init=false) === r
@test_approx_eq r sum(x .* w1, 1) .+ 1.0

r = ones(8, 1, 5)
@test wsum!(r, x, w2, 2; init=true) === r
@test_approx_eq r sum(x .* w2', 2)

r = ones(8, 1, 5)
@test wsum!(r, x, w2, 2; init=false) === r
@test_approx_eq r sum(x .* w2', 2) .+ 1.0

r = ones(8, 6)
@test wsum!(r, x, w3, 3; init=true) === r
@test_approx_eq r sum(x .* reshape(w3, (1, 1, 5)), 3)

r = ones(8, 6)
@test wsum!(r, x, w3, 3; init=false) === r
@test_approx_eq r sum(x .* reshape(w3, (1, 1, 5)), 3) .+ 1.0


## the sum and mean syntax

@test_approx_eq sum([1.0, 2.0, 3.0], weights([1.0, 0.5, 0.5])) 3.5
@test_approx_eq sum(1:3, weights([1.0, 1.0, 0.5])) 4.5

@test_approx_eq mean([1:3], weights([1.0, 1.0, 0.5])) 1.8
@test_approx_eq mean(1:3, weights([1.0, 1.0, 0.5])) 1.8

a = reshape(1.0:27.0, 3, 3, 3)
for wt in ([1.0, 1.0, 1.0], [1.0, 0.2, 0.0], [0.2, 0.0, 1.0])
	@test_approx_eq sum(a, weights(wt), 1) sum(a.*reshape(wt, length(wt), 1, 1), 1)
	@test_approx_eq sum(a, weights(wt), 2) sum(a.*reshape(wt, 1, length(wt), 1), 2)
	@test_approx_eq sum(a, weights(wt), 3) sum(a.*reshape(wt, 1, 1, length(wt)), 3)
	@test_approx_eq mean(a, weights(wt), 1) sum(a.*reshape(wt, length(wt), 1, 1), 1)/sum(wt)
	@test_approx_eq mean(a, weights(wt), 2) sum(a.*reshape(wt, 1, length(wt), 1), 2)/sum(wt)
	@test_approx_eq mean(a, weights(wt), 3) sum(a.*reshape(wt, 1, 1, length(wt)), 3)/sum(wt)
	@test_throws ErrorException mean(a, weights(wt), 4)
end

