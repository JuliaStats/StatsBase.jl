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

@test_approx_eq mean([1:3;], weights([1.0, 1.0, 0.5])) 1.8
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

# Weighted median tests
data = (
    [7, 1, 2, 4, 10],
    [7, 1, 2, 4, 10],
    [7, 1, 2, 4, 10, 15],
    [1, 2, 4, 7, 10, 15],
    [0, 10, 20, 30],
    [1, 2, 3, 4, 5],
    [1, 2, 3, 4, 5],
    [30, 40, 50, 60, 35],
    [2, 0.6, 1.3, 0.3, 0.3, 1.7, 0.7, 1.7, 0.4],
    [3.7, 3.3, 3.5, 2.8],
    [100, 125, 123, 60, 45, 56, 66],
    [2, 2, 2, 2, 2, 2],
    [2.3],
    [-2, -3, 1, 2, -10],
    [1, 2, 3, 4, 5],
    [5, 4, 3, 2, 1],
    [-2, 2, -1, 3, 6],
    [-10, 1, 1, -10, -10],
    [2, 4],
    [2, 2, 4, 4],    
    [2, 2, 2, 4]   
)
wt = (
    [1, 1/3, 1/3, 1/3, 1],
    [1, 1, 1, 1, 1],
    [1, 1/3, 1/3, 1/3, 1, 1],
    [1/3, 1/3, 1/3, 1, 1, 1],
    [30, 191, 9, 0],
    [10, 1, 1, 1, 9],
    [10, 1, 1, 1, 900],
    [1, 3, 5, 4, 2],
    [2, 2, 0, 1, 2, 2, 1, 6, 0],
    [5, 5, 4, 1],
    [30, 56, 144, 24, 55, 43, 67],
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    [12],
    [7, 1, 1, 1, 6],
    [1, 0, 0, 0, 2],
    [1, 2, -3, 4, -5],
    [0.1, 0.2, 0.3, -0.2, 0.1],
    [-1, -1, -1, -1, 1],
    [1, 1],
    [1, 1, 1, 1],  
    [1, 1, 1, 1]
)
median_answers = (7.0,   4.0,   8.5,
                  8.5,  10.0,   2.5,
                  5.0,  50.0,   1.7,
                  3.5, 100.0,   2.0,
                  2.3,  -2.0,   5.0,
                  2.0,  -1.0, -10.0,
                  3.0,   3.0,   2.0)
num_tests = length(data)
for i = 1:num_tests
    @test wmedian(data[i], wt[i]) == median_answers[i]
    @test wmedian(data[i], weights(wt[i])) == median_answers[i]
    @test median(data[i], weights(wt[i])) == median_answers[i]
    for j = 1:100
        # Make sure the weighted median does not change if the data
        # and weights are reordered.
        reorder = sortperm(rand(length(data[i])))
        @test median(data[i][reorder], weights(wt[i][reorder])) == median_answers[i]
    end
end
data = [4, 3, 2, 1]
wt = [0, 0, 0, 0]
@test_throws MethodError wmedian(data[1])
@test_throws ErrorException median(data, weights(wt))
@test_throws ErrorException wmedian(data, wt)
@test_throws ErrorException median((Float64)[], weights((Float64)[]))
wt = [1, 2, 3, 4, 5]
@test_throws ErrorException median(data, weights(wt))
@test_throws MethodError median([4 3 2 1 0], weights(wt))
@test_throws MethodError median([[1 2];[4 5];[7 8];[10 11];[13 14]], weights(wt))
data = [1, 3, 2, NaN, 2]
@test isnan(median(data, weights(wt)))
wt = [1, 2, NaN, 4, 5]
@test_throws ErrorException median(data, weights(wt))
data = [1, 3, 2, 1, 2]
@test_throws ErrorException median(data, weights(wt))
wt = [-1, -1, -1, -1, -1]
@test_throws ErrorException median(data, weights(wt))
wt = [-1, -1, -1, 0, 0]
@test_throws ErrorException median(data, weights(wt))
