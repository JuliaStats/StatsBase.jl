using ArrayViews
using StatsBase
using Base.Test

@test isa(weights([1, 2, 3]), WeightVec{Int})
@test isa(weights([1., 2., 3.]), WeightVec{Float64})

@test isa(weights([1 2 3; 4 5 6]), WeightVec{Int})

@test isempty(weights(Float64[]))

w  = [1., 2., 3.]
wv = weights(w)
@test eltype(wv) == Float64
@test length(wv) == 3
@test values(wv) === w
@test sum(wv) === 6.0
@test !isempty(wv)

b  = trues(3)
bv = weights(b)
@test eltype(bv) == Bool
@test length(bv) == 3
@test values(bv) == b
@test sum(bv)    == 3
@test !isempty(bv)

c  = [1, 4, 3, 5]; d = [4, 7, 18, 9]; e = [1, 2, -5, 3];
s  = sparse(c, d, e)
cv = weights(c)
@test eltype(cv) == Int64
@test length(cv) == 4
@test values(cv) == c
@test sum(cv)    === 13
@test !isempty(cv)

## wsum

x = [6., 8., 9.]
w = [2., 3., 4.]
p = [1. 2. ; 3. 4.]
q = [1., 2., 3., 4.]

@test wsum(Float64[], Float64[]) === 0.0
@test wsum(x, w) == 72.0
@test wsum(p, q) == 29.0

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


# Weighted quantile tests
data = (
    [7, 1, 2, 4, 10],
    [7, 1, 2, 4, 10],
    [7, 1, 2, 4, 10, 15],
    [1, 2, 4, 7, 10, 15],
    [0, 10, 20, 30],
    [1, 2, 3, 4, 5],
    [1, 2, 3, 4, 5],
    [30, 40, 50, 60, 35],
    [2, 0.6, 1.3, 0.3, 0.3, 1.7, 0.7, 1.7],
    [1, 2, 2],
    [3.7, 3.3, 3.5, 2.8],
    [100, 125, 123, 60, 45, 56, 66],
    [2, 2, 2, 2, 2, 2],
    [2.3],
    [-2, -3, 1, 2, -10],
    [1, 2, 3, 4, 5],
    [5, 4, 3, 2, 1],
    [-2, 2, -1, 3, 6],
    [-10, 1, 1, -10, -10],
)
wt = (
    weights([1, 1/3, 1/3, 1/3, 1]),
    weights([1, 1, 1, 1, 1]),
    weights([1, 1/3, 1/3, 1/3, 1, 1]),
    weights([1/3, 1/3, 1/3, 1, 1, 1]),
    weights([30, 191, 9, 0]),
    weights([10, 1, 1, 1, 9]),
    weights([10, 1, 1, 1, 900]),
    weights([1, 3, 5, 4, 2]),
    weights([2, 2, 5, 1, 2, 2, 1, 6]),
    weights([0.1, 0.1, 0.8]),
    weights([5, 5, 4, 1]),
    weights([30, 56, 144, 24, 55, 43, 67]),
    weights([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]),
    weights([12]),
    weights([7, 1, 1, 1, 6]),
    weights([1, 0, 0, 0, 2]),
    weights([1, 2, 3, 4, 5]),
    weights([0.1, 0.2, 0.3, 0.2, 0.1]),
    weights([1, 1, 1, 1, 1]),
)
quantile_answers = (
    [1.0,3.6000000000000005,6.181818181818182,8.2,10.0],
    [1.0,2.0,4.0,7.0,10.0],
    [1.0,4.75,8.0,10.833333333333334,15.0],
    [1.0,4.75,8.0,10.833333333333334,15.0],
    [0.0,6.1387900355871885,11.600000000000001,15.912500000000001,30.0],
    [1.0,1.5365853658536586,2.5999999999999996,4.405405405405405,5.0],
    [1.0,4.239377950569287,4.492918633712858,4.746459316856429,5.0],
    [30.0,38.75,45.714285714285715,52.85714285714286,60.0],
    [0.3,0.6903846153846154,1.484,1.7,2.0],
    [1.0,2.0,2.0,2.0,2.0],
    [2.8,3.3361111111111112,3.4611111111111112,3.581578947368421,3.7],
    [45.0,59.88593155893536,100.08846153846153,118.62115384615385,125.0],
    [2.0,2.0,2.0,2.0,2.0],
    [2.3,2.3,2.3,2.3,2.3],
    [-10.0,-5.52,-2.5882352941176467,-0.9411764705882351,2.0],
    [1.0,1.75,4.25,4.625,5.0],
    [1.0,1.625,2.3333333333333335,3.25,5.0],
    [-2.0,-0.5384615384615388,1.5384615384615383,2.6999999999999997,6.0],
    [-10.0,-10.0,-10.0,1.0,1.0]
)
p = [0.0, 0.25, 0.5, 0.75, 1.0]

srand(10)
for i = 1:length(data)
    @test_approx_eq quantile(data[i], wt[i], p) quantile_answers[i]
    for j = 1:10
        # order of p does not matter
        reorder = sortperm(rand(length(p)))
        @test_approx_eq quantile(data[i], wt[i], p[reorder]) quantile_answers[i][reorder]
    end
    for j = 1:10
        # order of w does not matter
        reorder = sortperm(rand(length(data[i])))
        @test_approx_eq quantile(data[i][reorder], weights(wt[i][reorder]), p) quantile_answers[i]
    end
end
# w = 1 corresponds to base quantile
for i = 1:length(data)
    @test_approx_eq quantile(data[i], weights(ones(Int64, length(data[i]))), p) quantile(data[i], p)
    for j = 1:10
        prandom = rand(4)
        @test_approx_eq quantile(data[i], weights(ones(Int64, length(data[i]))),  prandom) quantile(data[i], prandom)
    end
end

# other syntaxes
v = [7, 1, 2, 4, 10]
w = [1, 1/3, 1/3, 1/3, 1]
answer = 6.181818181818182
@test_approx_eq  quantile(data[1], weights(w), 0.5) answer
@test_approx_eq  wquantile(data[1], weights(w), [0.5]) [answer]
@test_approx_eq  wquantile(data[1], weights(w), 0.5) answer
@test_approx_eq  wquantile(data[1], w, [0.5]) [answer]
@test_approx_eq  wquantile(data[1], w, 0.5)  answer
