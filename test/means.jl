using StatsBase
using Base.Test

@test_approx_eq mean([1, 2, 3]) 2.0
@test_approx_eq mean(1:3) 2.0

@test_approx_eq geomean([1, 2, 3]) (6.0)^(1/3)
@test_approx_eq geomean(1:3) (6.0)^(1/3)
@test_approx_eq geomean([2, 8]) 4.0
@test_approx_eq geomean([4, 1, 1/32]) 0.5
@test geomean([1, 0, 2]) == 0.0

@test_approx_eq harmmean([1, 2, 3]) 3 / (1 + 1/2 + 1/3)
@test_approx_eq harmmean(1:3) 3 / (1 + 1/2 + 1/3)
@test_approx_eq harmmean([1, 2, 4]) 12 / 7

@test_approx_eq trimmean([-100, 2, 3, 7, 200], 0.0) 22.4
@test_approx_eq trimmean([-100, 2, 3, 7, 200], 0.4) 4.0
@test_approx_eq trimmean([-100, 2, 3, 7, 200], 0.8) 3.0

@test_approx_eq sum([1.0, 2.0, 3.0], weights([1/3, 1/3, 1/3])) 2.0
@test_approx_eq sum([1.0, 2.0, 3.0], weights([1.0, 0.0, 0.0])) 1.0
@test_approx_eq sum([1.0, 2.0, 3.0], weights([0.0, 1.0, 0.0])) 2.0
@test_approx_eq sum([1.0, 2.0, 3.0], weights([0.0, 0.0, 1.0])) 3.0
@test_approx_eq sum([1.0, 2.0, 3.0], weights([0.5, 0.0, 0.5])) 2.0
@test_approx_eq sum([1.0, 2.0, 3.0], weights([0.5, 0.5, 0.0])) 1.5
@test_approx_eq sum([1.0, 2.0, 3.0], weights([0.0, 0.5, 0.5])) 2.5

@test_approx_eq sum(1:3, weights([1/3, 1/3, 1/3])) 2.0
@test_approx_eq sum(1:3, weights([1.0, 0.0, 0.0])) 1.0
@test_approx_eq sum(1:3, weights([0.0, 1.0, 0.0])) 2.0
@test_approx_eq sum(1:3, weights([0.0, 0.0, 1.0])) 3.0
@test_approx_eq sum(1:3, weights([0.5, 0.0, 0.5])) 2.0
@test_approx_eq sum(1:3, weights([0.5, 0.5, 0.0])) 1.5
@test_approx_eq sum(1:3, weights([0.0, 0.5, 0.5])) 2.5
@test_approx_eq sum(1:3, weights([1.0, 1.0, 0.5])) 4.5
@test_approx_eq mean(1:3, weights([1.0, 1.0, 0.5])) 1.8

a = [1. 2. 3.; 4. 5. 6.]

@test size(mean(a, weights(ones(2)), 1)) == (1, 3)
@test_approx_eq sum(a, weights([1.0, 1.0]), 1) [5.0, 7.0, 9.0]
@test_approx_eq mean(a, weights([1.0, 1.0]), 1) [2.5, 3.5, 4.5]
@test_approx_eq sum(a, weights([1.0, 0.0]), 1) [1.0, 2.0, 3.0]
@test_approx_eq sum(a, weights([0.0, 1.0]), 1) [4.0, 5.0, 6.0]

@test size(mean(a, weights(ones(3)), 2)) == (2, 1)
@test_approx_eq sum(a, weights([1.0, 1.0, 1.0]), 2) [6.0 15.0]
@test_approx_eq mean(a, weights([1.0, 1.0, 1.0]), 2) [2.0 5.0]
@test_approx_eq sum(a, weights([1.0, 0.0, 0.0]), 2) [1.0 4.0]
@test_approx_eq sum(a, weights([0.0, 0.0, 1.0]), 2) [3.0 6.0]

@test_throws ErrorException mean(a, weights(ones(3)), 3)
@test_throws DimensionMismatch mean(a, weights(ones(2)), 2)
@test_throws DimensionMismatch mean!(ones(1, 1), a, weights(ones(3)), 2)

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
