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

@test_approx_eq mean([1.0, 2.0, 3.0], weights([1/3, 1/3, 1/3])) 2.0
@test_approx_eq mean([1.0, 2.0, 3.0], weights([1.0, 0.0, 0.0])) 1.0
@test_approx_eq mean([1.0, 2.0, 3.0], weights([0.0, 1.0, 0.0])) 2.0
@test_approx_eq mean([1.0, 2.0, 3.0], weights([0.0, 0.0, 1.0])) 3.0
@test_approx_eq mean([1.0, 2.0, 3.0], weights([0.5, 0.0, 0.5])) 2.0
@test_approx_eq mean([1.0, 2.0, 3.0], weights([0.5, 0.5, 0.0])) 1.5
@test_approx_eq mean([1.0, 2.0, 3.0], weights([0.0, 0.5, 0.5])) 2.5

@test_approx_eq mean(1:3, weights([1/3, 1/3, 1/3])) 2.0
@test_approx_eq mean(1:3, weights([1.0, 0.0, 0.0])) 1.0
@test_approx_eq mean(1:3, weights([0.0, 1.0, 0.0])) 2.0
@test_approx_eq mean(1:3, weights([0.0, 0.0, 1.0])) 3.0
@test_approx_eq mean(1:3, weights([0.5, 0.0, 0.5])) 2.0
@test_approx_eq mean(1:3, weights([0.5, 0.5, 0.0])) 1.5
@test_approx_eq mean(1:3, weights([0.0, 0.5, 0.5])) 2.5
@test_approx_eq mean(1:3, weights([1.0, 1.0, 0.5])) 1.8

