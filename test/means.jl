using Stats
using Base.Test

@test_approx_eq mean([1, 2, 3]) 2.0
@test_approx_eq mean(1:3) 2.0

@test_approx_eq gmean([1, 2, 3]) (6.0)^(1/3)
@test_approx_eq gmean(1:3) (6.0)^(1/3)
@test_approx_eq gmean([2, 8]) 4.0
@test_approx_eq gmean([4, 1, 1/32]) 0.5
@test gmean([1, 0, 2]) == 0.0

@test_approx_eq hmean([1, 2, 3]) 3 / (1 + 1/2 + 1/3)
@test_approx_eq hmean(1:3) 3 / (1 + 1/2 + 1/3)
@test_approx_eq hmean([1, 2, 4]) 12 / 7

@test_approx_eq wmean([1.0, 2.0, 3.0], [1/3, 1/3, 1/3]) 2.0
@test_approx_eq wmean([1.0, 2.0, 3.0], [1.0, 0.0, 0.0]) 1.0
@test_approx_eq wmean([1.0, 2.0, 3.0], [0.0, 1.0, 0.0]) 2.0
@test_approx_eq wmean([1.0, 2.0, 3.0], [0.0, 0.0, 1.0]) 3.0
@test_approx_eq wmean([1.0, 2.0, 3.0], [0.5, 0.0, 0.5]) 2.0
@test_approx_eq wmean([1.0, 2.0, 3.0], [0.5, 0.5, 0.0]) 1.5
@test_approx_eq wmean([1.0, 2.0, 3.0], [0.0, 0.5, 0.5]) 2.5

@test_approx_eq wmean(1:3, [1/3, 1/3, 1/3]) 2.0
@test_approx_eq wmean(1:3, [1.0, 0.0, 0.0]) 1.0
@test_approx_eq wmean(1:3, [0.0, 1.0, 0.0]) 2.0
@test_approx_eq wmean(1:3, [0.0, 0.0, 1.0]) 3.0
@test_approx_eq wmean(1:3, [0.5, 0.0, 0.5]) 2.0
@test_approx_eq wmean(1:3, [0.5, 0.5, 0.0]) 1.5
@test_approx_eq wmean(1:3, [0.0, 0.5, 0.5]) 2.5
