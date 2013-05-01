using Stats
using Base.Test

@test invlogit(-800) > 0.0
@test invlogit(800) <= 1.0
@test isinf(logit(0.0))
@test isinf(logit(1.0))

@test_approx_eq invlogit(logit(0.01)) 0.01
@test_approx_eq invlogit(logit(0.50)) 0.50
@test_approx_eq invlogit(logit(0.99)) 0.99

@test_approx_eq logsumexp([1., 2., 3.]) log(sum(exp([1., 2., 3.])))
@test isfinite(logsumexp([1000., 1000.]))
@test_approx_eq logsumexp([1000., 1000.]) 1000. + log(2.)
