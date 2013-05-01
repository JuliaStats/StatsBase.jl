using Stats
using Base.Test

@test invlogit(-800) > 0.0
@test invlogit(800) <= 1.0
@test isinf(logit(0.0))
@test isinf(logit(1.0))

@test_approx_eq invlogit(logit(0.01)) 0.01
@test_approx_eq invlogit(logit(0.50)) 0.50
@test_approx_eq invlogit(logit(0.99)) 0.99

# logsumexp

@test_approx_eq logsumexp([1., 2., 3.]) log(sum(exp([1., 2., 3.])))
@test isfinite(logsumexp([1000., 1000.]))
@test_approx_eq logsumexp([1000., 1000.]) 1000. + log(2.)

a = [1. 2. 3.; 4. 5. 6.]
r = zeros(3)
for i = 1 : 3
    r[i] = logsumexp(a[:,i])
end
@test_approx_eq logsumexp(a) r

# entropy

@test_approx_eq pventropy([0.5, 0.5]) log(2.)
@test_approx_eq pventropy([0.25, 0.5, 0., 0.25, 0.]) log(8.)/2
a = [0.0 0.2 0.4 0.6 0.8 1.0; 1.0 0.8 0.6 0.4 0.2 0.0]
r = zeros(6)
for i = 1 : 6
    r[i] = pventropy(a[:,i])
end 
@test_approx_eq pventropy(a) r

# softmax

r = exp([1., 2., 3.,]) / sum(exp([1., 2., 3.]))
@test_approx_eq softmax([1., 2., 3.]) r
@test_approx_eq softmax([10001., 10002., 10003.]) r

a = [1. 2. 3.; 4. 5. 6.]
r = zeros(size(a))
for i = 1 : 3
    r[:,i] = softmax(a[:,i])
end
@test_approx_eq softmax(a) r
