using StatsBase
using Base.Test

X = randn(3, 8)

Z1 = X .- mean(X, 1)
Z2 = X .- mean(X, 2)

w1 = rand(3)
w2 = rand(8)

wv1 = weights(w1)
wv2 = weights(w2)

Z1w = X .- mean(X, wv1, 1)
Z2w = X .- mean(X, wv2, 2)

## reference results

S1 = Z1'Z1
S2 = Z2 * Z2'

Sz1 = X'X
Sz2 = X * X'

S1w = Z1w' * diagm(w1) * Z1w
S2w = Z2w * diagm(w2) * Z2w'

Sz1w = X' * diagm(w1) * X
Sz2w = X * diagm(w2) * X'

## scattermat

@test_approx_eq scattermat(X) S1
@test_approx_eq scattermat(X; vardim=2) S2

@test_approx_eq scattermat(X; mean=0) Sz1
@test_approx_eq scattermat(X; mean=0, vardim=2) Sz2

@test_approx_eq scattermat(X; mean=mean(X,1)) S1
@test_approx_eq scattermat(X; mean=mean(X,2), vardim=2) S2

@test_approx_eq scattermat(X; mean=zeros(1,8)) Sz1
@test_approx_eq scattermat(X; mean=zeros(3), vardim=2) Sz2

## weighted scatter mat

@test_approx_eq scattermat(X, wv1) S1w
@test_approx_eq scattermat(X, wv2; vardim=2) S2w

@test_approx_eq scattermat(X, wv1; mean=0) Sz1w
@test_approx_eq scattermat(X, wv2; mean=0, vardim=2) Sz2w

@test_approx_eq scattermat(X, wv1; mean=mean(X, wv1, 1)) S1w
@test_approx_eq scattermat(X, wv2; mean=mean(X, wv2, 2), vardim=2) S2w

@test_approx_eq scattermat(X, wv1; mean=zeros(1,8)) Sz1w
@test_approx_eq scattermat(X, wv2; mean=zeros(3), vardim=2) Sz2w

# weighted covariance

@test_approx_eq cov(X, wv1) S1w ./ sum(wv1)
@test_approx_eq cov(X, wv2; vardim=2) S2w ./ sum(wv2)

@test_approx_eq cov(X, wv1; mean=0) Sz1w ./ sum(wv1)
@test_approx_eq cov(X, wv2; mean=0, vardim=2) Sz2w ./ sum(wv2)

@test_approx_eq cov(X, wv1; mean=mean(X, wv1, 1)) S1w ./ sum(wv1)
@test_approx_eq cov(X, wv2; mean=mean(X, wv2, 2), vardim=2) S2w ./ sum(wv2)

@test_approx_eq cov(X, wv1; mean=zeros(1,8)) Sz1w ./ sum(wv1)
@test_approx_eq cov(X, wv2; mean=zeros(3), vardim=2) Sz2w ./ sum(wv2)

