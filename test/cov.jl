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

if VERSION < v"0.5.0-dev+679"
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
else
    @test_approx_eq scattermat(X) S1
    @test_approx_eq scattermat(X, 2) S2

    @test_approx_eq StatsBase.scattermatm(X, 0) Sz1
    @test_approx_eq StatsBase.scattermatm(X, 0, 2) Sz2

    @test_approx_eq StatsBase.scattermatm(X, mean(X,1)) S1
    @test_approx_eq StatsBase.scattermatm(X, mean(X,2), 2) S2

    @test_approx_eq StatsBase.scattermatm(X, zeros(1,8)) Sz1
    @test_approx_eq StatsBase.scattermatm(X, zeros(3), 2) Sz2

    ## weighted scatter mat

    @test_approx_eq scattermat(X, wv1) S1w
    @test_approx_eq scattermat(X, wv2, 2) S2w

    @test_approx_eq StatsBase.scattermatm(X, 0, wv1) Sz1w
    @test_approx_eq StatsBase.scattermatm(X, 0, wv2, 2) Sz2w

    @test_approx_eq StatsBase.scattermatm(X, mean(X, wv1, 1), wv1) S1w
    @test_approx_eq StatsBase.scattermatm(X, mean(X, wv2, 2), wv2, 2) S2w

    @test_approx_eq StatsBase.scattermatm(X, zeros(1,8), wv1) Sz1w
    @test_approx_eq StatsBase.scattermatm(X, zeros(3), wv2, 2) Sz2w
end

# weighted covariance

if VERSION < v"0.5.0-dev+679"
    @test_approx_eq cov(X, wv1) S1w ./ sum(wv1)
    @test_approx_eq cov(X, wv2; vardim=2) S2w ./ sum(wv2)

    @test_approx_eq cov(X, wv1; mean=0) Sz1w ./ sum(wv1)
    @test_approx_eq cov(X, wv2; mean=0, vardim=2) Sz2w ./ sum(wv2)

    @test_approx_eq cov(X, wv1; mean=mean(X, wv1, 1)) S1w ./ sum(wv1)
    @test_approx_eq cov(X, wv2; mean=mean(X, wv2, 2), vardim=2) S2w ./ sum(wv2)

    @test_approx_eq cov(X, wv1; mean=zeros(1,8)) Sz1w ./ sum(wv1)
    @test_approx_eq cov(X, wv2; mean=zeros(3), vardim=2) Sz2w ./ sum(wv2)
else
    @test_approx_eq cov(X, wv1) S1w ./ sum(wv1)
    @test_approx_eq cov(X, wv2, 2) S2w ./ sum(wv2)

    @test_approx_eq Base.covm(X, 0, wv1) Sz1w ./ sum(wv1)
    @test_approx_eq Base.covm(X, 0, wv2, 2) Sz2w ./ sum(wv2)

    @test_approx_eq Base.covm(X, mean(X, wv1, 1), wv1) S1w ./ sum(wv1)
    @test_approx_eq Base.covm(X, mean(X, wv2, 2), wv2, 2) S2w ./ sum(wv2)

    @test_approx_eq Base.covm(X, zeros(1,8), wv1) Sz1w ./ sum(wv1)
    @test_approx_eq Base.covm(X, zeros(3), wv2, 2) Sz2w ./ sum(wv2)
end

# mean_and_cov
if VERSION < v"0.5.0-dev+679"
    (m, C) = mean_and_cov(X; vardim=1)
    @test m == mean(X, 1)
    @test C == cov(X; vardim=1)

    (m, C) = mean_and_cov(X; vardim=2)
    @test m == mean(X, 2)
    @test C == cov(X; vardim=2)

    (m, C) = mean_and_cov(X, wv1; vardim=1)
    @test m == mean(X, wv1, 1)
    @test C == cov(X, wv1; vardim=1)

    (m, C) = mean_and_cov(X, wv2; vardim=2)
    @test m == mean(X, wv2, 2)
    @test C == cov(X, wv2; vardim=2)
else
    (m, C) = mean_and_cov(X, 1)
    @test m == mean(X, 1)
    @test C == cov(X, 1)

    (m, C) = mean_and_cov(X, 2)
    @test m == mean(X, 2)
    @test C == cov(X, 2)

    (m, C) = mean_and_cov(X, wv1, 1)
    @test m == mean(X, wv1, 1)
    @test C == cov(X, wv1, 1)

    (m, C) = mean_and_cov(X, wv2, 2)
    @test m == mean(X, wv2, 2)
    @test C == cov(X, wv2, 2)
end
