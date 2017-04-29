using StatsBase
using Base.Test

X = randn(3, 8)

Z1 = X .- mean(X, 1)
Z2 = X .- mean(X, 2)

w1 = rand(3)
w2 = rand(8)

wv1 = fweights(w1)
wv2 = fweights(w2)

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


@test scattermat(X)    ≈ S1
@test scattermat(X, 2) ≈ S2

@test StatsBase.scattermatm(X, 0)    ≈ Sz1
@test StatsBase.scattermatm(X, 0, 2) ≈ Sz2

@test StatsBase.scattermatm(X, mean(X,1))    ≈ S1
@test StatsBase.scattermatm(X, mean(X,2), 2) ≈ S2

@test StatsBase.scattermatm(X, zeros(1,8))  ≈ Sz1
@test StatsBase.scattermatm(X, zeros(3), 2) ≈ Sz2

## weighted scatter mat

@test scattermat(X, wv1)    ≈ S1w
@test scattermat(X, wv2, 2) ≈ S2w

@test StatsBase.scattermatm(X, 0, wv1)    ≈ Sz1w
@test StatsBase.scattermatm(X, 0, wv2, 2) ≈ Sz2w

@test StatsBase.scattermatm(X, mean(X, wv1, 1), wv1)    ≈ S1w
@test StatsBase.scattermatm(X, mean(X, wv2, 2), wv2, 2) ≈ S2w

@test StatsBase.scattermatm(X, zeros(1,8), wv1)  ≈ Sz1w
@test StatsBase.scattermatm(X, zeros(3), wv2, 2) ≈ Sz2w

# weighted covariance

@test cov(X, wv1, false)    ≈ S1w ./ sum(wv1)
@test cov(X, wv2, 2, false) ≈ S2w ./ sum(wv2)

@test Base.covm(X, 0, wv1, 1, false) ≈ Sz1w ./ sum(wv1)
@test Base.covm(X, 0, wv2, 2, false) ≈ Sz2w ./ sum(wv2)

@test Base.covm(X, mean(X, wv1, 1), wv1, 1, false) ≈ S1w ./ sum(wv1)
@test Base.covm(X, mean(X, wv2, 2), wv2, 2, false) ≈ S2w ./ sum(wv2)

@test Base.covm(X, zeros(1,8), wv1, 1, false) ≈ Sz1w ./ sum(wv1)
@test Base.covm(X, zeros(3), wv2, 2, false)   ≈ Sz2w ./ sum(wv2)

# mean_and_cov

(m, C) = mean_and_cov(X, 1, false)
@test m == mean(X, 1)
@test C == cov(X, 1, false)

(m, C) = mean_and_cov(X, 2, false)
@test m == mean(X, 2)
@test C == cov(X, 2, false)

(m, C) = mean_and_cov(X, wv1, 1, false)
@test m == mean(X, wv1, 1)
@test C == cov(X, wv1, 1, false)

(m, C) = mean_and_cov(X, wv2, 2, false)
@test m == mean(X, wv2, 2)
@test C == cov(X, wv2, 2, false)
