# Scatter Matrix and Covariance

This package implements functions for computing scatter matrix, as well as weighted covariance matrix.

```@docs
scattermat
cov
cov(::CovarianceEstimator, ::AbstractVector)
cov(::CovarianceEstimator, ::AbstractVector, ::AbstractVector)
cov(::CovarianceEstimator, ::AbstractMatrix)
var(::CovarianceEstimator, ::AbstractVector)
std(::CovarianceEstimator, ::AbstractVector)
cor
mean_and_cov
cov2cor
StatsBase.cov2cor!
cor2cov
StatsBase.cor2cov!
CovarianceEstimator
SimpleCovariance
```
