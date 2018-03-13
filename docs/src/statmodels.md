# Abstraction for Statistical Models

This package defines an abstract type `StatisticalModel`, and an abstract subtype `RegressionModel`.

Particularly, instances of `StatisticalModel` implement the following methods.

```@docs
adjr2
aic
aicc
bic
coef
coefnames
coeftable
confint
deviance
dof
eim
fit
fit!
isfitted
islinear
leverage
loglikelihood
mss
nobs
nulldeviance
oim
r2
rss
score
stderr
vcov
weights(::StatisticalModel)
```

`RegressionModel` extends `StatisticalModel` by implementing the following additional methods.
```@docs
dof_residual
fitted
meanresponse
modelmatrix
response
predict
predict!
residuals
```
