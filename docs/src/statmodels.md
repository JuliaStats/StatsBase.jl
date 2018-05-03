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
fit
fit!
informationmatrix
isfitted
islinear
loglikelihood
mss
nobs
nulldeviance
r2
rss
score
stderror
vcov
weights(::StatisticalModel)
```

`RegressionModel` extends `StatisticalModel` by implementing the following additional methods.
```@docs
dof_residual
fitted
leverage
meanresponse
modelmatrix
response
predict
predict!
residuals
```
