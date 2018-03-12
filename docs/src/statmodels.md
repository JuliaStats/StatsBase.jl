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
ess
fit
fit!
isfitted
leverage
linear
loglikelihood
nobs
nulldeviance
r2
rss
stderr
tss
vcov
```

`RegressionModel` extends `StatisticalModel` by implementing the following additional methods.
```@docs
dof_residual
fitted
meanresponse
modelmatrix
response
modelweights
predict
predict!
residuals
```
