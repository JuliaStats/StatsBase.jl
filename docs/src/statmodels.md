# Abstraction for Statistical Models

This package defines an abstract type `StatisticalModel`, and an abstract subtype `RegressionModel`.

Particularly, instances of `StatisticalModel` implement the following methods.

```@docs
adjr2
aic
aicc
bic
coef
coeftable
confint
deviance
dof
fit
fit!
loglikelihood
nobs
nulldeviance
r2
stderr
vcov
```

`RegressionModel` extends `StatisticalModel` by implementing the following additional methods.
```@docs
dof_residual
fitted
model_response
predict
predict!
residuals
```
