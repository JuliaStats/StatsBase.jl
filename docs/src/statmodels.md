# Abstraction for Statistical Models

[StatsAPI.jl](https://github.com/JuliaStats/StatsAPI.jl) defines an abstract type `StatisticalModel`,
and an abstract subtype `RegressionModel`. They are both extended by StatsBase, and documented here.

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
nullloglikelihood
r2
rss
score
stderror
vcov
weights
```

`RegressionModel` extends `StatisticalModel` by implementing the following additional methods.
```@docs
crossmodelmatrix
dof_residual
fitted
leverage
cooksdistance
meanresponse
modelmatrix
response
responsename
predict
predict!
residuals
```

An exception type is provided to signal convergence failures during model estimation:
```@docs
ConvergenceException
```
