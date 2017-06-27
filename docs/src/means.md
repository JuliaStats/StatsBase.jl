# Mean Functions

The package provides functions to compute means of different kinds.

```@docs
geomean
harmmean
genmean
```

The `mean` and `mean!` functions are also extended to accept a weight vector of type
`AbstractWeights` to compute weighted mean.

```@docs
Base.mean(A::AbstractArray, w::AbstractWeights)
Base.mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights, dim::Int)
```
