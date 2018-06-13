# Empirical Estimation

## Histograms

The `Histogram` type represents data that has been tabulated into intervals
(known as *bins*) along the real line, or in higher dimensions, over the real
plane.

Histograms can be fitted to data using the `fit` method.

```@docs
fit(::Type{Histogram}, args...; kwargs...)
```

Additional methods
```@docs
merge!
merge
norm
normalize
normalize!
zero
```

## Empirical Cumulative Distribution Function

```@docs
ecdf
```
