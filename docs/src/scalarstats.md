# Scalar Statistics

The package implements functions for computing various statistics over an array of scalar real numbers.

## Weighted sum and mean

```@docs
sum
sum!
wsum
wsum!
mean
mean!
```

## Means

The package provides functions to compute means of different kinds.

```@docs
geomean
harmmean
genmean
```

## Moments and cumulants

```@docs
var
std
mean_and_var
mean_and_std
skewness
kurtosis
moment
cumulant
```

## Measurements of Variation

```@docs
span
variation
sem
mad
mad!
```

## Z-scores

```@docs
zscore
zscore!
```

## Entropy and Related Functions

```@docs
entropy
renyientropy
crossentropy
kldivergence
```

## Quantile and Related Functions

```@docs
percentile
iqr
nquantile
quantile
Statistics.median(v::AbstractVector{<:Real}, w::AbstractWeights{<:Real})
quantilerank
percentilerank
```

## Mode and Modes

```@docs
mode
modes
```

## Summary Statistics

```@docs
summarystats
describe
```

## Reliability Measures

```@docs
cronbachalpha
```
