# Data Transformations

In general, data transformations change raw feature vectors into
a representation that is more suitable for various estimators.

## Standardization

**Standardization** of dataset is a common requirement for many machine
learning techniques. These techniques might perform poorly if the individual
features do not more or less look like standard normally distributed data.

Standardization transforms data points into corresponding standard scores
by removing mean and scaling to unit variance.

The **standard score** is the signed number of standard deviations by which
the value of an observation or data point is above the mean value of what
is being observed or measured.

Standardization can be performed using `fit(ZScoreTransform, ...)`.

```@docs
fit(::Type{ZScoreTransform}, X::AbstractArray{<:Real,2}; center::Bool=true, scale::Bool=true)
```

## Unit range normalization

**Unit range normalization* is an alternative data transformation which scales features
to lie in the interval `[0; 1]`.

Unit range normalization can be performed using `fit(UnitRangeTransform, ...)`.

```@docs
fit(::Type{UnitRangeTransform}, X::AbstractArray{<:Real,2}; unit::Bool=true)
```

## Additional methods
```@docs
StatsBase.transform
StatsBase.transform!
StatsBase.reconstruct
StatsBase.reconstruct!
standardize
```
