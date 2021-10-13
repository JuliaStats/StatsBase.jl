# Data Transformations

In general, data transformations change raw feature vectors into
a representation that is more suitable for various estimators.

## Standardization a.k.a Z-score Normalization

**Standardization**, also known as Z-score normalization, is a common requirement
for many machine learning techniques. These techniques might perform poorly
if the individual features do not more or less look like standard normally
distributed data.

Standardization transforms data points into corresponding standard scores
by subtracting mean and scaling to unit variance.

The **standard score**, also known as Z-score, is the signed number of
standard deviations by which the value of an observation or data point
is above the mean value of what is being observed or measured.

Standardization can be performed using `t = fit(ZScoreTransform, ...)`
followed by `StatsBase.transform(t, ...)` or `StatsBase.transform!(t, ...)`.
`standardize(ZScoreTransform, ...)` is a shorthand to perform both operations
in a single call.

```@docs
fit(::Type{ZScoreTransform}, X::AbstractArray{<:Real,2}; center::Bool=true, scale::Bool=true)
```

## Unit Range Normalization

**Unit range normalization**, also known as min-max scaling, is an alternative
data transformation which scales features to lie in the interval `[0; 1]`.

Unit range normalization can be performed using `t = fit(UnitRangeTransform, ...)`
followed by `StatsBase.transform(t, ...)` or `StatsBase.transform!(t, ...)`.
`standardize(UnitRangeTransform, ...)` is a shorthand to perform both operations
in a single call.

```@docs
fit(::Type{UnitRangeTransform}, X::AbstractArray{<:Real,2}; unit::Bool=true)
```

## Methods
```@docs
StatsBase.transform
StatsBase.transform!
StatsBase.reconstruct
StatsBase.reconstruct!
standardize
```

## Types
```@docs
UnitRangeTransform
ZScoreTransform
```