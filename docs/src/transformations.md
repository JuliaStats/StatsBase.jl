# Data Transformations

In general, data transformations change raw feature vectors into
a representation that is more suitable for various estimators.

Many learning algorithms benefit from standardization of the data set, and if
some outliers are present in the dataset, robust scalers or transformers are
more appropriate.

## Standardization

**Standardization** of dataset is a common requirement for many machine
learning techniques. These techniques might perform poorly if the individual
features do not more or less look like standard normally distributed data.

Standardization can be perform by using `ZScoreTransform` in the `fit` method.

```@docs
fit(::Type{ZScoreTransform}, X::AbstractArray{<:Real,2}; center::Bool=true, scale::Bool=true)
```

## Scaling

An alternative standardization is scaling features to lie in the unit interval.

Unit range scaling can be perform by using `UnitRangeTransform` in the `fit` method.

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
