# Empirical Estimation

## Histograms

Histograms are often times viewed and talked about in the graphical sense. To understand how to create one graphically, check out [this notebook](https://github.com/mitmath/18338/blob/master/notebooks/0.0%20Histogramming.ipynb). In the below examples, we will explain how to create a histogram non-graphically and manipulate the data within it. 

```@docs
Histogram
```

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
