## Stats.jl

Functions for statistics in Julia. [![Build Status](https://travis-ci.org/JuliaStats/Stats.jl.png?branch=master)](https://travis-ci.org/JuliaStats/Stats.jl)

## List of Functions

The following is a list of functions by categories.

#### Weight vectors

This package introduce a type ``WeightVec``, which is used to wrap a real vector into a *weight vector*. This is a shallow wrapper, introduced for two purposes:

* indicator that an input vector is a weight vector when it is input as an argument
* maintain the sum of weights, which is usually useful in statistical routines

For example, one can write:
```julia
r = mean(x, weights(w))  # the inline function weights wraps w into a weight vector
```

#### Means

* **geomean**(x)

  Compute the geometric mean of ``x``.

* **harmean**(x)

  Compute the harmonic mean of ``x``.

* **trimean**(x, p)

  Compute the trimmed mean of ``x``, with fraction ``p`` of elements ignored.

  For example, ``trimmean(x, 0.2)`` computes the mean using *80%* of the numbers in ``x``, while ignoring the *10%* largest and *10%* smallest values.

* **wmean**(x, w)  

  Compute weighted mean of ``x`` using weights in ``w``.

* **mean(x, w)**

  The ``mean`` function is also extended to accept a weight vector. Now, one can write ``mean(x, weights(w))`` to compute a weighted mean.

#### Scalar Statistics

The following functions are for computing statistics over an array of scalar real numbers.  

* **skewness**(x)

  Compute the (standardized) skewness of ``x``. [Wikipedia](http://en.wikipedia.org/wiki/Skewness)

* **kurtosis**(x)

  Compute the (excessive) kurtosis of ``x``. [Wikipedia](http://en.wikipedia.org/wiki/Kurtosis)

* **variation**(x)

  Compute the variation of ``x``, *i.e.* ratio of standard deviation to mean.

* **sem**(x)

  Compute the standard error of the mean for ``x``, *i.e.* ``sqrt(var(x) / length(x))``.

* **mad**(x)

  Compute the median absolute deviation of ``x``. [Wikipedia](http://en.wikipedia.org/wiki/Median_absolute_deviation)

* **minmax**(x)

  Compute both maximum and minimum of ``x`` in a single pass, and return them in a pair as ``(min, max)``.

* **middle**(a, b)

  Compute the middle between ``a`` and ``b``, *i.e.* ``(a + b) / 2``.

* **midrange**(x)

  Compute the middle between ``minimum(x)`` and ``maximum(x)``.  

* **range**(x)

  Compute the difference between ``maximum(x)`` and ``minimum(x)``.

* **prctile**(x, p)

  Compute quantiles using percentage, *i.e.* ``quantile(x, p / 100)``.

* **iqr**(x)

  Compute the interquartile range, *i.e.* ``quantile(x, 0.75) - quantile(x, 0.25)``. [Wikipedia](http://en.wikipedia.org/wiki/Interquartile_range)

* **nquantile**(x)

  Compute quantiles at ``[0:n]/n``. For example, ``nquantiles(x, 5)`` returns a vector of quantiles, respectively at ``0.0, 0.2, 0.4, 0.6, 0.8, 1.0``.

* **quantile**(x)    

  Extended method of *quantile*. Equivalent to ``nquantile(x, 4)``, which returns a vector of quantiles at ``0.0, 0.25, 0.50, 0.75, 1.0``. 

* **mode**(x)  

  Return the mode of ``x``, one of the numbers that appear the most times in ``x``. [Wikipedia](http://en.wikipedia.org/wiki/Mode_(statistics))

* **modes**(x)

  Return a vector of all modes in ``x``. Even if ``x`` has only a single mode, it returns a vector that contains that mode.

* **summarystats**(x)

  Compute a set of statistics over ``x`` and return a struct of type ``SummaryStats`` defined as below:

  ```julia
  immutable SummaryStats{T<:FloatingPoint}
      mean::T
      min::T
      q25::T    
      median::T    
      q75::T
      max::T
  end
  ```
  
* **describe**(x)  

  Print a summary of stats of ``x``. 


#### Old documents


* `acf(A, l)`: Compute the autocorrelation of an array `A` at lag(s) `l`.
* `ccf(A, B, l)`: Compute the crosscorrelation of array `A` and `B` at lag(s) `l`.
* `cor_spearman(x, y)`: Compute Spearman's rank order correlation between `x` and `y`.
* `cov_spearman(x, y)`: Compute Spearman's rank order covariance between `x` and `y`.
* `decile(a)`: Compute the deciles of `a`.
* `describe(a)`: Print out basic summary statistics for `a`.
* `durbin(r)`: Solve the positive definite Toeplitz system Tr=-r.
* `ecdf(a)`: Return the empirical CDF of `a` as a function that can be evaluated anywhere on the real line.
* `findat!(indices, a, b):` Find the indices at which elements of `a` occur in `b`. Uses `0` indices for elements of `a` that do not occur in `b`. This occurs in-place by mutating indices.
* `findat(a, b)`: Find the indices at which elements of `a` occur in `b`. Uses `0` indices for elements of `a` that do not occur in `b`.
* `gmean(a)`: Compute the geometric mean of `a`.
* `hmean(a)`: Compute the harmonic mean of `a`.
* `inverse_rle(r)`: Translate a run length encoding, `r`, into a dense vector.
* `iqr(a)`: Compute the interquantile range of `a`.
* `kurtosis(a)`: Compute the excess kurtosis of `a`.
* `levinson(r,b)`: Solve the postitive definite Toeplitz system Tr=b.
* `mad(a)`: Compute the median absolute deviation of `a` with a correction factor, which ensures that the MAD will be a consistent estimator of the mean for normally distributed data.
* `midrange(a)`: Compute the mid point of the range of `a` (e.g `(max(a) + min(a) / 2)`).
* `modes(a)`: Compute all modes of `a`. Be warned that every element of an array with no repeated elements is considered a mode.
* `indicators(a)`: Encode categories using one-hot scheme aka one-of-C encoding, indicator matrix or dummy variables. Optionally, you can provide a list of possible values, e.g. ["A", "B, "C"] or [1:3].
* `pacf(A, l, method)`: Compute partial autocorrelation of an array `A` at lag(s) `l`. The computational method acn either be :regression (default) or :yulewalker.
* `percentile(a)`: Compute the percentiles (0%, 10%, ..., 100%) of `a`.
* `quantile(a)`: Compute any desired quantile of `a`.
* `quartile(a): Compute the quartiles of `a`.
* `quintile(a)`: Compute the quintiles of `a`.
* `minmax(a)`: Compute the min and max of `a`.
* `range(a)`: Compute the length of range, i.e. ``max(a) - min(a)``.
* `rle(a)`: Compute a run-length encoding of `a`.
* `sem(a)`: Compute the standard error of the mean of `a`.
* `skewness(a)`: Compute the skewness of `a`.
* `table(a): Produce a hash table containing counts of the unique elements of `a`.
* `tiedrank(a)`: Compute the rank of `a`.
* `variation(a)`: Compute the coefficient of variation of `a`.
* `weighted_mean(a, w)`: Compute the weighted mean of `a` using weights `w`.
