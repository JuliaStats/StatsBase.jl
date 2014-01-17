## StatsBase.jl

Functions for statistics in Julia. [![Build Status](https://travis-ci.org/JuliaStats/StatsBase.jl.png?branch=master)](https://travis-ci.org/JuliaStats/StatsBase.jl)

## List of Functions

The following is a list of functions by categories.

* [Weight Vectors](https://github.com/JuliaStats/StatsBase.jl#weight-vectors)
* [Means](https://github.com/JuliaStats/StatsBase.jl#means)
* [Scalar Statistics](https://github.com/JuliaStats/StatsBase.jl#scalar-statistics)
* [Counts](https://github.com/JuliaStats/StatsBase.jl#counts)
* [Ranking](https://github.com/JuliaStats/StatsBase.jl#ranking)
* [Covariances and Correlations](https://github.com/JuliaStats/StatsBase.jl#covariances-and-correlations)
* [Empirical Estimation](https://github.com/JuliaStats/StatsBase.jl#empirical-estimation)
* [Miscelleneous Functions](https://github.com/JuliaStats/StatsBase.jl#miscelleneous-functions)


#### Weight Vectors

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

* **skewness**(x[, wv])

  Compute the (standardized) skewness of ``x``. [Wikipedia](http://en.wikipedia.org/wiki/Skewness)

  One can optionally supply a weight vector of type ``WeightVec``, as ``weights(w)``.

* **kurtosis**(x[, wv])

  Compute the (excessive) kurtosis of ``x``. [Wikipedia](http://en.wikipedia.org/wiki/Kurtosis)

  One can optionally supply a weight vector of type ``WeightVec``, as ``weights(w)``.

* **variation**(x)

  Compute the variation of ``x``, *i.e.* ratio of standard deviation to mean.

* **sem**(x)

  Compute the standard error of the mean for ``x``, *i.e.* ``sqrt(var(x) / length(x))``.

* **mad**(x)

  Compute the median absolute deviation of ``x``. [Wikipedia](http://en.wikipedia.org/wiki/Median_absolute_deviation)

* **middle**(a, b)

  Compute the middle between ``a`` and ``b``, *i.e.* ``(a + b) / 2``.

* **midrange**(x)

  Compute the middle between ``minimum(x)`` and ``maximum(x)``.  

* **range**(x)

  Compute the difference between ``maximum(x)`` and ``minimum(x)``.

* **percentile**(x, p)

  Compute quantiles using percentage, *i.e.* ``quantile(x, p / 100)``.

* **iqr**(x)

  Compute the interquartile range, *i.e.* ``quantile(x, 0.75) - quantile(x, 0.25)``. [Wikipedia](http://en.wikipedia.org/wiki/Interquartile_range)

* **nquantile**(x)

  Compute quantiles at ``[0:n]/n``. For example, ``nquantiles(x, 5)`` returns a vector of quantiles, respectively at ``0.0, 0.2, 0.4, 0.6, 0.8, 1.0``.

* **quantile**(x)    

  Extended method of *quantile*. Equivalent to ``nquantile(x, 4)``, which returns a vector of quantiles at ``0.0, 0.25, 0.50, 0.75, 1.0``. 

* **mode**(x)  

  Return the mode of ``x``, one of the numbers that appear the most times in ``x``. 

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


#### Counts

* **counts**(x, a:b[, weights(w)])

  Count the number of times (or total weights if a weight vector is given) values in ``a:b`` appear in array ``x``. 

  For example, ``counts([1, 2, 2, 2, 3, 3], 1:3)`` returns ``[1, 3, 2]``. 

  Convenient forms:

  - ``counts(x, k)`` is equivalent to ``counts(x, 1:k)`` when ``k`` is an integer.
  - ``counts(x, k, w)`` is equivalent to ``counts(x, 1:k, w)``.


* **counts**(x, y, (a:b, c:d)[, weights(w)])

  Count the number of times (or total weights if a weight vector is given) pairs of values in ``a:b`` and ``c:d`` that appear in arrays ``x`` and ``y``.

  For example, ``counts([1, 2, 1, 1, 2], [1, 2, 3, 3, 2], (1:2, 1:3))`` returns an array ``r`` of size ``(2, 3)`` as 
  ```julia
  1 0 2
  0 0 2
  ```
  In particular, we have this invariant: ``r[i, j] == sum(x .== (a:b)[i] & y .== (c:d)[j])``.

  This function is useful in various applications, *e.g.* computing confusion matrix in evaluating the performance of a classifier.

  Convenient forms:

  - ``counts(x, y, a:b)`` is equivalent to ``counts(x, y, (a:b, a:b))``.
  - ``counts(x, y, a:b, w)`` is equivalent to ``counts(x, y, (a:b, a:b), w)``.
  - ``counts(x, y, (kx, ky))`` is equivalent to ``counts(x, y, (1:kx, 1:ky))``.
  - ``counts(x, y, (kx, ky), w)`` is equivalent to ``counts(x, y, (1:kx, 1:ky), w)``.
  - ``counts(x, y, k)`` is equivalent to ``counts(x, y, (1:k, 1:k)``.
  - ``counts(x, y, k, w)`` is equivalent to ``counts(x, y, (1:k, 1:k), w)``.


* **proportions**(x, a:b[, weights(w)])  

  Compute the proportions of values in ``a:b`` with respect to ``x``. Equivalent to ``counts(x, a:b) / length(x)``. 

* **proportions**(x, y, (a:b, c:d)[, weights(w)])

  Equivalent to ``counts(x, y, (a:b, c:d)) / length(x)``.

**Note:** all convenient forms for the function ``counts`` are also available to the function ``proportions``.  

* **addcounts!**(r, x, a:b[, weights(w)])

  Adds the counts of values in ``x`` to an accumulating array ``r``.

* **addcounts!**(r, x, y, (a:b, c:d)[, weights(w)])

  Adds the counts of pairs in ``x`` and ``y`` to an accumulating matrix ``r``.  

**Note:** For the functions above, when they encounter an value (or pair) that does not fall in the specified range, they simply ignore it (without raising an error or warning).

* **countmap**(x[, weights(w)])

  Return a dictionary that maps distinct values in ``x`` to their counts (or total weights).

* **proportionmap**(x[, weights(w)])

  Return a dictionary that maps distinct values in ``x`` to their proportions. 

* **addcounts!**(dict, x[, weights(w)])

  Add counts based on ``x`` to a count map. New entries will be added if new values come up.


#### Ranking

This package implements various stragies for computing ranks. Please refer to [Wikipedia](http://en.wikipedia.org/wiki/Ranking) for details of these ranking strategies.

* **ordinalrank**(x)

  Compute ordinal ranking (*1 2 3 4* ranking) for ``x``. 

* **competerank**(x)

  Compute competition ranking (*1 2 2 4* ranking) for ``x``.

* **denserank**(x)

  Compute dense ranking (*1 2 2 3* ranking) for ``x``.

* **tiedrank**(x)

  Compute tied (fractional) ranking (*1 2.5 2.5 4* ranking) for ``x``.


#### Covariances and Correlations

* **autocov**(x, lags[; demean={true}/false])

  Compute auto-covariance of ``x`` at specified lags. If ``x`` is a vector, it returns a vector of th same length of ``lags``. If ``x`` is a matrix, it returns a matrix of size ``(length(lags), size(x,2))``, where each column in the result corresponding to a column in ``x``. 

  Here, ``demean`` is a keyword argument (default value is ``true``), which means that the function will subtract each ``x`` from its mean before computing the results. Otherwise, ``x`` is considered as having been centered. 

* **autocov**(x[; demean={true}/false])

  Compute auto-covariance at default lags.  

* **autocov!**(r, x, lags[; demean={true}/false])

  Write the computed auto-covariance to ``r``.   
 
* **autocor**(x, lags[; demean={true}/false])

  Compute auto-correlation of ``x`` at specified lags. 

* **autocor**(x[; demean={true}/false])

  Compute auto-correlation at default lags.  

* **autocor!**(r, x, lags[; demean={true}/false])

  Write the computed auto-correlation to ``r``.   

* **crosscov**(x, y, lags[; demean={true}/false])

  Compute cross-covariance of ``x`` at specified lags. If both ``x`` and ``y`` are vectors, it returns a vector of th same length of ``lags``. Otherwise, it computes cross covariances between each pairs of columns in ``x`` and ``y``.

* **crosscov**(x, y[; demean={true}/false])

  Compute cross-covariance at default lags.  

* **crosscov!**(r, x, y, lags[; demean={true}/false])

  Write the computed cross-covariance to ``r``.   
 
* **crosscor**(x, y, lags[; demean={true}/false])

  Compute cross-correlation of ``x`` at specified lags. 

* **crosscor**(x, y[; demean={true}/false])

  Compute cross-correlation at default lags.  

* **crosscor!**(r, x, y, lags[; demean={true}/false])

  Write the computed cross-correlation to ``r``.   

* **pacf**(x, lags[; method={:regression}/:yulewalker])

  Compute partial auto-correlation of ``x`` at specified lags. If ``x`` is a vector, it returns a vector of th same length of ``lags``. If ``x`` is a matrix, it returns a matrix of size ``(length(lags), size(x,2))``, where each column in the result corresponding to a column in ``x``. 

  Here, ``method`` is a keyword argument to specify the choice of algorithm, which can be either ``:regresion`` or ``:yulewalker``. The default value is ``:regression``. 


* **pacf!**(r, x, lags[; method={:regression}/:yulewalker])

  Write the computed partial auto-correlation to ``r``.


* **corspearman**(x, y)

  Spearman's rank correlation. Here, ``x`` and ``y`` can be either real vectors or matrices. When ``xx`` and ``y`` are matrices, it computes the correlations between their columns (pairwisely).

* **corkendall**(x, y)  	

  Kendall's rank correlation. Here, ``x`` and ``y`` can be either real vectors or matrices. When ``xx`` and ``y`` are matrices, it computes the correlations between their columns (pairwisely).


#### Empirical Estimation

* **ecdf**(x)

  Return an empirical cumulative distribution function based on a vector of samples given in ``x``. 

  **Note:** this is a higher-level function that returns a function, which can then be applied to evaluate CDF values on other samples.

* **kde**(data[; width=NaN, npoints=2048]) 

  Kernel density estimation. This function returns an instance of ``UnivariateKDE`` defined as below:
  ```julia
  immutable UnivariateKDE
      x::Vector{Float64}
      density::Vector{Float64}
  end
  ```

  It has two keyword arguments: ``width`` and ``npoints``, respectively for specifying the kernel width and the number of sample points. When ``width`` is set to ``NaN`` (which is the default case), the function will choose its value automatically.


* **kde**(x, y[; width=NaN, resolution=25])

  Two-dimensional kernel density estimation. This function returns an instance of ``BivariateKDE`` defined as below:
  ```julia
  immutable BivariateKDE
      x::Vector{Float64}
      y::Vector{Float64}
      density::Matrix{Float64}
  end
  ```

#### Miscelleneous Functions

* **rle**(x)

  Run-length encoding of ``x``. It returns ``(vals, lens)``, a sequence of values and their corresponding chunk length. [Wikipedia](http://en.wikipedia.org/wiki/Run-length_encoding).

  Example:
  ```julia
  julia> rle([1,1,1,2,2,3,3,3,3,2,2,2])
  ([1,2,3,2],[3,2,4,3])
  ```

* **inverse_rle**(vals, lens)

  Inversed run-length encoding. It takes the results of ``rle`` and reconstructs the original sequence. 

* **levelsmap**(x)

  Construct a dictionary that maps each of the ``n`` distinct values in ``x`` to a number between ``1`` and ``n``.

* **indexmap**(x)

  Construct a dictionary that maps each distinct value in ``x`` to its first index.

* **findat**(a, x)

  For each element in ``x``, find its first index in ``a``. If the value does not appear in ``a``, the correspinding index is ``0``. 

  Example:
  ```julia
  julia> findat([2,4,6], [2,3,4])
  3-element Array{Int64,1}:
   1
   0
   2
  ```

* **findat!**(r, a, x)

  Write the results of ``findat(a, x)`` to a pre-allocated array ``r``.

* **indicatormat**(x, k[; sparse=false])  

  Construct a boolean matrix ``r`` of size ``(k, length(x))`` such that ``r[x[i], i] = true`` and all other elements are set to ``false``.

  The keyword argument ``sparse`` controls whether to construct a sparse matrix. By default, it is false. 

  Example:
  ```julia
  julia> indicatormat([1 2 2], 2)
  2x3 Array{Bool,2}:
    true  false  false
   false   true   true
  ```

* **indicatormat**(x, c[; sparse=false])

  Construct a boolean matrix ``r`` of size ``(length(c), length(x))``. Let ``ci`` be the index of ``x[i]`` in ``c``, then ``r[ci, i] = true`` and all other elements are zero. 

  The keyword argument ``sparse`` controls whether to construct a sparse matrix. By default, it is false. 

* **indicatormat**(x[; sparse=false])

  Equivalent to ``indicatormap(x, sort(unique(x)); sparse=...)``. 


