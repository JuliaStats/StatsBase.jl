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

* **counts**(x, y, [weights(w),] a:b, c:d)

  Count the number of times (or total weights if a weight vector is given) pairs of values in ``a:b`` and ``c:d`` that appear in arrays ``x`` and ``y``.

  For example, ``counts([1, 2, 1, 1, 2], [1, 2, 3, 3, 2], 1:2, 1:3)`` returns an array ``r`` of size ``(2, 3)`` as 
  ```julia
  1 0 2
  0 0 2
  ```
  In particular, we have this invariant: ``r[i, j] == sum(x .== (a:b)[i] & y .== (c:d)[j])``.

  This function is useful in various applications, *e.g.* computing confusion matrix in evaluating the performance of a classifier.

* **proportions**(x, a:b[, weights(w)])  

  Compute the proportions of values in ``a:b`` with respect to ``x``. Equivalent to ``counts(x, a:b) / length(x)``. 

* **proportions**(x, y, [weights(w),] a:b, c:d)

  Equivalent to ``counts(x, y, a:b, c:d) / length(x)``.

* **addcounts!**(r, x, a:b[, weights(w)])

  Adds the counts of values in ``x`` to an accumulating array ``r``.

* **addcounts!**(r, x, y, [weights(w),] a:b, c:d)

  Adds the counts of pairs in ``x`` and ``y`` to an accumulating matrix ``r``.  

**Note:** For the functions above, when they encounter an value (or pair) that does not fall in the specified range, they simply ignore it (without raising an error or warning).

* **countmap**(x[, weights(w)])

  Return a dictionary that maps distinct values in ``x`` to their counts (or total weights).

* **proportionmap**(x[, weights(w)])

  Return a dictionary that maps distinct values in ``x`` to their proportions. 

* **addcounts!**(dict, x[, weights(w)])

  Add counts based on ``x`` to a count map. New entries will be added if new values come up.
  



