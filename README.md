## Stats.jl

Basic statistics functions for Julia

## List of Functions

* `autocor(A, l)`: Compute the autocorrelation of an array `A` at lag(s) `l`.
* `cor_spearman(x, y)`: Compute Spearman's rank order correlation between `x` and `y`.
* `cov_spearman(x, y)`: Compute Spearman's rank order covariance between `x` and `y`.
* `decile(a)`: Compute the deciles of `a`.
* `describe(a)`: Print out basic summary statistics for `a`.
* `ecdf(a)`: Return the empirical CDF of `a` as a function that can be evaluated anywhere on the real line.
* `findat!(indices, a, b):` Find the indices at which elements of `a` occur in `b`. Uses `0` indices for elements of `a` that do not occur in `b`. This occurs in-place by mutating indices.
* `findat(a, b)`: Find the indices at which elements of `a` occur in `b`. Uses `0` indices for elements of `a` that do not occur in `b`.
* `gmean(a)`: Compute the geometric mean of `a`.
* `hmean(a)`: Compute the harmonic mean of `a`.
* `inverse_rle(r)`: Translate a run length encoding, `r`, into a dense vector.
* `iqr(a)`: Compute the interquantile range of `a`.
* `kurtosis(a)`: Compute the excess kurtosis of `a`.
* `mad(a)`: Compute the median absolute deviation of `a` with a correction factor, which ensures that the MAD will be a consistent estimator of the mean for normally distributed data.
* `midrange(a)`: Compute the mid point of the range of `a` (e.g `(max(a) + min(a) / 2)`).
* `modes(a)`: Compute all modes of `a`. Be warned that every element of an array with no repeated elements is considered a mode.
* `indicators(a)`: Encode categories using one-hot scheme aka one-of-C encoding, indicator matrix or dummy variables. Optionally, you can provide either a range min:max for categories that are encoded as integers (so that min will be encoded as [1 0 ...] and max as [... 0 1]) or a list of possible values, e.g. ["A", "B, "C"].
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
