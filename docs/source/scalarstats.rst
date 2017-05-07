Scalar Statistics
===================

The package implements functions for computing various statistics over an array of scalar real numbers.

Moments
---------

.. function:: var(x, w, [dim][; mean=..., corrected=...])

  Compute the variance of a real-valued array ``x``, optionally over a dimension ``dim``.
  Observations in ``x`` are weighted using weight vector ``w``.
  The uncorrected (when ``corrected=false``) sample variance is defined as:

  :math:`\frac{1}{\sum{w}} \sum_{i=1}^n {w_i\left({x_i - m}\right)^2 }`

  where ``n`` is the length of the input and ``m`` is the mean.
  The unbiased estimate (when ``corrected=true``) of the population variance is computed by
  replacing :math:`\frac{1}{\sum{w}}` with a factor dependent on the type of weights used:

    * ``AnalyticWeights``: :math:`\frac{1}{\sum w - \sum {w^2} / \sum w}`
    * ``FrequencyWeights``: :math:`\frac{1}{\sum{w} - 1}`
    * ``ProbabilityWeights``: :math:`\frac{n}{(n - 1) \sum w}` where ``n`` equals ``count(!iszero, w)``
    * ``Weights``: ``ArgumentError`` (bias correction not supported)

.. function:: std(v, w, [dim][; mean=..., corrected=...])

  Compute the standard deviation of a real-valued array ``x``, optionally over a dimension ``dim``.
  Observations in ``x`` are weighted using weight vector ``w``.
  The uncorrected (when ``corrected=false``) sample standard deviation is defined as:

  :math:`\sqrt{\frac{1}{\sum{w}} \sum_{i=1}^n {w_i\left({x_i - m}\right)^2 }}`

  where ``n`` is the length of the input and ``m`` is the mean.
  The unbiased estimate (when ``corrected=true``) of the population standard deviation is
  computed by replacing :math:`\frac{1}{\sum{w}}` with a factor dependent on the type of
  weights used:

    * ``AnalyticWeights``: :math:`\frac{1}{\sum w - \sum {w^2} / \sum w}`
    * ``FrequencyWeights``: :math:`\frac{1}{\sum{w} - 1}`
    * ``ProbabilityWeights``: :math:`\frac{n}{(n - 1) \sum w}` where ``n`` equals ``count(!iszero, w)``
    * ``Weights``: ``ArgumentError`` (bias correction not supported)

.. function:: mean_and_var(x[, w][, dim][; corrected=...])

  Jointly compute the mean and variance of a real-valued array ``x``, optionally over a dimension ``dim``, as a tuple.
  Observations in ``x`` can be weighted using weight vector ``w``.
  Finally, bias correction is be applied to the variance calculation if ``corrected=true``.
  See ``var`` documentation for more details.

.. function:: mean_and_std(x[, w][, dim][; corrected=...])

  Jointly compute the mean and standard deviation of a real-valued array ``x``, optionally over a dimension ``dim``, as a tuple.
  A weighting vector ``w`` can be specified to weight the estimates.
  Finally, bias correction is applied to the standard deviation calculation if ``corrected=true``.
  See ``std`` documentation for more details.

.. function:: skewness(x[, wv])

  Compute the (standardized) `skewness <http://en.wikipedia.org/wiki/Skewness>`_ of ``x``.

  One can optionally supply a weight vector of type ``AbstractWeights`` (see :ref:`weightvec`).

.. function:: kurtosis(x[, wv])

  Compute the (excessive) `kurtosis <http://en.wikipedia.org/wiki/Kurtosis>`_ of ``x``.

  One can optionally supply a weight vector of type ``AbstractWeights`` (see :ref:`weightvec`).

.. function:: moment(x, k[, m][, wv])

  Compute the ``k``-th order central moment of the values in `x`. It is the sample mean of
  ``(x - mean(x)).^k``.

  One can optionally supply the center ``m``, and/or a weight vector of type ``AbstractWeights`` (see :ref:`weightvec`).


Measurements of Variation
---------------------------

.. function:: span(x)

  Get the range ``minimum(x):maximum(x)``.

  **Note:** Here, the minimum and maximum of ``x`` are computed in one-pass using ``extrema``.

.. function:: variation(x)

  Compute the variation of ``x``, *i.e.* the ratio of standard deviation to mean.

.. function:: sem(x)

  Compute the standard error of the mean for ``x``, *i.e.* ``sqrt(var(x) / length(x))``.

.. function:: mad(x[, center][; constant=1.4826])

  Compute the `median absolute deviation <http://en.wikipedia.org/wiki/Median_absolute_deviation>`_ of ``x``.

  One can optionally supply the ``center``. By default, ``constant=1.4826`` for consistent estimation of the standard deviation of a normal distribution.


Z-scores
----------

.. function:: zscore(X, μ, σ)

    Compute the Z-scores, given the mean ``μ`` and standard deviation ``σ``, which is defined as ``(x - μ) / σ``.

    This function returns an array ``Z`` of the same size as ``X``.

    Here, ``μ`` and ``σ`` should be both scalars or both arrays. The computation is broadcasting.
    In particular, when ``μ`` and ``σ`` are arrays, they should have the same size, and
    ``size(μ, i) == 1  || size(μ, i) == size(X, i)`` for each dimension.

.. function:: zscore!(X, μ, σ)

    Compute the Z-scores inplace, given the mean ``μ`` and standard deviation ``σ``.

.. function:: zscore!(Z, X, μ, σ)

    Compute the Z-scores, given the mean ``μ`` and standard deviation ``σ``, and write the results to a pre-allocated
    array ``Z``. Here, ``Z`` and ``X`` should have the same size.

.. function:: zscore(X)

    Compute the Z-scores for ``X``. The mean and standard deviation will be computed within the function.

.. function:: zscore(X, dim)

    Compute the Z-scores for ``X`` along a specific dimension. The mean and standard deviation will be computed within the function.



Entropy and Friends
---------------------

.. function:: entropy(p)

   Compute the entropy of the probability vector ``p`` using natural logarithms (units of nats).

.. function:: entropy(p, b)

   Compute the entropy of the probability vector ``p`` using logarithms of base ``b`` (e.g. ``entropy(p,2)`` returns the entropy in bits).

.. function:: crossentropy(p, q)

   Compute the cross entropy between two probability vectors ``p`` and ``q``.

.. function:: crossentropy(p, q, b)

   Compute the cross entropy between ``p`` and ``q`` using logarithms of base ``b``.

.. function:: kldivergence(p, q)

   Compute the Kullback-Leibler divergence between ``p`` and ``q``.

.. function:: kldivergence(p, q, b)

   Compute the Kullback-Leibler divergence between ``p`` and ``q`` using logarithms of base ``b``.

Quantile and Friends
---------------------

.. function:: percentile(x, p)

  Compute quantiles using percentage, *i.e.* ``quantile(x, p / 100)``.

.. function:: iqr(x)

  Compute the `interquartile range <http://en.wikipedia.org/wiki/Interquartile_range>`_ of ``x``, *i.e.* ``quantile(x, 0.75) - quantile(x, 0.25)``.

.. function:: nquantile(x, n)

  Compute quantiles at ``[0:n]/n``. For example, ``nquantiles(x, 5)`` returns a vector of quantiles, respectively at ``0.0, 0.2, 0.4, 0.6, 0.8, 1.0``.

.. function:: quantile(x)

  Extended method of *quantile*. Equivalent to ``nquantile(x, 4)``, which returns a vector of quantiles at ``0.0, 0.25, 0.50, 0.75, 1.0``.

.. function:: median(x, w)

  Compute the weighted median of ``x``, using weights given by a weight vector ``w`` (of type ``AbstractWeights``).  The weight and data vectors must have the same length.  The weighted median :math:`x_k` is the element of ``x`` that satisfies :math:`\sum_{x_i < x_k} w_i \le \frac{1}{2} \sum_{j} w_j` and :math:`\sum_{x_i > x_k} w_i \le \frac{1}{2} \sum_{j} w_j`.  If a weight has value zero, then its associated data point is ignored.  If none of the weights are positive, an error is thrown.  ``NaN`` is returned if ``x`` contains any ``NaN`` values.  An error is raised if ``w`` contains any ``NaN`` values.

  **Examples:**

  .. code-block:: julia

    w = rand(n)
    xk = median(x, weights(w))

.. function:: quantile(x, w, p)

  Compute the weighted quantiles of a vector ``x`` at a specified set of probability values ``p``, using weights given by a weight vector ``w`` (of type ``AbstractWeights``).  Weights must not be negative. The weights and data vectors must have the same length. The quantile for :math:`p` is defined as follows.  Denoting :math:`S_k = (k-1)w_k + (n-1) \sum_{i<k}w_i`, define :math:`x_{k+1}` the smallest element of ``x`` such that :math:`S_{k+1}/S_{n}` is strictly superior to :math:`p`. The function returns :math:`(1-\gamma) x_k + \gamma x_{k+1}` with  :math:`\gamma = (pS_n- S_k)/(S_{k+1}-S_k)`. This corresponds to  R-7, Excel, SciPy-(1,1), Maple-6 when ``w`` is one (see https://en.wikipedia.org/wiki/Quantile).

Mode and Modes
---------------

.. function:: mode(x)

  Return the mode of ``x``, one of the numbers that appear the most times in ``x``.

.. function:: modes(x)

  Return a vector of all modes in ``x``. Even if ``x`` has only a single mode, it returns a vector that contains that mode.


Summary of Statistics
-----------------------

.. function:: summarystats(x)

  Compute a set of statistics over ``x`` and return a struct of type ``SummaryStats`` defined as below:

  .. code-block:: julia

    immutable SummaryStats{T<:AbstractFloat}
        mean::T
        min::T
        q25::T
        median::T
        q75::T
        max::T
    end

.. function:: describe(x)

  Print a summary of stats of ``x``.
