Scalar Statistics
===================

The package implements functions for computing various statistics over an array of scalar real numbers.  

Moments
---------

.. function:: var(x, wv[; mean=...])

  Compute weighted variance. 

  One can set the keyword argument ``mean``, which can be either ``nothing`` (to compute the mean value within the function), ``0``, or a pre-computed mean value.

  **Note:** the result is normalized by ``sum(wv)`` without correction.

.. function:: var(x, wv, dim[; mean=...])

  Weighted variance along a specific dimension.

.. function:: std(x, wv[; mean=...])

  Compute weighted standard deviation. 

  One can set the keyword argument ``mean``, which can be either ``nothing`` (to compute the mean value within the function), ``0``, or a pre-computed mean value.

.. function:: std(x, wv, dim[; mean=...])

  Weighted standard deviation along a specific dimension.

.. function:: mean_and_var(x[, wv][, dim])

  Jointly compute the mean and variance of ``x``. 

.. function:: mean_and_std(x[, wv][, dim])

  Jointly compute the mean and standard deviation of ``x``.

.. function:: skewness(x[, wv])

  Compute the (standardized) `skewness <http://en.wikipedia.org/wiki/Skewness>`_ of ``x``. 

  One can optionally supply a weight vector of type ``WeightVec`` (see :ref:`weightvec`).

.. function:: kurtosis(x[, wv])

  Compute the (excessive) `kurtosis <http://en.wikipedia.org/wiki/Kurtosis>`_ of ``x``. 

  One can optionally supply a weight vector of type ``WeightVec`` (see :ref:`weightvec`).

.. function:: moment(x, k[, m][, wv])

  Compute the ``k``-th order central moment of the values in `x`. It is the sample mean of 
  ``(x - mean(x)).^k``. 

  One can optionally supply the center ``m``, and/or a weight vector of type ``WeightVec`` (see :ref:`weightvec`).


Measurements of Variation
---------------------------

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

   Compute the entropy of the probability vector ``p``.

.. function:: crossentropy(p, q)

   Compute the cross entropy between two probability vectors ``p`` and ``q``.

.. function:: kldivergence(p, q)

   Compute the Kullback-Leibler divergence between ``p`` and ``q``.


Quantile and Friends
---------------------

.. function:: percentile(x, p)

  Compute quantiles using percentage, *i.e.* ``quantile(x, p / 100)``.

.. function:: iqr(x)

  Compute the `interquartile range <http://en.wikipedia.org/wiki/Interquartile_range>`_ of ``x``, *i.e.* ``quantile(x, 0.75) - quantile(x, 0.25)``.

.. function:: nquantile(x)

  Compute quantiles at ``[0:n]/n``. For example, ``nquantiles(x, 5)`` returns a vector of quantiles, respectively at ``0.0, 0.2, 0.4, 0.6, 0.8, 1.0``.

.. function:: quantile(x)    

  Extended method of *quantile*. Equivalent to ``nquantile(x, 4)``, which returns a vector of quantiles at ``0.0, 0.25, 0.50, 0.75, 1.0``. 

.. function:: median(x, w[; checknan=true])

  Compute the weighted median of ``x``, using weights given by a weight vector ``w`` (of type ``WeightVec``).  The weight and data vectors must have the same length.  The weighted median :math:`x_k` is the element of ``x`` that satisfies :math:`\sum_{x_i < x_k} w_i \le \frac{1}{2} \sum_{j} w_j` and :math:`\sum_{x_i > x_k} w_i \le \frac{1}{2} \sum_{j} w_j`.  If a weight has value zero, then its associated data point is ignored.  If none of the weights are positive, an error is thrown.  ``NaN`` is returned if ``x`` contains any ``NaN`` values.  An error is raised if ``w`` contains any ``NaN`` values.

  **Examples:**

  .. code-block:: julia

    w = rand(n)
    xk = median(x, weights(w))


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

    immutable SummaryStats{T<:FloatingPoint}
        mean::T
        min::T
        q25::T    
        median::T    
        q75::T
        max::T
    end

.. function:: describe(x)  

  Print a summary of stats of ``x``. 

