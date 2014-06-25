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

.. function:: mad(x)

  Compute the `median absolute deviation <http://en.wikipedia.org/wiki/Median_absolute_deviation>`_ of ``x``.


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

