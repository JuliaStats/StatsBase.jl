Scalar Statistics
===================

The package implements functions for computing various statistics over an array of scalar real numbers.  

Moments
---------

.. py:function:: skewness(x[, wv])

  Compute the (standardized) `skewness <http://en.wikipedia.org/wiki/Skewness>`_ of ``x``. 

  One can optionally supply a weight vector of type ``WeightVec`` (see :ref:`weightvec`).

.. py:function:: kurtosis(x[, wv])

  Compute the (excessive) `kurtosis <http://en.wikipedia.org/wiki/Kurtosis>`_ of ``x``. 

  One can optionally supply a weight vector of type ``WeightVec`` (see :ref:`weightvec`).


Measurements of Variation
---------------------------

.. py:function:: variation(x)

  Compute the variation of ``x``, *i.e.* the ratio of standard deviation to mean.

.. py:function:: sem(x)

  Compute the standard error of the mean for ``x``, *i.e.* ``sqrt(var(x) / length(x))``.

.. py:function:: mad(x)

  Compute the `median absolute deviation <http://en.wikipedia.org/wiki/Median_absolute_deviation>`_ of ``x``.

.. py:function:: middle(a, b)

  Compute the middle between ``a`` and ``b``, *i.e.* ``(a + b) / 2``.

.. py:function:: midrange(x)

  Compute the middle between ``minimum(x)`` and ``maximum(x)``.  

.. py:function:: range(x)

  Compute the difference between ``maximum(x)`` and ``minimum(x)``.


Quantile and Friends
---------------------

.. py:function:: percentile(x, p)

  Compute quantiles using percentage, *i.e.* ``quantile(x, p / 100)``.

.. py:function:: iqr(x)

  Compute the `interquartile range <http://en.wikipedia.org/wiki/Interquartile_range>`_ of ``x``, *i.e.* ``quantile(x, 0.75) - quantile(x, 0.25)``.

.. py:function:: nquantile(x)

  Compute quantiles at ``[0:n]/n``. For example, ``nquantiles(x, 5)`` returns a vector of quantiles, respectively at ``0.0, 0.2, 0.4, 0.6, 0.8, 1.0``.

.. py:function:: quantile(x)    

  Extended method of *quantile*. Equivalent to ``nquantile(x, 4)``, which returns a vector of quantiles at ``0.0, 0.25, 0.50, 0.75, 1.0``. 


Mode and Modes
---------------

.. py:function:: mode(x)  

  Return the mode of ``x``, one of the numbers that appear the most times in ``x``. 

.. py:function:: modes(x)

  Return a vector of all modes in ``x``. Even if ``x`` has only a single mode, it returns a vector that contains that mode.


Summary of Statistics
-----------------------

.. py:function:: summarystats(x)

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

.. py:function:: describe(x)  

  Print a summary of stats of ``x``. 

