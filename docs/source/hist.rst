.. _hist:

Histograms
================

The ``Histogram`` type represents data that has been tabulated into intervals
(known as *bins*) along the real line, or in higher dimensions, over the real
plane.


Fitting
--------------

Histograms can be fitted to data using the ``fit`` method.

.. function:: fit(`Histogram`, data[, weight][, edges])

Arguments:

``data`` 
  is either a vector (for a 1-dimensional histogram), or a tuple of
  vectors of equal length (for an *n*-dimensional histogram).

``weight``
  is an optional ``:ref:`weightvec` WeightVec``` (of the same length as the
  data vectors), denoting the weight each observation contributes to the
  bin. If no weight vector is supples, each observation has weight 1.

``edges``
  is a vector (typically a `Range` object), or tuple of vectors, that gives
  the edges of the bins along each dimension. If no edges are provided, these
  are determined from the data.

Keyword arguments:

``closed=:left/:right`` 
  determines whether the bin intervals are left-closed [a,b), or right-closed
  (a,b] (default = ``:right``).

``nbins``
  if no ``edges`` argument is supplied, the approximate number of bins to use
  along each dimension (can be either a single integer, or a tuple of integers).


.. code-block:: julia

    # Univariate
    h = fit(Histogram, rand(100))
    h = fit(Histogram, rand(100), 0:0.1:1.0)
    h = fit(Histogram, rand(100), nbins=10)
    h = fit(Histogram, rand(100), weights(rand(100)), 0:0.1:1.0)
    h = fit(Histogram, [20], 0:20:100)
    h = fit(Histogram, [20], 0:20:100, closed=:left)
    
    # Multivariate
    h = fit(Histogram, (rand(100),rand(100)))
    h = fit(Histogram, (rand(100),rand(100)),nbins=10)

