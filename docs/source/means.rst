Mean Functions
===============

The package provides functions to compute means of different kinds.

.. function:: geomean(x)

  Compute the geometric mean of ``x``.

.. function:: harmmean(x)

  Compute the harmonic mean of ``x``.

.. function:: genmean(x, p)

  Compute the generalized/power mean of ``x`` with exponent ``p``,
  i.e. :math:`\left( \frac{1}{n} \sum_{i=1}^n x_i^p \right)^{\frac{1}{p}}`, where ``n = length(x)``.
  It is taken to be the geometric mean when ``p == 0``.

.. function:: trimmean(x, p)

  Compute the trimmed mean of ``x``, with fraction ``p`` of elements ignored.

  **Examples:**

  .. code-block:: julia

    # computes the mean using 80% of the numbers in ``x``,
    # ignoring 10% of the largest and 10% of the smallest values.

    trimmean(x, 0.2)

.. function:: mean(x, w)

  The ``mean`` function is also extended to accept a weight vector of type ``WeightVec`` (see :ref:`weightvec`) to compute weighted mean.

  **Examples:**

  .. code-block:: julia

    w = rand(n)
    x = mean(x, weights(w))

.. function:: mean(x, w, dim)

  Compute weighted means of ``x`` along a certain dimension (specified by an integer ``dim``). The weights are given by a weight vector ``w`` (of type ``WeightVec``).

.. function:: mean!(dst, x, w, dim)

  Compute weighted means along a certain dimension, and write results to a pre-allocated destination vector ``dst``.
