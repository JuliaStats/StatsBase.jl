Mean Functions
===============

The package provides functions to compute means of different kinds.

.. py:function:: geomean(x)

  Compute the geometric mean of ``x``.

.. py:function:: harmmean(x)

  Compute the harmonic mean of ``x``.

.. py:function:: trimmean(x, p)

  Compute the trimmed mean of ``x``, with fraction ``p`` of elements ignored.

  **Examples:**

  .. code-block:: julia

    # computes the mean using 80% of the numbers in ``x``, 
    # ignoring 10% of the largest and 10% of the smallest values.

    trimmean(x, 0.2)

.. py:function:: wmean(x, w[; wsum=NaN])  

  Compute weighted mean of ``x`` using weights in ``w``. Here, ``w`` is an ordinary vector. 

  One can use ``wsum`` to provide a pre-computed weight. By default, ``wsum`` is set to ``NaN``, which indicates to compute the sum of weights in the function.

.. py:function:: mean(x, w)

  The ``mean`` function is also extended to accept a weight vector of type ``WeightVec`` (see :ref:`weightvec`). 

  **Examples:**

  .. code-block:: julia

    w = rand(n)
    x = mean(x, weights(w))

