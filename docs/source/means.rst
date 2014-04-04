Mean Functions
===============

The package provides functions to compute means of different kinds.

- **geomean** (x)

  Compute the geometric mean of ``x``.

- **harmmean** (x)

  Compute the harmonic mean of ``x``.

- **trimmean** (x, p)

  Compute the trimmed mean of ``x``, with fraction ``p`` of elements ignored.

  **Examples:**

  .. code-block:: julia

    # computes the mean using 80% of the numbers in ``x``, 
    # ignoring 10% of the largest and 10% of the smallest values.

    trimmean(x, 0.2)

- **mean** (x, w)

  The ``mean`` function is also extended to accept a weight vector of type ``WeightVec`` (see :ref:`weightvec`) to compute weighted mean. 

  **Examples:**

  .. code-block:: julia

    w = rand(n)
    x = mean(x, weights(w))

