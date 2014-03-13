Counting Functions
=====================

The package provides functions to count the occurences of distinct values.

Counting over an Integer Range
--------------------------------

.. py:function:: counts(x, a:b[, wv])

    Count the number of times (or total weights if a weight vector ``wv`` is given) values in ``a:b`` appear in array ``x``. Here, the optional argument ``wv`` should be a weight vector of type ``WeightVec`` (see :ref:`weightvec`).

    This function returns a vector ``r`` of length ``n``, with ``n = length(a:b) = b-a+1``. In particular, we have

        ``r[i] == countnz(x .== a+(i-1))``

    **Convenient forms:**

    =========================  ======================
      convenient interface      equivalents
    -------------------------  ----------------------
     ``counts(x, k)``           ``counts(x, 1:k)``
     ``counts(x, k, w)``        ``counts(x, 1:k, w)``
    =========================  ======================

    **Examples:**

    .. code-block:: julia

        julia> counts([1, 2, 2, 2, 3, 3], 1:3)
        3-element Array{Int64,1}:
         1
         3
         2

    **Notes:** When the function encounters an value (or pair) that does not fall in the specified range, they simply ignore it (without raising an error or warning).

.. py:function:: counts(x, y, (a:b, c:d)[, wv])

    Count the number of times (or total weights if a weight vector ``wv`` is given) pairs of values in ``a:b`` and ``c:d`` that appear in arrays ``x`` and ``y``.

    This function returns a matrix ``r`` of size ``(nx, ny)``, where ``nx = b-a+1`` and ``ny=d-c+1``. In particular, we have 

        ``r[i,j] == countnz((x .== a+(i-1)) & (y .== b+(j-1)))``

    **Convenient forms:**

    ==============================  ===================================
      convenient interface            equivalents
    ------------------------------  -----------------------------------
    ``counts(x, y, a:b)``           ``counts(x, y, (a:b, a:b))``
    ``counts(x, y, a:b, w)``        ``counts(x, y, (a:b, a:b), w)``
    ``counts(x, y, (nx, ny))``      ``counts(x, y, (1:nx, 1:ny))``
    ``counts(x, y, (nx, ny), w)``   ``counts(x, y, (1:nx, 1:ny), w)``
    ``counts(x, y, n)``             ``counts(x, y, (1:n, 1:n)``
    ``counts(x, y, n, w)``          ``counts(x, y, (1:n, 1:n), w)``
    ==============================  ===================================

    **Examples:**

    .. code-block:: julia

        julia> counts([1, 2, 1, 1, 2], [1, 2, 3, 3, 2], (1:2, 1:3))  
        2x3 Array{Int64,2}:
         1  0  2
         0  2  0
        
    This function is useful in various applications, *e.g.* computing confusion matrix in evaluating the performance of a classifier.

.. py:function:: proportions(x, a:b[, wv])  

    Compute the proportions of values in ``a:b`` with respect to ``x``. Equivalent to ``counts(x, a:b) / length(x)``. 

.. py:function:: proportions(x, y, (a:b, c:d)[, wv])

    Equivalent to ``counts(x, y, (a:b, c:d)) / length(x)``.

    **Note:** all convenient forms for the function ``counts`` are also available to the function ``proportions``.  

.. py:function:: addcounts!(r, x, a:b[, wv])

    Adds the counts of values in ``x`` to an accumulating array ``r``.

.. py:function:: addcounts!(r, x, y, (a:b, c:d)[, wv])

    Adds the counts of pairs in ``x`` and ``y`` to an accumulating matrix ``r``.  


Counting over arbitrary distinct values
---------------------------------------- 

.. py:function:: countmap(x[, wv])

    Return a dictionary that maps distinct values in ``x`` to their counts (or total weights).

.. py:function:: proportionmap(x[, wv])

    Return a dictionary that maps distinct values in ``x`` to their proportions. 

.. py:function:: addcounts!(dict, x[, wv])

    Add counts based on ``x`` to a count map. New entries will be added if new values come up.


