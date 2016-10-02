Miscellaneous Functions
========================

.. function:: rle(x)

    Run-length encoding of ``x``. It returns ``(vals, lens)``, a sequence of values and their corresponding chunk length. [Wikipedia](http://en.wikipedia.org/wiki/Run-length_encoding).

    **Examples:**

    .. code-block:: julia

        julia> rle([1,1,1,2,2,3,3,3,3,2,2,2])
        ([1,2,3,2],[3,2,4,3])


.. function:: inverse_rle(vals, lens)

    Inversed run-length encoding. It takes the results of ``rle`` and reconstructs the original sequence. 

.. function:: levelsmap(x)

    Construct a dictionary that maps each of the ``n`` distinct values in ``x`` to a number between ``1`` and ``n``.

.. function:: indexmap(x)

    Construct a dictionary that maps each distinct value in ``x`` to its first index.

.. function:: indicatormat(x, k[; sparse=false])  

    Construct a boolean matrix ``r`` of size ``(k, length(x))`` such that ``r[x[i], i] = true`` and all other elements are set to ``false``.

    The keyword argument ``sparse`` controls whether to construct a sparse matrix. By default, it is false. 

    **Examples:**

    .. code-block:: julia

        julia> indicatormat([1 2 2], 2)
        2x3 Array{Bool,2}:
        true   false  false
        false  true   true

.. function:: indicatormat(x, c[; sparse=false])

    Construct a boolean matrix ``r`` of size ``(length(c), length(x))``. Let ``ci`` be the index of ``x[i]`` in ``c``, then ``r[ci, i] = true`` and all other elements are zero. 

    The keyword argument ``sparse`` controls whether to construct a sparse matrix. By default, it is set to ``false``. 

.. function:: indicatormat(x[; sparse=false])

    Equivalent to ``indicatormap(x, sort(unique(x)); sparse=...)``. 

