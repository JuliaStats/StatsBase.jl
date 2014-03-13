Sampling from Population
--------------------------

The package provides functions for sampling from a given population (with or without replacement).

.. py:function:: sample(a)

    Randomly draw an element from an array ``a``.

.. py:function:: sample(a, n[; replace=true, ordered=false])  

    Randomly draw ``n`` elements from ``a``. 

    **Keyword arguments**

    - ``replace``: indicates whether to have replacement (default = ``true``).
    - ``ordered``: indicates whether to arrange the samples in ascending order (default = ``false``).

.. py:function:: sample!(a, x[; replace=true, ordered=false])

    Draw ``length(x)`` elements from ``a`` and write them to a pre-allocated array ``x``.

.. py:function:: sample(wv) 

    Draw an integer in ``1:length(w)`` with probabilities proportional to the weights given in ``wv``. 

    Here, ``wv`` should be a weight vector of type ``WeightVec`` (see :ref:`weightvec`).

.. py:function:: sample(a, wv)

    Draw an element from ``a`` with probabilities proportional to the corresponding weights given in ``wv``.

.. py:function:: sample(a, wv, n[; replace=true, ordered=false])

    Draw ``n`` elements from ``a`` with probabilities proportional to the corresponding weights given in ``wv``.

    **Keyword arguments**

    - ``replace``: indicates whether to have replacement (default = ``true``).
    - ``ordered``: indicates whether to arrange the samples in ascending order (default = ``false``).    

.. py:function:: sample!(a, wv, x[; replace=true, ordered=false])

    Weighted sampling with results written to a pre-allocated array ``x``.

