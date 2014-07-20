Computing Deviations
=====================

This package provides functions to compute various deviations between arrays in a variety of ways:

.. function:: counteq(a, b)

    Count the number of equal pairs of elements in ``a`` and ``b``, *i.e* ``countnz(a .== b)``.

.. function:: countne(a, b)

    Count the number of non-equal pairs of elements in ``a`` and ``b``, *i.e* ``countnz(a .!= b)``.

.. function:: sqL2dist(a, b)

    Squared L2 distance between ``a`` and ``b``, as :math:`\sum_{i=1}^n |a_i - b_i|^2`.

.. function:: L2dist(a, b)

    L2 distance between ``a`` and ``b``, *i.e* ``sqrt(sqL2dist(a, b))``. 

.. function:: L1dist(a, b)

    L1 distance between ``a`` and ``b``, as :math:`\sum_{i=1}^n |a_i - b_i|`.

.. function:: Linfdist(a, b)

    Linf distance between ``a`` and ``b``, as :math:`\mathrm{\mathop{max}}_{i=1:n} |a_i - b_i|`. 

.. function:: gkldiv(a, b)

    Generalized Kullback-Leibler divergence between two arrays ``a`` and ``b``, defined as
    :math:`\sum_{i=1}^n (a_i * \log(a_i/b_i) - a_i + b_i)`. 

    **Note:** When ``sum(a) == 1`` and ``sum(b) == 1``, it reduces to the KL-divergence in standard sense.

.. function:: meanad(a, b)

    Mean absolute deviation between ``a`` and ``b``, *i.e.* ``mean(abs(a - b))``.

.. function:: maxad(a, b)

    Maximum absolute deviation between ``a`` and ``b``, *i.e.* ``maximum(abs(a - b))``.

.. function:: msd(a, b)

    Mean squared deviation between ``a`` and ``b``, *i.e.* ``mean(abs2(a - b))``.

.. function:: rmsd(a, b[; normalize={true|false}])

    Root mean squared deviation, *i.e.* ``sqrt(msd(a, b))``.

    The keyword argument ``normalize`` is default to ``false``. If it is set to ``true``, the result is normalized by ``(maximum(a) - minimum(a)``.

.. function:: psnr(a, b, maxv)

    Peak signal-to-noise ratio, *i.e.* ``10 * log10(maxv^2 / msd(a, b))``.

**Note:** all these functions are implemented in a reasonably efficient way without creating any temporary arrays in the middle.

