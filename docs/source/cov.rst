Scatter Matrix and Covariance
===============================

This package implements functions for computing scatter matrix, as well as weighted covariance matrix.

.. function:: scattermat(X[; vardim=..., mean=...])

    Compute *scatter matrix* for the variables contained in ``X``.

    A `scatter matrix <http://en.wikipedia.org/wiki/Scatter_matrix>`_ can be considered as a unnormalized version of the covariance matrix.

    By default, it considers each column as a variable (*i.e* each row as an observation), and subtract the mean from each vector. One may change this default behavior by setting the keyword arguments.

    **Keyword arguments**

    - ``vardim``: the dimension of variable, which can take either ``1`` or ``2``.

        - ``vardim == 1``: each column is a variable, each row is an observation (default).
        - ``vardim == 2``: each column is a observation, each row is a variable.

    - ``mean``: pre-computed mean vector. Default value of ``mean`` is set to ``nothing``, which indicates that the function would compute the mean internally. One can also set ``mean`` to ``0``, which indicates that the input ``X`` has already been centralized. Otherwise, the supplied mean will be subtracted from ``X``.

.. function:: scatter(X, wv[; vardim=..., mean=...])

    Weighted scatter matrix. The weights are given by a weight vector ``wv`` of type ``AbstractWeights`` (see :ref:`weightvec`).

.. function:: cov(X, wv[; vardim=..., mean=..., corrected=...])

    Weighted covariance matrix.

    **Note:** By default, the covariance is normalized by the sum of weights, that is, ``cov(X, wv)`` is equal to ``scatter(X, wv) / sum(wv)``. However, if ``corrected`` is set to ``true`` then the appropriate bias correction is used for that `wv`.

.. function:: mean_and_cov(x[, wv][; vardim=..., corrected=...])

  Jointly compute the mean and covariance of ``x``.
