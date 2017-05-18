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

.. function:: cov(X, w[; vardim=..., mean=..., corrected=...])

  Compute the weighted covariance matrix. Similar to ``var`` and ``std`` the biased covariance matrix (``corrected=false``) is computed by multiplying ``scattermat(X, w)`` by :math:`\frac{1}{\sum{w}}` to normalize.
  However, the unbiased covariance matrix (``corrected=true``) is dependent on the type of weights used:

    * ``AnalyticWeights``: :math:`\frac{1}{\sum w - \sum {w^2} / \sum w}`
    * ``FrequencyWeights``: :math:`\frac{1}{\sum{w} - 1}`
    * ``ProbabilityWeights``: :math:`\frac{n}{(n - 1) \sum w}` where ``n`` equals ``count(!iszero, w)``
    * ``Weights``: ``ArgumentError`` (bias correction not supported)

.. function:: cor(X, w[, vardim])

    Compute the Pearson correlation matrix of ``X`` along the dimension ``vardim`` with a weighting ``w``.

.. function:: mean_and_cov(x[, wv][; vardim=..., corrected=...])

  Jointly compute the mean and covariance matrix as a tuple.
  A weighting vector ``wv`` can be specified. ``vardim`` that designates whether the variables are columns in the matrix (``1``) or rows (``2``).
  Finally, bias correction is applied to the covariance calculation if ``corrected=true``.
  See ``cov`` documentation for more details.

.. function:: cov2cor(C, s)

    Compute the correlation matrix from the covariance matrix ``C`` and a vector of standard deviations ``s``. Use ``Base.cov2cor!`` for an in-place version.

.. function:: cor2cov(C, s)

    Compute the covariance matrix from the correlation matrix ``C`` and a vector of standard deviations ``s``. Use ``StatsBase.cor2cov!`` for an in-place version.
