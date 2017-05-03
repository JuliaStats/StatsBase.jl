.. _weightvec:

Weight Vectors
================

In statistical applications, it is not uncommon to assign weights to samples. To facilitate the use of weight vectors, we introduce the abstract type ``AbstractWeights`` for the purpose of representing weight vectors.

Construction
--------------

A generic weight vector instance can be constructed using the ``Weights`` constructor or the ``weights`` function:

.. code-block:: julia

    w = Weights([1., 2., 3.])
    w = Weights([1., 2., 3.], 6.)

    w = weights([1., 2., 3.])

**Note:**

- The weight vector is a light-weight wrapper of the input vector. The input vector is NOT copied during construction.

- The weight vector maintains the sum of weights, which is computed upon construction. If the value of the sum is pre-computed, one can supply it as the second argument to the constructor and save the time of computing the sum again.


Methods
---------

Let ``w`` be an instance of ``AbstractWeights``:

.. function:: eltype(w)

    Get the type of weight values.

.. function:: length(w)

    Get the length of the weight vector.

.. function:: isempty(w)

    Test whether ``w`` is empty, *i.e.* ``length(w) == 0``.

.. function:: values(w)

    Get the vector of weight values.

.. function:: sum(w)

    Get the sum of weights.

    :note: The sum of weights is maintained by the weight vector, and thus this function can immediately return the value in ``O(1)`` (without computation).


Why we want an AbstractWeights type
------------------------------------

The ``AbstractWeights`` type is introduced as the standard way to pass weights, which has two advantages:

- A different type ``AbstractWeights`` distinguishes the role of the weight vector from other data vectors in the input arguments.
- Statistical functions that utilize weights often need the sum of weights for various purposes. The weight vector maintains the sum of weights, so that it needn't be computed repeatedly each time the sum of weights is needed.


Other AbstractWeights types and bias correction
-------------------------------------------------

When computing the weighted uncorrected (when `corrected=false`) sample variance, standard deviation or covariance a of :math:`\frac{1}{\sum{w}}` is used instead of :math:`\frac{1}{n}` where `n` is the number of observations.

Example:

:math:`s^2 = \frac{1}{\sum{w}} \sum_{i=1}^n {w_i\left({x_i - m}\right)^2 }`

vs

:math:`s^2 = \frac{1}{n} \sum_{i=1}^n {\left({x_i - m}\right)^2 }`

The unbiased estimate of the population variance, standard deviation or covariance is computed by replacing :math:`\frac{1}{\sum{w}}` with a factor dependent on the type of weights used:

- ``AnalyticWeights``: :math:`\frac{1}{\sum w - \sum {w^2} / \sum w}`
- ``FrequencyWeights``: :math:`\frac{1}{\sum{w} - 1}`
- ``ProbabilityWeights``: :math:`\frac{n}{(n - 1) \sum w}` where ``n`` equals `count(!iszero, w)`

These weights can be created with the appropriate constructor (i.e., ``AnalyticWeights(a)``, ``FrequencyWeights(a)``, ``ProbabilityWeights(a)``) or the utility functions ``aweights(a)``, ``fweights(a)`` and ``pweights(a)``.
