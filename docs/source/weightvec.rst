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


Variance bias correction
-------------------------------------------------

When computing the weighted uncorrected (when ``corrected=false``) sample variance, standard deviation or covariance :math:`\frac{1}{\sum{w}}` is used instead of :math:`\frac{1}{n}` (where ``n`` is the number of observations).

Example:

:math:`s^2 = \frac{1}{\sum{w}} \sum_{i=1}^n {w_i\left({x_i - m}\right)^2 }`

vs

:math:`s^2 = \frac{1}{n} \sum_{i=1}^n {\left({x_i - m}\right)^2 }`

However, unbiased estimates (when ``corrected=true``) are dependent on the types of weights used. All weights presented here are a subtype of ``AbstractWeights``.

Weights
~~~~~~~

The `Weights` type describes a generic weights vector which does not support all operations possible for ``FrequencyWeights``, ``AnalyticWeights`` and ``ProbabilityWeights``.

- ``corrected=true``: ``ArgumentError``
- ``corrected=false``: :math:`\frac{1}{\sum{w}}`

AnalyticWeights
~~~~~~~~~~~~~~~~

Analytic weights describe a non-random relative importance (usually between 0 and 1) for each observation. These weights may also be referred to as reliability weights, precision weights or inverse variance weights. These are typically used when the observations being weighted are aggregate values (e.g., averages) with differing variances.

- ``corrected=true``: :math:`\frac{1}{\sum w - \sum {w^2} / \sum w}`
- ``corrected=false``: :math:`\frac{1}{\sum{w}}`

.. code-block:: julia

    w = AnalyticWeights([0.2, 0.1, 0.3])

    w = aweights([0.2, 0.1, 0.3])


FrequencyWeights
~~~~~~~~~~~~~~~~~

Frequency weights describe the number of times (or frequency) each observation
was observed. These weights may also be referred to as case weights or repeat weights.

- ``corrected=true``: :math:`\frac{1}{\sum{w} - 1}`
- ``corrected=false``: :math:`\frac{1}{\sum{w}}`

.. code-block:: julia

    w = FrequencyWeights([2, 1, 3])

    w = fweights([2, 1, 3])


ProbabilityWeights
~~~~~~~~~~~~~~~~~~~

Probability weights represent the inverse of the sampling probability for each observation, providing a correction mechanism for under- or over-sampling certain population groups. These weights may also be referred to as sampling weights.

- ``corrected=true``: :math:`\frac{n}{(n - 1) \sum w}` where ``n`` equals ``count(!iszero, w)``
- ``corrected=false``: :math:`\frac{1}{\sum{w}}`

.. code-block:: julia

    w = ProbabilityWeights([0.2, 0.1, 0.3])

    w = pweights([0.2, 0.1, 0.3])
