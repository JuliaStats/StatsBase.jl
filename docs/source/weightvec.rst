.. _weightvec:

Weight Vectors
================

In statistical applications, it is not uncommon to assign weights to samples. To facilitate the use of weight vectors, we introduce the abstract type ``AbstractWeights`` for the purpose of representing weight vectors.

**Note:**

- The weight vector is a light-weight wrapper of the input vector. The input vector is NOT copied during construction.

- The weight vector maintains the sum of weights, which is computed upon construction. If the value of the sum is pre-computed, one can supply it as the second argument to the constructor and save the time of computing the sum again.


Implementations
---------------

Several statistical weight types are provided which subtype ``AbstractWeights``. The choice of weights impacts how bias is corrected in several methods. See the ``var``, ``std`` and ``cov`` docstrings for more details.

``Weights``
~~~~~~~~~~~~

The ``Weights`` type describes a generic weights vector which does not support all operations possible for ``FrequencyWeights``, ``AnalyticWeights`` and ``ProbabilityWeights``.

.. code-block:: julia

    w = Weights([1., 2., 3.])
    w = Weights([1., 2., 3.,], 6.)

    w = weights([1., 2., 3.])


``AnalyticWeights``
~~~~~~~~~~~~~~~~~~~~

Analytic weights describe a non-random relative importance (usually between 0 and 1) for each observation. These weights may also be referred to as reliability weights, precision weights or inverse variance weights. These are typically used when the observations being weighted are aggregate values (e.g., averages) with differing variances.

.. code-block:: julia

    w = AnalyticWeights([0.2, 0.1, 0.3])
    w = aweights([0.2, 0.1, 0.3])


``FrequencyWeights``
~~~~~~~~~~~~~~~~~~~~~

Frequency weights describe the number of times (or frequency) each observation was observed. These weights may also be referred to as case weights or repeat weights.

.. code-block:: julia

    w = FrequencyWeights([2, 1, 3])
    w = fweights([2, 1, 3])


``ProbabilityWeights``
~~~~~~~~~~~~~~~~~~~~~~

Probability weights represent the inverse of the sampling probability for each observation, providing a correction mechanism for under- or over-sampling certain population groups. These weights may also be referred to as sampling weights.

.. code-block:: julia

    w = ProbabilityWeights([0.2, 0.1, 0.3])
    w = pweights([0.2, 0.1, 0.3])


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


Why we want an ``AbstractWeights`` type
----------------------------------------

The ``AbstractWeights`` type is introduced as the standard way to pass weights, which has two advantages:

- A different type ``AbstractWeights`` distinguishes the role of the weight vector from other data vectors in the input arguments.
- Statistical functions that utilize weights often need the sum of weights for various purposes. The weight vector maintains the sum of weights, so that it needn't be computed repeatedly each time the sum of weights is needed.
