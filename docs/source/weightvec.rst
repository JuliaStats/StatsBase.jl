.. _weightvec:

Weight Vectors
================

In statistical applications, it is not uncommon to assign weights to samples. To facilitate the use of weight vectors, we introduce a type ``WeightVec`` for the purpose of representing weight vectors.

Construction
--------------

A weight vector instance can be constructed using the ``WeightVec`` constructor or the ``weights`` function:

.. code-block:: julia

    w = WeightVec([1., 2., 3.])
    w = WeightVec([1., 2., 3.], 6.)
    
    w = weights([1., 2., 3.])

**Note:** 

- The weight vector is a light-weight wrapper of the input vector. The input vector is NOT copied during construction.

- The weight vector maintains the sum of weights, which is computed upon construction. If the value of the sum is pre-computed, one can supply it as the second argument to the constructor and save the time of computing the sum again.


Methods
---------

Let ``w`` be an instance of ``WeightVec``:

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


Why we want a WeightVec type
-----------------------------

The ``WeightVector`` type is introduced as the standard way to pass weights, which has two advantages:

- A different type ``WeightVec`` distinguishes the role of the weight vector from other data vectors in the input arguments.
- Statistical functions that utilize weights often need the sum of weights for various purposes. The weight vector maintains the sum of weights, so that it needn't be computed repeatedly each time the sum of weights is needed.

