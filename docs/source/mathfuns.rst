Statistics-related Math functions
===================================

The package provides a set of math functions that are related to statistical computation:

.. function:: xlogx(x)

    ``x * log(x)`` when ``x > 0``, or ``0`` otherwise.

.. function:: xlogy(x, y)

    ``x * log(y)`` when ``x > 0``, or ``0`` otherwise.

.. function:: logistic(x)

    Logistic function: ``1 / (1 + exp(-x))``.

.. function:: logit(x) 

    Logit function: ``log(x / (1 - x))``.

.. function:: softplus(x)

    Softplus function: ``log(1 + exp(x))``.

.. function:: invsoftplus(x)

    Inverse softplus function: ``log(exp(x) - 1)``.

**Note:** all functions listed above have vectorized versions.

.. function:: logsumexp(x, y)

    Logarithm of sum of exponents: ``log(exp(x) + exp(y))``, computed in a numerically stable way.

.. function:: logsumexp(x)

    Logarithm of sum of exponents of all elements in ``x``, *i.e.* ``log(sum(exp(x)))``. Here, ``x`` is an array.

.. function:: softmax(x)

    Softmax function: ``exp(x) ./ sum(exp(x))`` for a given array ``x``.

.. function:: softmax!(r, x)

    Compute the softmax function over ``x`` and write the results to ``r``. Here, ``r`` should have the same size as ``x``.

.. function:: softmax!(x)

    Compute the softmax function over ``x`` and write the results inplace to ``x``.
