Sampling from Population
=========================

Sampling API
--------------

The package provides functions for sampling from a given population (with or without replacement).

.. function:: sample([rng], a)

    Randomly draw an element from an array ``a``.

    Optionally specify a random number generator ``rng`` as the first argument (defaults to ``Base.GLOBAL_RNG``).

.. function:: sample([rng], a, n[; replace=true, ordered=false])  

    Randomly draw ``n`` elements from ``a``. 

    Optionally specify a random number generator ``rng`` as the first argument (defaults to ``Base.GLOBAL_RNG``).

    **Keyword arguments**

    - ``replace``: indicates whether to have replacement (default = ``true``).
    - ``ordered``: indicates whether to arrange the samples in ascending order (default = ``false``).

.. function:: sample!([rng], a, x[; replace=true, ordered=false])

    Draw ``length(x)`` elements from ``a`` and write them to a pre-allocated array ``x``.

    Optionally specify a random number generator ``rng`` as the first argument (defaults to ``Base.GLOBAL_RNG``).

.. function:: sample([rng], wv) 

    Draw an integer in ``1:length(wv)`` with probabilities proportional to the weights given in ``wv``. 

    Here, ``wv`` should be a weight vector of type ``WeightVec`` (see :ref:`weightvec`).

    Optionally specify a random number generator ``rng`` as the first argument (defaults to ``Base.GLOBAL_RNG``).

.. function:: sample([rng], a, wv)

    Draw an element from ``a`` with probabilities proportional to the corresponding weights given in ``wv``.

    Optionally specify a random number generator ``rng`` as the first argument (defaults to ``Base.GLOBAL_RNG``).

.. function:: sample([rng], a, wv, n[; replace=true, ordered=false])

    Draw ``n`` elements from ``a`` with probabilities proportional to the corresponding weights given in ``wv``.

    Optionally specify a random number generator ``rng`` as the first argument (defaults to ``Base.GLOBAL_RNG``).

    **Keyword arguments**

    - ``replace``: indicates whether to have replacement (default = ``true``).
    - ``ordered``: indicates whether to arrange the samples in ascending order (default = ``false``).    

.. function:: sample!([rng], a, wv, x[; replace=true, ordered=false])

    Weighted sampling with results written to a pre-allocated array ``x``.

    Optionally specify a random number generator ``rng`` as the first argument (defaults to ``Base.GLOBAL_RNG``).


Algorithms
-----------

Internally, this package implements multiple algorithms, and the ``sample`` (and ``sample!``) methods integrate them into a poly-algorithm, which chooses a specific algorithm based on inputs.

Note that the choices made in ``sample`` are decided based on extensive benchmarking (see ``perf/sampling.jl`` and ``perf/wsampling.jl``). It performs reasonably fast for most cases. That being said, if you know that a certain algorithm is particularly suitable for your context, directly calling an internal algorithm function might be slightly more efficient.

Here are a list of algorithms implemented in the package. The functions below are not exported (one can still import them from ``StatsBase`` using ``import`` though).

**Notations:**

- ``a``: source array representing the population
- ``x``: the destination array
- ``wv``: the weight vector (of type ``WeightVec``), for weighted sampling
- ``n``: the length of ``a``
- ``k``: the length of ``x``. For sampling without replacement, ``k`` must not exceed ``n``.
- ``rng``: optional random number generator (defaults to ``Base.GLOBAL_RNG``)

All following functions write results to ``x`` (pre-allocated) and return ``x``.


**Sampling Algorithms (Non-Weighted):**

.. function:: direct_sample!([rng], a, x)

    *Direct sampling.*

    For each ``j`` in ``1:k``, randomly pick ``i`` from ``1:n``, and set ``x[j] = a[i]``.

    This algorithm consumes ``k`` random numbers.

.. function:: samplepair([rng], a)

    Pick two elements at distinct positions from ``a``, and return them as a pair.

    This algorithm consumes exactly two random numbers.

.. function:: knuths_sample!([rng], a, x)

    *Knuth's Algorithm S* for random sampling without replacement.

    Reference: D. Knuth. *The Art of Computer Programming*. Vol 2, 3.4.2, p.142.

    This algorithm consumes ``n`` random numbers. It requires no additional memory space. Suitable for the case where memory is tight.

.. function:: fisher_yates_sample!([rng], a, x)

    *Fisher-Yates shuffling* (with early termination). 

    Pseudo-code ::

        create an array of index inds = [1:n]

        for i = 1:k
            swap inds[i] with a random one in inds[i:n]
            set x[i] = a[inds[i]]
        end
    

    This algorithm consumes ``k`` random numbers. It uses an integer array of length ``n`` internally to maintain the shuffled indices. It is considerably faster than Knuth's algorithm especially when ``n`` is greater than ``k``.

.. function:: self_avoid_sample!([rng], a, x)

    Use a set to maintain the index that has been sampled. Each time draw a new index, if the index has already been sampled, redraw until it draws an unsampled one. 

    This algorithm consumes about (or slightly more than) ``k`` random numbers, and requires ``O(k)`` memory to store the set of sampled indices. Very fast when ``n >> k``. 

    However, if ``k`` is large and approaches ``n``, the rejection rate would increase drastically, resulting in poorer performance.

.. function:: seqsample_a!([rng], a, x)

    *Algorithm A* described in the following paper (page 714).

    Jeffrey Scott Vitter. *Faster Methods for Random Sampling*. Communications of the ACM, 27 (7), July 1984.

    This algorithm consumes ``O(n)`` random numbers. The outputs are ordered.

.. function:: seqsample_c!([rng], a, x)

    *Algorithm C* described in the following paper (page 714).

    Jeffrey Scott Vitter. *Faster Methods for Random Sampling*. Communications of the ACM, 27 (7), July 1984.

    This algorithm consumes ``O(k^2)`` random numbers. The outputs are ordered.


**Weighted Sampling Algorithms:**

.. function:: direct_sample!([rng], a, wv, x)

    *Direct sampling.*

    Draw each sample by scanning the weight vector. 

    This algorithm: (1) consumes ``k`` random numbers; (2) has time complexity ``O(n k)``, as scanning the weight vector each time takes ``O(n)``; and (3) requires no additional memory space.

.. function:: alias_sample!([rng], a, wv, x)

    *Alias method.*

    Build a alias table, and sample therefrom.

    Reference: Walker, A. J. *An Efficient Method for Generating Discrete Random Variables with General Distributions.* ACM Transactions on Mathematical Software 3 (3): 253, 1977.

    This algorithm takes ``O(n log n)`` time for building the alias table, and then ``O(1)`` to draw each sample. It consumes ``2 k`` random numbers.

.. function:: naive_wsample_norep!([rng], a, wv, x)

    Naive implementation of weighted sampling without replacement.

    It makes a copy of the weight vector at initialization, and sets the weight to zero when the corresponding sample is picked.

    This algorithm consumes ``O(k)`` random numbers, and has overall time complexity ``O(n k)``. 

