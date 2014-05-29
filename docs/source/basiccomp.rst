Basic Computation
=====================

The package provides a set of functions to make computation (especially those related to statistics) convenient:


Inplace Arithmetics
---------------------

- **negate!** (x): Negate each element in ``x`` inplace.

- **add!** (y, x): Add ``x`` to ``y`` inplace, in a broadcasting manner. 

- **subtract!** (y, x): Subtract ``x`` from ``y`` inplace, in a broadcasting manner.

- **addscale!** (y, x, c): Add ``x * c`` to ``y`` inplace. Here, ``x`` and ``y`` should be arrays of the same size, and ``c`` be a scalar number. 

**Note:** No temporary arrays will be created within ``addscale!``. ``x * c`` won't be formed explicitly. 

- **abs!** (x): Apply ``abs`` inplace to ``x``.

- **abs2!** (x): Apply ``abs2`` inplace to ``x``.

- **sqrt!** (x): Apply ``sqrt`` inplace to ``x``.

- **exp!** (x): Apply ``exp`` inplace to ``x``.

- **log!** (x): Apply ``log`` inplace to ``x``.


Additional Statistics Functions
--------------------------------

- **sumabs** (x):  Sum of absolute values.

- **sumabs2** (x): Sum of squared absolute values.

- **maxabs** (x):  Maximum of absolute values.

- **sumabsdiff** (x, y): Sum of absolute differences. ``x`` and ``y`` can be either arrays or scalars.

- **sumabs2diff** (x, y): Squared sum of absolute differences. ``x`` and ``y`` can be either arrays or scalars.

- **maxabsdiff** (x, y): Maximum of absolute differences. ``x`` and ``y`` can be either arrays or scalars.


More Functions
-----------------

- **xlogx** (x): ``x * log(x)`` when ``x > 0``, or zero otherwise.

- **xlogy** (x, y): ``x * log(y)`` when ``x > 0``, or zero otherwise.

- **logistic** (x): ``1 / (1 + exp(-x))``.

- **logit** (x): ``log(x / (1 - x))``.

- **softplus** (x): ``log(1 + exp(x))``.

- **invsoftplus** (x): ``log(exp(x) - 1)``.

**Note:** all functions listed above have vectorized versions.

- **logsumexp** (x, y): ``log(exp(x) + exp(y))``.

- **logsumexp** (x): ``log(sum(exp(x)))`` when ``x`` is an array.

- **softmax** (x): ``exp(x) ./ sum(exp(x))`` for a given array ``x``.

- **softmax!** (r, x): write the results of ``softmax`` to ``r``.

- **softmax!** (x): write the results of ``softmax`` inplace to ``x``.
