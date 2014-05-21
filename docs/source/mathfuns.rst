Statistics-related Math functions
===================================

The package provides a set of math functions that are related to statistical computation:

- **xlogx(x)**(x): ``x * log(x)`` when ``x > 0``, or zero otherwise.

- **xlogy(x, y)**(x, y): ``x * log(y)`` when ``x > 0``, or zero otherwise.

- **logistic**(x): ``1 / (1 + exp(-x))``.

- **logit(x)**: ``log(x / (1 - x))``.

- **softplus(x)**: ``log(1 + exp(x))``.

- **invsoftplus(x)**: ``log(exp(x) - 1)``.

**Note:** all functions listed above have vectorized versions.

- **logsumexp(x, y)**: ``log(exp(x) + exp(y))``.

- **logsumexp(x)**: ``log(sum(exp(x)))`` when ``x`` is an array.

