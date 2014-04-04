Canonical Correlation
=====================

The package provides canonical correlation analysis.

- **canoncor** (X, Y)

    Compute the canonical correlation of ``X`` and ``Y``. Returns a tuple of form ``(A, B, r)`` where ``A`` and ``B`` are the canonical coefficients of ``X`` and ``Y`` and ``r`` are the canonical correlations. ``X`` and ``Y`` are scaled so that the canonical variables ``(X .- mean(X, 1))*A`` and ``(Y .- mean(Y, 1))*B`` have identity covariance. ``X`` and ``Y`` are assumed to be full rank.
