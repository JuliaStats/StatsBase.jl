Add hsm_mode: half-sample mode for continuous data
Summary

This PR adds hsm_mode(), an implementation of the half-sample mode (HSM), a robust estimator of the mode for continuous distributions. It is introduced as a separate function (not an overload of mode()) to preserve existing behavior while providing a statistically meaningful alternative for continuous data.

This addresses and closes issue #957.

Motivation

StatsBase.mode() is frequency-based and works well for discrete data. For continuous distributions, however, samples are usually unique, which makes frequency counts unstable and highly variable in practice.

Issue #957 documents this behavior, particularly for heavy-tailed distributions, where mode() can show extreme variance. This PR provides an estimator designed specifically for continuous data.

Approach

hsm_mode() implements the standard half-sample method described in the literature:

Non-finite values (NaN, Inf) are filtered

The data are sorted

The algorithm repeatedly selects the contiguous half-sample with the smallest width

Once ≤ 2 points remain, the midpoint of the final interval is returned

The midpoint may not be a sample value, but provides a stable estimate of the location of highest density.

Time complexity is dominated by sorting (O(n log n)); space complexity is O(n).
After sorting, the contraction loop operates on SubArray views to avoid allocations.

API Design

The estimator is exposed as a new function:

hsm_mode(x::AbstractVector)


It is not added as an overload of mode() in order to:

avoid changing existing semantics

clearly distinguish frequency-based and density-based estimation

let users choose the appropriate method explicitly

The return type is an AbstractFloat, promoted from the input element type (e.g. integers → Float64, Float32 → Float32).

Testing and Documentation

Tests cover basic correctness, edge cases, robustness to outliers, handling of non-finite values, and type behavior. All tests pass.

The docstring explains intended use cases, compares with mode(), documents complexity, provides examples, and cites the relevant literature.

References

Robertson & Cryer (1974), JASA

Bickel & Fruehwirth (2006), CSDA

Notes

This PR is intentionally small and focused. Extensions such as weighted HSM or support for missing values are left for future work.

Feedback on naming or API placement is welcome.