# Weight Vectors

In statistical applications, it is not uncommon to assign weights to samples. To facilitate the use of weight vectors, we introduce the abstract type `AbstractWeights` for the purpose of representing weight vectors, which has two advantages:

- A different type `AbstractWeights` distinguishes the role of the weight vector from other data vectors in the input arguments.
- Statistical functions that utilize weights often need the sum of weights for various purposes. The weight vector maintains the sum of weights, so that it needn't be computed repeatedly each time the sum of weights is needed.

!!! note
    - The weight vector is a light-weight wrapper of the input vector. The input vector is NOT copied during construction.
    - The weight vector maintains the sum of weights, which is computed upon construction. If the value of the sum is pre-computed, one can supply it as the second argument to the constructor and save the time of computing the sum again.


## Implementations

Several statistical weight types are provided which subtype `AbstractWeights`. The choice of weights impacts how bias is corrected in several methods. See the [`var`](@ref), [`std`](@ref) and [`cov`](@ref) docstrings for more details.

### `AnalyticWeights`

Analytic weights describe a non-random relative importance (usually between 0 and 1) for each observation. These weights may also be referred to as reliability weights, precision weights or inverse variance weights. These are typically used when the observations being weighted are aggregate values (e.g., averages) with differing variances.

```julia
w = AnalyticWeights([0.2, 0.1, 0.3])
w = aweights([0.2, 0.1, 0.3])
```

### `FrequencyWeights`

Frequency weights describe the number of times (or frequency) each observation was observed. These weights may also be referred to as case weights or repeat weights.

```julia
w = FrequencyWeights([2, 1, 3])
w = fweights([2, 1, 3])
```

### `ProbabilityWeights`

Probability weights represent the inverse of the sampling probability for each observation, providing a correction mechanism for under- or over-sampling certain population groups. These weights may also be referred to as sampling weights.

```julia
w = ProbabilityWeights([0.2, 0.1, 0.3])
w = pweights([0.2, 0.1, 0.3])
```

### `Weights`

The `Weights` type describes a generic weights vector which does not support all operations possible for `FrequencyWeights`, `AnalyticWeights` and `ProbabilityWeights`.

```julia
w = Weights([1., 2., 3.])
w = weights([1., 2., 3.])
```

## Methods

`AbstractWeights` implements the following methods:
```
eltype
length
isempty
values
sum
```

The following constructors are provided:
```@docs
AnalyticWeights
FrequencyWeights
ProbabilityWeights
Weights
aweights
fweights
pweights
weights
```