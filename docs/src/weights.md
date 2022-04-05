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

### `UnitWeights`

Unit weights are a special case in which all observations are given a weight equal to `1`. Using such weights is equivalent to computing unweighted statistics.

This type can notably be used when implementing an algorithm so that a only a weighted variant has to be written. The unweighted variant is then obtained by passing a `UnitWeights` object. This is very efficient since no weights vector is actually allocated.

```julia
w = uweights(3)
w = uweights(Float64, 3)
```

### `Weights`

The `Weights` type describes a generic weights vector which does not support all operations possible for `FrequencyWeights`, `AnalyticWeights`, `ProbabilityWeights` and `UnitWeights`.

```julia
w = Weights([1., 2., 3.])
w = weights([1., 2., 3.])
```

### Exponential weights: `eweights`

Exponential weights are a common form of temporal weights which assign exponentially decreasing
weights to past observations.

If `t` is a vector of temporal indices then for each index `i` we compute the weight as:

``λ (1 - λ)^{1 - i}``

``λ`` is a smoothing factor or rate parameter such that ``0 < λ ≤ 1``.
As this value approaches 0, the resulting weights will be almost equal,
while values closer to 1 will put greater weight on the tail elements of the vector.

For example, the following call generates exponential weights for ten observations with ``λ = 0.3``.
```julia-repl
julia> eweights(1:10, 0.3)
10-element Weights{Float64,Float64,Array{Float64,1}}:
 0.3
 0.42857142857142855
 0.6122448979591837
 0.8746355685131197
 1.249479383590171
 1.7849705479859588
 2.549957925694227
 3.642797036706039
 5.203995766722913
 7.434279666747019
```

Simply passing the number of observations `n` is equivalent to passing in `1:n`.

```julia-repl
julia> eweights(10, 0.3)
10-element Weights{Float64,Float64,Array{Float64,1}}:
 0.3
 0.42857142857142855
 0.6122448979591837
 0.8746355685131197
 1.249479383590171
 1.7849705479859588
 2.549957925694227
 3.642797036706039
 5.203995766722913
 7.434279666747019
```

Finally, you can construct exponential weights from an arbitrary subset of timestamps within a larger range.

```julia-repl
julia> t
2019-01-01T01:00:00:2 hours:2019-01-01T05:00:00

julia> r
2019-01-01T01:00:00:1 hour:2019-01-02T01:00:00

julia> eweights(t, r, 0.3)
3-element Weights{Float64,Float64,Array{Float64,1}}:
 0.3
 0.6122448979591837
 1.249479383590171
```

NOTE: This is equivalent to `eweights(something.(indexin(t, r)), 0.3)`, which is saying that for each value in `t` return the corresponding index for that value in `r`.
Since `indexin` returns `nothing` if there is no corresponding value from `t` in `r` we use `something` to eliminate that possibility.

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
UnitWeights
Weights
aweights
fweights
pweights
eweights
uweights
weights(vs::AbstractArray{<:Real})
```
