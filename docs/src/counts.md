# Counting Functions

The package provides functions to count the occurrences of distinct values.

## Counting over an Integer Range

```@docs
counts
proportions
addcounts!(r::AbstractArray, x::StatsBase.IntegerArray, levels::StatsBase.IntUnitRange)
```

## Counting over arbitrary distinct values

```@docs
countmap
proportionmap
addcounts!(cm::Dict, x::Any)
```
