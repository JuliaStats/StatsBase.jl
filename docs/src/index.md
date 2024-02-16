# Getting Started

```@meta
CurrentModule = StatsBase
DocTestSetup = quote
    using Statistics
    using Random
end
```

*StatsBase.jl* is a Julia package that provides basic support for statistics. Particularly, it implements a variety of statistics-related functions, such as scalar statistics, high-order moment computation, counting, ranking, covariances, sampling, and empirical density estimation.

## Installation

To install StatsBase through the Julia REPL, you can type `] add StatsBase` or:
```julia
using Pkg
Pkg.add("StatsBase")
```

To load the package, use the command:
```
using StatsBase
```

## Available Features

```@contents
Pages = ["weights.md", "scalarstats.md", "robust.md", "deviation.md", "cov.md", "counts.md", "ranking.md", "sampling.md", "empirical.md", "signalcorr.md", "misc.md", "statmodels.md", "transformations.md"]
Depth = 2
```


