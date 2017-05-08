# Contributing

Thanks for taking a look at StatsBase!

In general this guide summarizes and extends the [Julia Contributing Guidelines](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md#general-formatting-guidelines-for-julia-code-contributions) and [Julia Style Guide](https://docs.julialang.org/en/latest/manual/style-guide.html).

## Reporting Issues

* It's always helpful to start by reviewing existing issues is and if necessary referencing back to previous discussions.
* Including minimal sample code (either to reproduce a bug or demonstrating desired behaviour) is greatly appreciated.
* For bugs, it's helpful to checkout the latest development version (`Pkg.checkout("StatsBase")`) and post the output from `versioninfo()`.

## Pull Requests (PR)s

* Try to keep changes in PRs constrained if possible. This will reduce the amount of work for the reviewer and likely speed up the review process.
* Using `git diff upstream/master` and `git diff --stat upstream/master` (assuming you've run `git remote add upstream https://github.com/JuliaStats/StatsBase.jl`) can be useful commands for comparing your branch against the upstream master before pushing.
* Running tests locally with `Pkg.test("StatsBase")` before committing can be a helpful way of reducing noise.
* [Coverage](https://github.com/JuliaCI/Coverage.jl) can be useful for checking if you've lowered the overall test coverage against [master](https://coveralls.io/r/JuliaStats/StatsBase.jl?branch=master)

## Style

### Summary

* Use 4 spaces per indentation level, no tabs.
* Try to adhere to the 92 character line length limit.
* Use upper camel case convention for [modules](http://julia.readthedocs.org/en/latest/manual/modules/) and [types](http://julia.readthedocs.org/en/latest/manual/types/).
* Use lower case for method names, preferably without underscores.
* Comments are a good way to explain the intentions of your code and can aid in Pull Request (PR) reviews.
* No whitespace at the end of a line (trailing whitespace).
* Avoid padding brackets with spaces. ex. `Int64(value)` preferred over `Int64( value )`.
* Include spaces after commas and around operators

### Methods Definitions

Short form method definitions when the method body only contains 1 statement:

```julia
# Yes
Base.stdm(v::RealArray, w::AbstractWeights, m::Real; corrected::DepBool=nothing) =
    sqrt(varm(v, w, m, corrected=depcheck(:stdm, corrected)))

# No
function Base.stdm(v::RealArray, w::AbstractWeights, m::Real; corrected::DepBool=nothing)
    sqrt(varm(v, w, m, corrected=depcheck(:stdm, corrected)))
end
```

For long method definitions (and calls) our preferences is to align arguments on new lines with the first argument and avoid adding an extra newline after the last argument:

```julia
# Yes
function Base.var(v::RealArray, w::AbstractWeights; mean=nothing,
                  corrected::DepBool=nothing)
    ...
end

# No
function Base.var(
    v::RealArray, w::AbstractWeights; mean=nothing, corrected::DepBool=nothing
)
    ...
end

# No
function Base.var(
    v::RealArray,
    w::AbstractWeights;
    mean=nothing,
    corrected::DepBool=nothing
)
    ...
end

# Yes
scale!(_wsum_centralize!(R, abs2, A, values(w), M, dim, true),
       varcorrection(w, corrected))

# No
scale!(
    _wsum_centralize!(R, abs2, A, values(w), M, dim, true),
    varcorrection(w, corrected)
)
```

### `return`

For simple methods, we try to avoid using the `return` keyword:

```julia
# Yes
function Base.var(v::RealArray, w::AbstractWeights; mean=nothing,
                  corrected::DepBool=nothing)
    corrected = depcheck(:var, corrected)

    if mean == nothing
        varm(v, w, Base.mean(v, w); corrected=corrected)
    else
        varm(v, w, mean; corrected=corrected)
    end
end

# No
function Base.var(v::RealArray, w::AbstractWeights; mean=nothing,
                  corrected::DepBool=nothing)
    corrected = depcheck(:var, corrected)

    if mean == nothing
        return varm(v, w, Base.mean(v, w); corrected=corrected)
    else
        return varm(v, w, mean; corrected=corrected)
    end
end
```

However, for more complicated methods it's helpful to be more explicit:

```julia
# Yes
function counteq(a::AbstractArray, b::AbstractArray)
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent lengths."))
    c = 0
    for i = 1:n
        @inbounds if a[i] == b[i]
            c += 1
        end
    end
    return c
end

# No
function counteq(a::AbstractArray, b::AbstractArray)
    n = length(a)
    length(b) == n || throw(DimensionMismatch("Inconsistent lengths."))
    c = 0
    for i = 1:n
        @inbounds if a[i] == b[i]
            c += 1
        end
    end
    c   # The `c` return is less obvious here
end
```

### Documentation

All exported functions should have a docstring associated with it.
A docstring should start with the method signature(s) followed by a description
of what it does. If a function has many arguments it may be advisable to add a
`# Arguments` section with a description of each argument.
Extra whitespace and newlines in docstrings should generally be avoided unless
it is required to separate logical components.

Examples)
```julia
"""
    weights(vs)

Construct a `Weights` vector from array `vs`.
See the documentation for [`Weights`](@ref) for more details.
"""
weights(vs::RealVector) = Weights(vs)
weights(vs::RealArray) = Weights(vec(vs))
```
```julia
"""
    var(x, w::AbstractWeights, [dim]; mean=nothing, corrected=false)

Compute the variance of a real-valued array `x`, optionally over a dimension `dim`.
Observations in `x` are weighted using weight vector `w`.
The uncorrected (when `corrected=false`) sample variance is defined as:
```math
\\frac{1}{\\sum{w}} \\sum_{i=1}^n {w_i\\left({x_i - μ}\\right)^2 }
```
where ``n`` is the length of the input and ``μ`` is the mean.
The unbiased estimate (when `corrected=true`) of the population variance is computed by
replacing ``\\frac{1}{\\sum{w}}`` with a factor dependent on the type of weights used:
* `AnalyticWeights`: ``\\frac{1}{\\sum w - \\sum {w^2} / \\sum w}``
* `FrequencyWeights`: ``\\frac{1}{\\sum{w} - 1}``
* `ProbabilityWeights`: ``\\frac{n}{(n - 1) \\sum w}`` where ``n`` equals `count(!iszero, w)`
* `Weights`: `ArgumentError` (bias correction not supported)
"""
function Base.var(v::RealArray, w::AbstractWeights; mean=nothing,
                  corrected::DepBool=nothing)
    ...
end
```
