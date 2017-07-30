# Contributing

Thanks for taking a look at StatsBase!

In general this guide summarizes and extends the [Julia Contributing Guidelines](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md#general-formatting-guidelines-for-julia-code-contributions) and [Julia Style Guide](https://docs.julialang.org/en/latest/manual/style-guide.html).

## Reporting Issues

* It's always helpful to start by reviewing existing issues and if necessary referencing back to previous discussions.
* Including minimal sample code (either to reproduce a bug or demonstrating desired behaviour) is greatly appreciated.
* For bugs, it's helpful to checkout the latest development version (`Pkg.checkout("StatsBase")`) and post the output from `versioninfo()`.

## Pull Requests (PRs)

* Try to keep changes in PRs focused on the specific issue at hand. This will reduce the amount of work for the reviewer and likely speed up the review process.
* PRs are expected to have passing tests and not reduce the overall coverage.
* Tests are run automatically and a coverage report will be posted in the PR discussion when tests pass.

### Checking Changes Locally

* Using the `git diff` and `git diff --stat` commands can be useful for comparing your feature branch against the upstream master before pushing.
* Tests can be run locally with `Pkg.test("StatsBase")`.
* [Coverage](https://github.com/JuliaCI/Coverage.jl) can be useful for checking if you've lowered the overall test coverage against [master](https://coveralls.io/r/JuliaStats/StatsBase.jl?branch=master).

## Style

### Summary

* Use 4 spaces per indentation level, no tabs.
* Try to adhere to the 92 character line length limit.
* Use upper camel case convention for [modules](https://docs.julialang.org/en/stable/manual/modules/) and [types](https://docs.julialang.org/en/stable/manual/types/).
* Use lower case for method names, preferably without underscores.
* Comments are a good way to explain the intentions of your code and can aid in Pull Request (PR) reviews.
* No whitespace at the end of a line (trailing whitespace).
* Avoid padding brackets with spaces. ex. `Int64(value)` preferred over `Int64( value )`.
* Include spaces after commas and around operators

### Methods Definitions

Short form method definitions when the method body only contains 1 statement:

```julia
# Yes
Base.stdm(v::RealArray, w::AbstractWeights, m::Real; corrected::Bool=true) =
    sqrt(varm(v, w, m, corrected=corrected))

# No
function Base.stdm(v::RealArray, w::AbstractWeights, m::Real; corrected::Bool=true)
    sqrt(varm(v, w, m, corrected=corrected))
end
```

For long method definitions (and calls) our preferences is to align arguments on new lines with the first argument and avoid adding an extra newline after the last argument:

```julia
# Yes
function Base.var(v::RealArray, w::AbstractWeights; mean=nothing,
                  corrected::Bool=true)
    ...
end

# No
function Base.var(
    v::RealArray, w::AbstractWeights; mean=nothing, corrected::Bool=true
)
    ...
end

# No
function Base.var(
    v::RealArray,
    w::AbstractWeights;
    mean=nothing,
    corrected::Bool=true
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
                  corrected::Bool=true)
    if mean == nothing
        varm(v, w, Base.mean(v, w); corrected=corrected)
    else
        varm(v, w, mean; corrected=corrected)
    end
end

# No
function Base.var(v::RealArray, w::AbstractWeights; mean=nothing,
                  corrected::Bool=true)
    if mean == nothing
        return varm(v, w, Base.mean(v, w); corrected=corrected)
    else
        return varm(v, w, mean; corrected=corrected)
    end
end
```

### Documentation

All exported functions must have a docstring associated with it.
For details on writing docstrings please review the
[documentation section](https://docs.julialang.org/en/stable/manual/documentation/)
of the julia manual.

Our hosted documentation is generated with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl),
to generate the html files locally install Documenter (e.g., `Pkg.add("Documenter")`) and run
`include(make.jl)` from the docs directory.
You can then inspect the generated html packages in your
browser to confirm that your docstrings are displayed correctly.
This is particularly useful if your docstrings include LaTeX equations.
