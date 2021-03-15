# Statistical Models

abstract type StatisticalModel end

"""
    coef(model::StatisticalModel)

Return the coefficients of the model.
"""
coef(model::StatisticalModel) = error("coef is not defined for $(typeof(model)).")

"""
    coefnames(model::StatisticalModel)

Return the names of the coefficients.
"""
coefnames(model::StatisticalModel) = error("coefnames is not defined for $(typeof(model)).")

"""
    coeftable(model::StatisticalModel; level::Real=0.95)

Return a table with coefficients and related statistics of the model.
`level` determines the level for confidence intervals (by default, 95%).

The returned `CoefTable` object implements the
[Tables.jl](https://github.com/JuliaData/Tables.jl/) interface, and can be
converted e.g. to a `DataFrame` via `using DataFrames; DataFrame(coeftable(model))`.
"""
coeftable(model::StatisticalModel) = error("coeftable is not defined for $(typeof(model)).")

"""
    confint(model::StatisticalModel; level::Real=0.95)

Compute confidence intervals for coefficients, with confidence level `level` (by default 95%).
"""
confint(model::StatisticalModel) = error("confint is not defined for $(typeof(model)).")

"""
    deviance(model::StatisticalModel)

Return the deviance of the model relative to a reference, which is usually when applicable
the saturated model. It is equal, *up to a constant*, to ``-2 \\log L``, with ``L``
the likelihood of the model.
"""
deviance(model::StatisticalModel) = error("deviance is not defined for $(typeof(model)).")

"""
    islinear(model::StatisticalModel)

Indicate whether the model is linear.
"""
islinear(model::StatisticalModel) = error("islinear is not defined for $(typeof(model)).")

"""
    nulldeviance(model::StatisticalModel)

Return the deviance of the null model, that is the one including only the intercept.
"""
nulldeviance(model::StatisticalModel) =
    error("nulldeviance is not defined for $(typeof(model)).")

"""
    loglikelihood(model::StatisticalModel)

Return the log-likelihood of the model.
"""
loglikelihood(model::StatisticalModel) =
    error("loglikelihood is not defined for $(typeof(model)).")

"""
    loglikelihood(model::StatisticalModel)

Return the log-likelihood of the null model corresponding to `model`.
This is usually the model containing only the intercept.
"""
nullloglikelihood(model::StatisticalModel) =
    error("nullloglikelihood is not defined for $(typeof(model)).")

"""
    score(model::StatisticalModel)

Return the score of the model, that is the gradient of the
log-likelihood with respect to the coefficients.
"""
score(model::StatisticalModel) = error("score is not defined for $(typeof(model)).")

"""
    nobs(model::StatisticalModel)

Return the number of independent observations on which the model was fitted. Be careful
when using this information, as the definition of an independent observation may vary
depending on the model, on the format used to pass the data, on the sampling plan
(if specified), etc.
"""
nobs(model::StatisticalModel) = error("nobs is not defined for $(typeof(model)).")

"""
    dof(model::StatisticalModel)

Return the number of degrees of freedom consumed in the model, including
when applicable the intercept and the distribution's dispersion parameter.
"""
dof(model::StatisticalModel) = error("dof is not defined for $(typeof(model)).")

"""
    mss(model::StatisticalModel)

Return the model sum of squares.
"""
mss(model::StatisticalModel) = error("mss is not defined for $(typeof(model)).")

"""
    rss(model::StatisticalModel)

Return the residual sum of squares of the model.
"""
rss(model::StatisticalModel) = error("rss is not defined for $(typeof(model)).")

"""
    informationmatrix(model::StatisticalModel; expected::Bool = true)

Return the information matrix of the model. By default the Fisher information matrix
is returned, while the observed information matrix can be requested with `expected = false`.
"""
informationmatrix(model::StatisticalModel; expected::Bool = true) =
    error("informationmatrix is not defined for $(typeof(model)).")

"""
    stderror(model::StatisticalModel)

Return the standard errors for the coefficients of the model.
"""
stderror(model::StatisticalModel) = sqrt.(diag(vcov(model)))

"""
    vcov(model::StatisticalModel)

Return the variance-covariance matrix for the coefficients of the model.
"""
vcov(model::StatisticalModel) = error("vcov is not defined for $(typeof(model)).")

"""
    weights(model::StatisticalModel)

Return the weights used in the model.
"""
weights(model::StatisticalModel) = error("weights is not defined for $(typeof(model)).")

"""
    isfitted(model::StatisticalModel)

Indicate whether the model has been fitted.
"""
isfitted(model::StatisticalModel) = error("isfitted is not defined for $(typeof(model)).")

"""
Fit a statistical model.
"""
fit(model::StatisticalModel, args...) = error("fit is not defined for $(typeof(model)).")

"""
Fit a statistical model in-place.
"""
fit!(model::StatisticalModel, args...) = error("fit! is not defined for $(typeof(model)).")

"""
    aic(model::StatisticalModel)

Akaike's Information Criterion, defined as ``-2 \\log L + 2k``, with ``L`` the likelihood
of the model, and `k` its number of consumed degrees of freedom
(as returned by [`dof`](@ref)).
"""
aic(model::StatisticalModel) = -2loglikelihood(model) + 2dof(model)

"""
    aicc(model::StatisticalModel)

Corrected Akaike's Information Criterion for small sample sizes (Hurvich and Tsai 1989),
defined as ``-2 \\log L + 2k + 2k(k-1)/(n-k-1)``, with ``L`` the likelihood of the model,
``k`` its number of consumed degrees of freedom (as returned by [`dof`](@ref)),
and ``n`` the number of observations (as returned by [`nobs`](@ref)).
"""
function aicc(model::StatisticalModel)
    k = dof(model)
    n = nobs(model)
    -2loglikelihood(model) + 2k + 2k*(k+1)/(n-k-1)
end

"""
    bic(model::StatisticalModel)

Bayesian Information Criterion, defined as ``-2 \\log L + k \\log n``, with ``L``
the likelihood of the model,  ``k`` its number of consumed degrees of freedom
(as returned by [`dof`](@ref)), and ``n`` the number of observations
(as returned by [`nobs`](@ref)).
"""
bic(model::StatisticalModel) = -2loglikelihood(model) + dof(model)*log(nobs(model))

"""
    r2(model::StatisticalModel)
    r²(model::StatisticalModel)

Coefficient of determination (R-squared).

For a linear model, the R² is defined as ``ESS/TSS``, with ``ESS`` the explained sum of squares
and ``TSS`` the total sum of squares.
"""
function r2(model::StatisticalModel)
    Base.depwarn("The default r² method for linear models is deprecated. " *
                 "Packages should define their own methods.", :r2)

    mss(model) / deviance(model)
end

"""
    r2(model::StatisticalModel, variant::Symbol)
    r²(model::StatisticalModel, variant::Symbol)

Pseudo-coefficient of determination (pseudo R-squared).

For nonlinear models, one of several pseudo R² definitions must be chosen via `variant`.
Supported variants are:
- `:MacFadden` (a.k.a. likelihood ratio index), defined as ``1 - \\log (L)/\\log (L_0)``;
- `:CoxSnell`, defined as ``1 - (L_0/L)^{2/n}``;
- `:Nagelkerke`, defined as ``(1 - (L_0/L)^{2/n})/(1 - L_0^{2/n})``.
- `:devianceratio`, defined as ``1 - D/D_0``.

In the above formulas, ``L`` is the likelihood of the model,
``L_0`` is the likelihood of the null model (the model with only an intercept),
``D`` is the deviance of the model (from the saturated model),
``D_0`` is the deviance of the null model,
``n`` is the number of observations (given by [`nobs`](@ref)).

The Cox-Snell and the deviance ratio variants both match the classical definition of R²
for linear models.
"""
function r2(model::StatisticalModel, variant::Symbol)
    loglikbased = (:McFadden, :CoxSnell, :Nagelkerke)
    if variant in loglikbased
        ll = loglikelihood(model)
        ll0 = nullloglikelihood(model)
        if variant == :McFadden
            1 - ll/ll0
        elseif variant == :CoxSnell
            1 - exp(2 * (ll0 - ll) / nobs(model))
        elseif variant == :Nagelkerke
            (1 - exp(2 * (ll0 - ll) / nobs(model))) / (1 - exp(2 * ll0 / nobs(model)))
        end
    elseif variant == :devianceratio
        dev  = deviance(model)
        dev0 = nulldeviance(model)
        1 - dev/dev0
    else
        error("variant must be one of $(join(loglikbased, ", ")) or :devianceratio")
    end
end

const r² = r2

"""
    adjr2(model::StatisticalModel)
    adjr²(model::StatisticalModel)

Adjusted coefficient of determination (adjusted R-squared).

For linear models, the adjusted R² is defined as ``1 - (1 - (1-R^2)(n-1)/(n-p))``, with ``R^2``
the coefficient of determination, ``n`` the number of observations, and ``p`` the number of
coefficients (including the intercept). This definition is generally known as the Wherry Formula I.
"""
adjr2(model::StatisticalModel) = error("adjr2 is not defined for $(typeof(model)).")

"""
    adjr2(model::StatisticalModel, variant::Symbol)
    adjr²(model::StatisticalModel, variant::Symbol)

Adjusted pseudo-coefficient of determination (adjusted pseudo R-squared).

For nonlinear models, one of the several pseudo R² definitions must be chosen via `variant`.
The only currently supported variants are `:MacFadden`, defined as ``1 - (\\log (L) - k)/\\log (L0)`` and
`:devianceratio`, defined as ``1 - (D/(n-k))/(D_0/(n-1))``.
In these formulas, ``L`` is the likelihood of the model, ``L0`` that of the null model
(the model including only the intercept), ``D`` is the deviance of the model,
``D_0`` is the deviance of the null model, ``n`` is the number of observations (given by [`nobs`](@ref)) and
``k`` is the number of consumed degrees of freedom of the model (as returned by [`dof`](@ref)).
"""
function adjr2(model::StatisticalModel, variant::Symbol)
    k = dof(model)
    if variant == :McFadden
        ll = loglikelihood(model)
        ll0 = nullloglikelihood(model)
        1 - (ll - k)/ll0
    elseif variant == :devianceratio
        n = nobs(model)
        dev  = deviance(model)
        dev0 = nulldeviance(model)
        1 - (dev*(n-1))/(dev0*(n-k))
    else
        error("variant must be one of :McFadden or :devianceratio")
    end
end

const adjr² = adjr2

abstract type RegressionModel <: StatisticalModel end

"""
    fitted(model::RegressionModel)

Return the fitted values of the model.
"""
fitted(model::RegressionModel) = error("fitted is not defined for $(typeof(model)).")

"""
    response(model::RegressionModel)

Return the model response (a.k.a. the dependent variable).
"""
response(model::RegressionModel) = error("response is not defined for $(typeof(model)).")

"""
    responsename(model::RegressionModel)

Return the name of the model response (a.k.a. the dependent variable).
"""
responsename(model::RegressionModel) = error("responsename is not defined for $(typeof(model)).")

"""
    meanresponse(model::RegressionModel)

Return the mean of the response.
"""
meanresponse(model::RegressionModel) = error("meanresponse is not defined for $(typeof(model)).")

"""
    modelmatrix(model::RegressionModel)

Return the model matrix (a.k.a. the design matrix).
"""
modelmatrix(model::RegressionModel) = error("modelmatrix is not defined for $(typeof(model)).")

"""
    crossmodelmatrix(model::RegressionModel)

Return `X'X` where `X` is the model matrix of `model`.
This function will return a pre-computed matrix stored in `model` if possible.
"""
crossmodelmatrix(model::RegressionModel) = (x = modelmatrix(model); Symmetric(x' * x))

"""
    leverage(model::RegressionModel)

Return the diagonal of the projection matrix of the model.
"""
leverage(model::RegressionModel) = error("leverage is not defined for $(typeof(model)).")

"""
    residuals(model::RegressionModel)

Return the residuals of the model.
"""
residuals(model::RegressionModel) = error("residuals is not defined for $(typeof(model)).")

"""
    predict(model::RegressionModel, [newX])

Form the predicted response of `model`. An object with new covariate values `newX` can be supplied,
which should have the same type and structure as that used to fit `model`; e.g. for a GLM
it would generally be a `DataFrame` with the same variable names as the original predictors.
"""
function predict end

predict(model::RegressionModel) = error("predict is not defined for $(typeof(model)).")

"""
    predict!

In-place version of [`predict`](@ref).
"""
function predict! end

predict!(model::RegressionModel) = error("predict! is not defined for $(typeof(model)).")

"""
    dof_residual(model::RegressionModel)

Return the residual degrees of freedom of the model.
"""
dof_residual(model::RegressionModel) = error("dof_residual is not defined for $(typeof(model)).")

"""
    params(model)

Return all parameters of a model.
"""
params(model) = error("params is not defined for $(typeof(model))")
function params! end

## coefficient tables with specialized show method

mutable struct CoefTable
    cols::Vector
    colnms::Vector
    rownms::Vector
    pvalcol::Int
    teststatcol::Int
    function CoefTable(cols::Vector,colnms::Vector,rownms::Vector,
                       pvalcol::Int=0,teststatcol::Int=0)
        nc = length(cols)
        nrs = map(length,cols)
        nr = nrs[1]
        length(colnms) in [0,nc] || throw(ArgumentError("colnms should have length 0 or $nc"))
        length(rownms) in [0,nr] || throw(ArgumentError("rownms should have length 0 or $nr"))
        all(nrs .== nr) || throw(ArgumentError("Elements of cols should have equal lengths, but got $nrs"))
        pvalcol in 0:nc || throw(ArgumentError("pvalcol should be between 0 and $nc"))
        teststatcol in 0:nc || throw(ArgumentError("teststatcol should be between 0 and $nc"))
        new(cols,colnms,rownms,pvalcol,teststatcol)
    end

    function CoefTable(mat::Matrix,colnms::Vector,rownms::Vector,
                       pvalcol::Int=0,teststatcol::Int=0)
        nc = size(mat,2)
        cols = Any[mat[:, i] for i in 1:nc]
        CoefTable(cols,colnms,rownms,pvalcol,teststatcol)
    end
end

Base.length(ct::CoefTable) = length(ct.cols[1])
function Base.eltype(ct::CoefTable)
    names = isempty(ct.rownms) ?
        tuple(Symbol.(ct.colnms)...) :
        tuple(Symbol("Name"), Symbol.(ct.colnms)...)
    types = isempty(ct.rownms) ?
        Tuple{eltype.(ct.cols)...} :
        Tuple{eltype(ct.rownms), eltype.(ct.cols)...}
    NamedTuple{names, types}
end

function Base.iterate(ct::CoefTable, i::Integer=1)
    if i in 1:length(ct)
        cols = getindex.(ct.cols, Ref(i))
        nt = isempty(ct.rownms) ?
            eltype(ct)(tuple(cols...)) :
            eltype(ct)(tuple(ct.rownms[i], cols...))
        (nt, i+1)
    else
        nothing
    end
end

"""
Show a p-value using 6 characters, either using the standard 0.XXXX
representation or as <Xe-YY.
"""
struct PValue <: Real
    v::Real
    function PValue(v::Real)
        0 <= v <= 1 || isnan(v) || error("p-values must be in [0; 1]")
        new(v)
    end
end
PValue(p::PValue) = p

function show(io::IO, pv::PValue)
    v = pv.v
    if isnan(v)
        @printf(io,"%d", v)
    elseif v >= 1e-4
        @printf(io,"%.4f", v)
    else
        @printf(io,"<1e%2.2d", ceil(Integer, max(nextfloat(log10(v)), -99)))
    end
end

"""Show a test statistic using 2 decimal digits"""
struct TestStat <: Real
    v::Real
end

show(io::IO, x::TestStat) = @printf(io, "%.2f", x.v)
TestStat(x::TestStat) = x

float(x::Union{TestStat, PValue}) = float(x.v)

for op in [:(==), :<, :≤, :>, :≥, :(isless), :(isequal)] # isless and < to place nice with NaN
    @eval begin
        Base.$op(x::Union{TestStat, PValue}, y::Real) = $op(x.v, y)
        Base.$op(y::Real, x::Union{TestStat, PValue}) = $op(y, x.v)
        Base.$op(x1::Union{TestStat, PValue}, x2::Union{TestStat, PValue}) = $op(x1.v, x2.v)
    end
end

# necessary to avoid a method ambiguity with isless(::TestStat, NaN)
Base.isless(x::Union{TestStat, PValue}, y::AbstractFloat) = isless(x.v, y)
Base.isless(y::AbstractFloat, x::Union{TestStat, PValue},) = isless(y, x.v)
Base.isequal(y::AbstractFloat, x::Union{TestStat, PValue}) = isequal(y, x.v)
Base.isequal(x::Union{TestStat, PValue}, y::AbstractFloat) = isequal(x.v, y)


Base.isapprox(x::Union{TestStat, PValue}, y::Real; kwargs...) = isapprox(x.v, y; kwargs...)
Base.isapprox(y::Real, x::Union{TestStat, PValue}; kwargs...) = isapprox(y, x.v; kwargs...)
Base.isapprox(x1::Union{TestStat, PValue}, x2::Union{TestStat, PValue}; kwargs...) = isapprox(x1.v, x2.v; kwargs...)


"""Wrap a string so that show omits quotes"""
struct NoQuote
    s::String
end

show(io::IO, n::NoQuote) = print(io, n.s)

function show(io::IO, ct::CoefTable)
    cols = ct.cols; rownms = ct.rownms; colnms = ct.colnms;
    nc = length(cols)
    nr = length(cols[1])
    if length(rownms) == 0
        rownms = [lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    mat = [j == 1 ? NoQuote(rownms[i]) :
           j-1 == ct.pvalcol ? NoQuote(sprint(show, PValue(cols[j-1][i]))) :
           j-1 in ct.teststatcol ? TestStat(cols[j-1][i]) :
           cols[j-1][i] isa AbstractString ? NoQuote(cols[j-1][i]) : cols[j-1][i]
           for i in 1:nr, j in 1:nc+1]
    # Code inspired by print_matrix in Base
    io = IOContext(io, :compact=>true, :limit=>false)
    A = Base.alignment(io, mat, 1:size(mat, 1), 1:size(mat, 2),
                       typemax(Int), typemax(Int), 3)
    nmswidths = pushfirst!(length.(colnms), 0)
    A = [nmswidths[i] > sum(A[i]) ? (A[i][1]+nmswidths[i]-sum(A[i]), A[i][2]) : A[i]
         for i in 1:length(A)]
    totwidth = sum(sum.(A)) + 2 * (length(A) - 1)
    println(io, repeat('─', totwidth))
    print(io, repeat(' ', sum(A[1])))
    for j in 1:length(colnms)
        print(io, "  ", lpad(colnms[j], sum(A[j+1])))
    end
    println(io, '\n', repeat('─', totwidth))
    for i in 1:size(mat, 1)
        Base.print_matrix_row(io, mat, A, i, 1:size(mat, 2), "  ")
        i != size(mat, 1) && println(io)
    end
    print(io, '\n', repeat('─', totwidth))
    nothing
end

show(mime::MIME"text/markdown", ct::CoefTable) = show(Base.stdout, mime, ct)

function show(io::IO, ::MIME"text/markdown", ct::CoefTable)
    cols = ct.cols; rownms = ct.rownms; colnms = ct.colnms;
    nc = length(cols)
    nr = length(cols[1])
    if length(rownms) == 0
        rownms = [lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    mat = [j == 1 ? NoQuote(rownms[i]) :
           j-1 == ct.pvalcol ? PValue(cols[j-1][i]) :
           j-1 in ct.teststatcol ? TestStat(cols[j-1][i]) :
           cols[j-1][i] isa AbstractString ? NoQuote(cols[j-1][i]) : cols[j-1][i]
           for i in 1:nr, j in 1:nc+1]
    # Code inspired by print_matrix in Base
    io = IOContext(io, :compact=>true, :limit=>false)
    A = Base.alignment(io, mat, 1:size(mat, 1), 1:size(mat, 2),
                       typemax(Int), typemax(Int), 3)
    nmswidths = pushfirst!(length.(colnms), 0)
    A = [nmswidths[i] > sum(A[i]) ? (A[i][1]+nmswidths[i]-sum(A[i]), A[i][2]) : A[i]
         for i in 1:length(A)]
    totwidth = sum(sum.(A)) + 2 * (length(A) - 1)

    # not using Markdown stdlib here because that won't give us nice decimal
    # alignment (even if that is lost when rendering to HTML, it's still nice
    # when looking at the markdown itself)

    print(io, '|', ' '^(sum(A[1])+1))
    for j in 1:length(colnms)
        print(io, " | ", lpad(colnms[j], sum(A[j+1])))
    end

    println(io, " |")
    print(io, '|', rpad(':', sum(A[1])+2, '-'))
    for j in 1:length(colnms)
        _pad = j-1 in [ct.teststatcol; ct.pvalcol] ? rpad : lpad
        print(io, '|', _pad(':', sum(A[j+1])+2, '-'))
    end
    println(io, '|')

    for i in 1:size(mat, 1)
        print(io, "| ")
        Base.print_matrix_row(io, mat, A, i, 1:size(mat, 2), " | ")
        print(io, " |")
        i != size(mat, 1) && println(io)
    end

    nothing
end

"""
    ConvergenceException(iters::Int, lastchange::Real=NaN, tol::Real=NaN)

The fitting procedure failed to converge in `iters` number of iterations,
i.e. the `lastchange` between the cost of the final and penultimate iteration was greater than
specified tolerance `tol`.
"""
struct ConvergenceException{T<:Real} <: Exception
    iters::Int
    lastchange::T
    tol::T
    msg::String
    function ConvergenceException{T}(iters, lastchange::T, tol::T, msg::String) where T<:Real
        if tol > lastchange
            throw(ArgumentError("Change must be greater than tol."))
        else
            new(iters, lastchange, tol, msg)
        end
    end
end

ConvergenceException(iters, lastchange::T=NaN, tol::T=NaN,
                     msg::AbstractString="") where {T<:Real} =
    ConvergenceException{T}(iters, lastchange, tol, String(msg))

function Base.showerror(io::IO, ce::ConvergenceException)
    print(io, "failure to converge after $(ce.iters) iterations.")
    if !isnan(ce.lastchange)
        print(io, " Last change ($(ce.lastchange)) was greater than tolerance ($(ce.tol)).")
    end
    if !isempty(ce.msg)
        print(io, ' ', ce.msg)
    end
end
