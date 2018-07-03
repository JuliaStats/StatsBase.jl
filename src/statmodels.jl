# Statistical Models

abstract type StatisticalModel end

"""
    coef(obj::StatisticalModel)

Return the coefficients of the model.
"""
coef(obj::StatisticalModel) = error("coef is not defined for $(typeof(obj)).")

"""
    coefnames(obj::StatisticalModel)

Return the names of the coefficients.
"""
coefnames(obj::StatisticalModel) = error("coefnames is not defined for $(typeof(obj)).")

"""
    coeftable(obj::StatisticalModel)

Return a table of class `CoefTable` with coefficients and related statistics.
"""
coeftable(obj::StatisticalModel) = error("coeftable is not defined for $(typeof(obj)).")

"""
    confint(obj::StatisticalModel)

Compute confidence intervals for coefficients.
"""
confint(obj::StatisticalModel) = error("coefint is not defined for $(typeof(obj)).")

"""
    deviance(obj::StatisticalModel)

Return the deviance of the model relative to a reference, which is usually when applicable
the saturated model. It is equal, *up to a constant*, to ``-2 \\log L``, with ``L``
the likelihood of the model.
"""
deviance(obj::StatisticalModel) = error("deviance is not defined for $(typeof(obj)).")

"""
    islinear(obj::StatisticalModel)

Indicate whether the model is linear.
"""
islinear(obj::StatisticalModel) = error("islinear is not defined for $(typeof(obj)).")

"""
    nulldeviance(obj::StatisticalModel)

Return the deviance of the null model, that is the one including only the intercept.
"""
nulldeviance(obj::StatisticalModel) = error("nulldeviance is not defined for $(typeof(obj)).")

"""
    loglikelihood(obj::StatisticalModel)

Return the log-likelihood of the model.
"""
loglikelihood(obj::StatisticalModel) = error("loglikelihood is not defined for $(typeof(obj)).")

"""
    loglikelihood(obj::StatisticalModel)

Return the log-likelihood of the null model corresponding to model `obj`.
This is usually the model containing only the intercept.
"""
nullloglikelihood(obj::StatisticalModel) = error("nullloglikelihood is not defined for $(typeof(obj)).")

"""
    score(obj::StatisticalModel)

Return the score of the statistical model. The score is the gradient of the
log-likelihood with respect to the coefficients.
"""
score(obj::StatisticalModel) = error("score is not defined for $(typeof(obj)).")

"""
    nobs(obj::StatisticalModel)

Return the number of independent observations on which the model was fitted. Be careful
when using this information, as the definition of an independent observation may vary
depending on the model, on the format used to pass the data, on the sampling plan
(if specified), etc.
"""
nobs(obj::StatisticalModel) = error("nobs is not defined for $(typeof(obj)).")

"""
    dof(obj::StatisticalModel)

Return the number of degrees of freedom consumed in the model, including
when applicable the intercept and the distribution's dispersion parameter.
"""
dof(obj::StatisticalModel) = error("dof is not defined for $(typeof(obj)).")

"""
    mss(obj::StatisticalModel)

Return the model sum of squares.
"""
mss(obj::StatisticalModel) = error("mss is not defined for $(typeof(obj)).")

"""
    rss(obj::StatisticalModel)

Return the residual sum of squares.
"""
rss(obj::StatisticalModel) = error("rss is not defined for $(typeof(obj)).")

"""
    informationmatrix(model::StatisticalModel; expected::Bool = true)

Return the information matrix. By default the Fisher information matrix is returned,
while the observed information matrix can be requested with `expected = false`.
"""
informationmatrix(model::StatisticalModel; expected::Bool = true) =
    error("informationmatrix is not defined for $(typeof(obj)).")

"""
    stderror(obj::StatisticalModel)

Return the standard errors for the coefficients of the model.
"""
stderror(obj::StatisticalModel) = sqrt.(diag(vcov(obj)))

"""
    vcov(obj::StatisticalModel)

Return the variance-covariance matrix for the coefficients of the model.
"""
vcov(obj::StatisticalModel) = error("vcov is not defined for $(typeof(obj)).")

"""
    weights(obj::StatisticalModel)

Return the weights used in the model.
"""
weights(obj::StatisticalModel) = error("weights is not defined for $(typeof(obj)).")

"""
    isfitted(obj::StatisticalModel)

Indicate whether the model has been fitted.
"""
isfitted(obj::StatisticalModel) = error("isfitted is not defined for $(typeof(obj)).")

"""
Fit a statistical model.
"""
fit(obj::StatisticalModel, args...) = error("fit is not defined for $(typeof(obj)).")

"""
Fit a statistical model in-place.
"""
fit!(obj::StatisticalModel, args...) = error("fit! is not defined for $(typeof(obj)).")

"""
    aic(obj::StatisticalModel)

Akaike's Information Criterion, defined as ``-2 \\log L + 2k``, with ``L`` the likelihood
of the model, and `k` its number of consumed degrees of freedom
(as returned by [`dof`](@ref)).
"""
aic(obj::StatisticalModel) = -2loglikelihood(obj) + 2dof(obj)

"""
    aicc(obj::StatisticalModel)

Corrected Akaike's Information Criterion for small sample sizes (Hurvich and Tsai 1989),
defined as ``-2 \\log L + 2k + 2k(k-1)/(n-k-1)``, with ``L`` the likelihood of the model,
``k`` its number of consumed degrees of freedom (as returned by [`dof`](@ref)),
and ``n`` the number of observations (as returned by [`nobs`](@ref)).
"""
function aicc(obj::StatisticalModel)
    k = dof(obj)
    n = nobs(obj)
    -2loglikelihood(obj) + 2k + 2k*(k+1)/(n-k-1)
end

"""
    bic(obj::StatisticalModel)

Bayesian Information Criterion, defined as ``-2 \\log L + k \\log n``, with ``L``
the likelihood of the model,  ``k`` its number of consumed degrees of freedom
(as returned by [`dof`](@ref)), and ``n`` the number of observations
(as returned by [`nobs`](@ref)).
"""
bic(obj::StatisticalModel) = -2loglikelihood(obj) + dof(obj)*log(nobs(obj))

"""
    r2(obj::StatisticalModel)
    r²(obj::StatisticalModel)

Coefficient of determination (R-squared).

For a linear model, the R² is defined as ``ESS/TSS``, with ``ESS`` the explained sum of squares
and ``TSS`` the total sum of squares.
"""
r2(obj::StatisticalModel) = mss(obj) / deviance(obj)

"""
    r2(obj::StatisticalModel, variant::Symbol)
    r²(obj::StatisticalModel, variant::Symbol)

Pseudo-coefficient of determination (pseudo R-squared).

For nonlinear models, one of several pseudo R² definitions must be chosen via `variant`.
Supported variants are:
- `:MacFadden` (a.k.a. likelihood ratio index), defined as ``1 - \\log L/\\log L0``.
- `:CoxSnell`, defined as ``1 - (L0/L)^{2/n}``
- `:Nagelkerke`, defined as ``(1 - (L0/L)^{2/n})/(1 - L0^{2/n})``, with ``n`` the number
of observations (as returned by [`nobs`](@ref)).

In the above formulas, ``L`` is the likelihood of the model, ``L0`` that of the null model
(the model including only the intercept). These two quantities are taken to be minus half
`deviance` of the corresponding models.
"""
function r2(obj::StatisticalModel, variant::Symbol)
    ll = -deviance(obj)/2
    ll0 = -nulldeviance(obj)/2

    if variant == :McFadden
        1 - ll/ll0
    elseif variant == :CoxSnell
        1 - exp(2/nobs(obj) * (ll0 - ll))
    elseif variant == :Nagelkerke
        (1 - exp(2/nobs(obj) * (ll0 - ll)))/(1 - exp(2/nobs(obj) * ll0))
    else
        error("variant must be one of :McFadden, :CoxSnell or :Nagelkerke")
    end
end

const r² = r2

"""
    adjr2(obj::StatisticalModel)
    adjr²(obj::StatisticalModel)

Adjusted coefficient of determination (adjusted R-squared).

For linear models, the adjusted R² is defined as ``1 - (1 - (1-R^2)(n-1)/(n-p))``, with ``R^2``
the coefficient of determination, ``n`` the number of observations, and ``p`` the number of
coefficients (including the intercept). This definition is generally known as the Wherry Formula I.
"""
adjr2(obj::StatisticalModel) = error("adjr2 is not defined for $(typeof(obj)).")

"""
    adjr2(obj::StatisticalModel, variant::Symbol)
    adjr²(obj::StatisticalModel, variant::Symbol)

Adjusted pseudo-coefficient of determination (adjusted pseudo R-squared).

For nonlinear models, one of the several pseudo R² definitions must be chosen via `variant`.
The only currently supported variant is `:MacFadden`, defined as ``1 - (\\log L - k)/\\log L0``.
In this formula, ``L`` is the likelihood of the model, ``L0`` that of the null model
(the model including only the intercept). These two quantities are taken to be minus half
`deviance` of the corresponding models. ``k`` is the number of consumed degrees of freedom
of the model (as returned by [`dof`](@ref)).
"""
function adjr2(obj::StatisticalModel, variant::Symbol)
    ll = -deviance(obj)/2
    ll0 = -nulldeviance(obj)/2
    k = dof(obj)

    if variant == :McFadden
        1 - (ll - k)/ll0
    else
        error(":McFadden is the only currently supported variant")
    end
end

const adjr² = adjr2

abstract type RegressionModel <: StatisticalModel end

"""
    fitted(obj::RegressionModel)

Return the fitted values of the model.
"""
fitted(obj::RegressionModel) = error("fitted is not defined for $(typeof(obj)).")

"""
    response(obj::RegressionModel)

Return the model response (a.k.a. the dependent variable).
"""
response(obj::RegressionModel) = error("response is not defined for $(typeof(obj)).")

"""
    meanresponse(obj::RegressionModel)

Return the mean of the response.
"""
meanresponse(obj::RegressionModel) = error("meanresponse is not defined for $(typeof(obj)).")

"""
    modelmatrix(obj::RegressionModel)

Return the model matrix (a.k.a. the design matrix).
"""
modelmatrix(obj::RegressionModel) = error("modelmatrix is not defined for $(typeof(obj)).")

"""
    leverage(obj::RegressionModel)

Return the diagonal of the projection matrix.
"""
leverage(obj::RegressionModel) = error("leverage is not defined for $(typeof(obj)).")

"""
    residuals(obj::RegressionModel)

Return the residuals of the model.
"""
residuals(obj::RegressionModel) = error("residuals is not defined for $(typeof(obj)).")

"""
    predict(obj::RegressionModel, [newX])

Form the predicted response of model `obj`. An object with new covariate values `newX` can be supplied,
which should have the same type and structure as that used to fit `obj`; e.g. for a GLM
it would generally be a `DataFrame` with the same variable names as the original predictors.
"""
function predict end

predict(obj::RegressionModel) = error("predict is not defined for $(typeof(obj)).")

"""
    predict!

In-place version of [`predict`](@ref).
"""
function predict! end

predict!(obj::RegressionModel) = error("predict! is not defined for $(typeof(obj)).")

"""
    dof_residual(obj::RegressionModel)

Return the residual degrees of freedom of the model.
"""
dof_residual(obj::RegressionModel) = error("dof_residual is not defined for $(typeof(obj)).")

"""
    params(obj)

Return all parameters of a model.
"""
params(obj) = error("params is not defined for $(typeof(obj))")
function params! end

## coefficient tables with specialized show method

## Nms are the coefficient names, corresponding to rows in the table
mutable struct CoefTable
    cols::Vector
    colnms::Vector
    rownms::Vector
    function CoefTable(cols::Vector,colnms::Vector,rownms::Vector)
        nc = length(cols)
        nrs = map(length,cols)
        nr = nrs[1]
        length(colnms) in [0,nc] || error("colnms should have length 0 or $nc")
        length(rownms) in [0,nr] || error("rownms should have length 0 or $nr")
        all(nrs .== nr) || error("Elements of cols should have equal lengths, but got $nrs")
        new(cols,colnms,rownms)
    end

    function CoefTable(mat::Matrix,colnms::Vector,rownms::Vector,pvalcol::Int=0)
        nc = size(mat,2)
        cols = Any[mat[:, i] for i in 1:nc]
        if pvalcol != 0                         # format the p-values column
            cols[pvalcol] = [PValue(cols[pvalcol][j])
                            for j in eachindex(cols[pvalcol])]
        end
        CoefTable(cols,colnms,rownms)
    end
end

mutable struct PValue
    v::Number
    function PValue(v::Number)
        0. <= v <= 1. || isnan(v) || error("p-values must be in [0.,1.]")
        new(v)
    end
end

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

function show(io::IO, ct::CoefTable)
    cols = ct.cols; rownms = ct.rownms; colnms = ct.colnms;
    nc = length(cols)
    nr = length(cols[1])
    if length(rownms) == 0
        rownms = [lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    rnwidth = max(4,maximum([length(nm) for nm in rownms]) + 1)
    rownms = [rpad(nm,rnwidth) for nm in rownms]
    widths = [length(cn)::Int for cn in colnms]
    str = String[isa(cols[j][i], AbstractString) ? cols[j][i] :
        sprint(show, cols[j][i], context=:compact=>true) for i in 1:nr, j in 1:nc]
    for j in 1:nc
        for i in 1:nr
            lij = length(str[i,j])
            if lij > widths[j]
                widths[j] = lij
            end
        end
    end
    widths .+= 1
    println(io," " ^ rnwidth *
            join([lpad(string(colnms[i]), widths[i]) for i = 1:nc], ""))
    for i = 1:nr
        print(io, rownms[i])
        for j in 1:nc
            print(io, lpad(str[i,j],widths[j]))
        end
        println(io)
    end
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
    function ConvergenceException{T}(iters, lastchange::T, tol::T) where T<:Real
        if tol > lastchange
            throw(ArgumentError("Change must be greater than tol."))
        else
            new(iters, lastchange, tol)
        end
    end
end

ConvergenceException(iters, lastchange::T=NaN, tol::T=NaN) where {T<:Real} =
    ConvergenceException{T}(iters, lastchange, tol)

function Base.showerror(io::IO, ce::ConvergenceException)
    print(io, "failure to converge after $(ce.iters) iterations.")
    if !isnan(ce.lastchange)
        print(io, " Last change ($(ce.lastchange)) was greater than tolerance ($(ce.tol)).")
    end
end
