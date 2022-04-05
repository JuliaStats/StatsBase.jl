using StatsBase
using StatsBase: PValue, TestStat
using Test, Random, StatsAPI, LinearAlgebra

v1 = [1.45666, -23.14, 1.56734e-13]
v2 = ["Good", "Great", "Bad"]
v3 = [1, 56, 2]
v4 = [-12.56, 0.1326, 2.68e-16]
v5 = [0.12, 0.3467, 1.345e-16]
ct = CoefTable(Any[v1, v2, v3, v4, v5],
               ["Estimate", "Comments", "df", "t", "p"],
               ["x1", "x2", "x3"], 5, 4)
ct_noname = CoefTable(Any[v1, v2, v3, v4, v5],
                      ["Estimate", "Comments", "df", "t", "p"],
                      [], 5, 4)

@test sprint(show, ct) == """
───────────────────────────────────────────────
         Estimate  Comments  df       t       p
───────────────────────────────────────────────
x1    1.45666         Good    1  -12.56  0.1200
x2  -23.14            Great  56    0.13  0.3467
x3    1.56734e-13     Bad     2    0.00  <1e-15
───────────────────────────────────────────────"""

@test sprint(show, ct_noname) == """
────────────────────────────────────────────────
          Estimate  Comments  df       t       p
────────────────────────────────────────────────
[1]    1.45666         Good    1  -12.56  0.1200
[2]  -23.14            Great  56    0.13  0.3467
[3]    1.56734e-13     Bad     2    0.00  <1e-15
────────────────────────────────────────────────"""

@test sprint(show, MIME"text/markdown"(), ct) == """
|    |      Estimate | Comments | df |      t |      p |
|:---|--------------:|---------:|---:|-------:|:-------|
| x1 |   1.45666     |    Good  |  1 | -12.56 | 0.1200 |
| x2 | -23.14        |    Great | 56 |   0.13 | 0.3467 |
| x3 |   1.56734e-13 |    Bad   |  2 |   0.00 | <1e-15 |"""

@test sprint(show, MIME"text/markdown"(), ct_noname) == """
|     |      Estimate | Comments | df |      t |      p |
|:----|--------------:|---------:|---:|-------:|:-------|
| [1] |   1.45666     |    Good  |  1 | -12.56 | 0.1200 |
| [2] | -23.14        |    Great | 56 |   0.13 | 0.3467 |
| [3] |   1.56734e-13 |    Bad   |  2 |   0.00 | <1e-15 |"""

@test length(ct) === 3
@test eltype(ct) ==
    NamedTuple{(:Name, :Estimate, :Comments, :df, :t, :p),
               Tuple{String,Float64,String,Int,Float64,Float64}}
@test collect(ct) == [
    (Name = "x1", Estimate = 1.45666, Comments = "Good", df = 1, t = -12.56, p = 0.12)
    (Name = "x2", Estimate = -23.14, Comments = "Great", df = 56, t = 0.1326, p = 0.3467)
    (Name = "x3", Estimate = 1.56734e-13, Comments = "Bad", df = 2, t = 2.68e-16, p = 1.345e-16)
]


m = [0.11258244478647295  0.05664544616214151  0.38181274408522614  0.8197779704008801
     0.36831406658084287  0.12078054506961555  0.8151038332483567   0.6699313951612162
     0.3444540231363058   0.17957407667101322  0.2422083248151139   0.4530583319523316]
ct = CoefTable(m, ["Estimate", "Stderror", "df", "p"], [], 4)
@test sprint(show, ct) == """
──────────────────────────────────────────
     Estimate   Stderror        df       p
──────────────────────────────────────────
[1]  0.112582  0.0566454  0.381813  0.8198
[2]  0.368314  0.120781   0.815104  0.6699
[3]  0.344454  0.179574   0.242208  0.4531
──────────────────────────────────────────"""
@test length(ct) === 3
@test eltype(ct) ==
    NamedTuple{(:Estimate, :Stderror, :df, :p),
               Tuple{Float64,Float64,Float64,Float64}}
@test collect(ct) == [
    (Estimate = 0.11258244478647295, Stderror = 0.05664544616214151,
     df = 0.38181274408522614, p = 0.8197779704008801)
    (Estimate = 0.36831406658084287, Stderror = 0.12078054506961555,
     df = 0.8151038332483567, p = 0.6699313951612162)
    (Estimate = 0.3444540231363058, Stderror = 0.17957407667101322,
     df = 0.2422083248151139, p = 0.4530583319523316)
]

@test sprint(show, PValue(1.0)) == "1.0000"
@test sprint(show, PValue(1e-1)) == "0.1000"
if VERSION > v"1.6.0-DEV"
    @test sprint(show, PValue(1e-5)) == "<1e-04"
else
    @test sprint(show, PValue(1e-5)) == "<1e-4"
end
@test sprint(show, PValue(NaN)) == "NaN"
@test_throws ErrorException PValue(-0.1)
@test_throws ErrorException PValue(1.1)

@test sprint(show, TestStat(NaN)) == "NaN"
@test sprint(show, TestStat(1e-1)) == "0.10"
@test sprint(show, TestStat(1e-5)) == "0.00"
@test sprint(show, TestStat(π)) == "3.14"

@testset "Union{PValue, TestStat} is Real" begin
    vals = [0.0, Rational(1,3), NaN]
    for T in [PValue, TestStat],
        f in (==, <, ≤, >, ≥, isless, isequal),
        lhs in vals, rhs in vals
        # make sure that T behaves like a Real,
        # regardless of whether it's on the LHS, RHS or both
        @test f(T(lhs), T(rhs)) == f(lhs, rhs)
        @test f(lhs, T(rhs)) == f(lhs, rhs)
        @test f(T(lhs), rhs) == f(lhs, rhs)
    end

    # the (approximate) equality operators get a bit more attention
    for T in [PValue, TestStat]
        @test T(Rational(1,3)) ≈ T(1/3)
        @test Rational(1,3) ≈ T(1/3) atol=0.01
        @test T(Rational(1,3)) isa Real
        @test T(T(0.05)) === T(0.05)
        @test hash(T(0.05)) == hash(0.05)
        @test hash(T(0.05), UInt(42)) == hash(0.05, UInt(42))
    end
end

@test sprint(showerror, ConvergenceException(10)) == "failure to converge after 10 iterations."

@test sprint(showerror, ConvergenceException(10, 0.2, 0.1)) ==
    "failure to converge after 10 iterations. Last change (0.2) was greater than tolerance (0.1)."

@test sprint(showerror, ConvergenceException(10, 0.2, 0.1, "Try changing maxIter.")) ==
    "failure to converge after 10 iterations. Last change (0.2) was greater than tolerance (0.1). Try changing maxIter."

err = @test_throws ArgumentError ConvergenceException(10,.1,.2)
@test err.value.msg == "Change must be greater than tol."

struct MyStatisticalModel <: StatisticalModel
end

StatsAPI.vcov(::MyStatisticalModel) = [1 2; 3 4]
StatsAPI.loglikelihood(::MyStatisticalModel) = 3
StatsAPI.nullloglikelihood(::MyStatisticalModel) = 4
StatsAPI.deviance(::MyStatisticalModel) = 25
StatsAPI.nulldeviance(::MyStatisticalModel) = 40
StatsAPI.dof(::MyStatisticalModel) = 5
StatsAPI.nobs(::MyStatisticalModel) = 100

@testset "StatisticalModel" begin
    m = MyStatisticalModel()

    @test stderror(m) == [1, 2]
    @test aic(m) == 4
    @test aicc(m) ≈ 4.638297872340425
    @test bic(m) ≈ 17.02585092994046
    @test r2(m, :McFadden) ≈ 0.25
    @test r2(m, :CoxSnell) ≈ -0.020201340026755776
    @test r2(m, :Nagelkerke) ≈ 0.24255074155803877
    @test r2(m, :devianceratio) ≈ 0.375

    @test_throws Union{ErrorException, ArgumentError} r2(m, :err)
    @test_throws MethodError r2(m)
    @test adjr2(m, :McFadden) ≈ 1.5
    @test adjr2(m, :devianceratio) ≈ 0.3486842105263158
    @test_throws Union{ErrorException, ArgumentError} adjr2(m, :err)

    @test r2 === r²
    @test adjr2 === adjr²
end

struct MyRegressionModel <: RegressionModel
end

StatsAPI.modelmatrix(::MyRegressionModel) = [1 2; 3 4]

@testset "TestRegressionModel" begin
    m = MyRegressionModel()

    @test crossmodelmatrix(m) == [10 14; 14 20]
    @test crossmodelmatrix(m) isa Symmetric
end

@testset "StatsAPI model reexports" begin
    for f in (fitted, response, responsename, meanresponse,
              modelmatrix, crossmodelmatrix, leverage, cooksdistance, residuals,
              predict, predict!, dof_residual, coef, coefnames, coeftable, confint,
              deviance, islinear, nulldeviance, loglikelihood, nullloglikelihood,
              loglikelihood, loglikelihood, score, nobs, dof, mss, rss,
              informationmatrix, stderror, vcov, weights, isfitted, fit, fit!,
              aic, aicc, bic, r2, r², adjr2, adjr²)
        @test f isa Function
    end
    # Defined but not reexported
    @test StatsBase.params isa Function
    @test StatsBase.params! isa Function
end