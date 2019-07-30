using StatsBase
using Test, Random

v1 = [1.45666, -23.14, 1.56734e-13]
v2 = ["Good", "Great", "Bad"]
v3 = [1, 56, 2]
v4 = [-12.56, 0.1326, 2.68e-16]
v5 = [0.12, 0.3467, 1.345e-16]
@test sprint(show, CoefTable(Any[v1, v2, v3, v4, v5],
    ["Estimate", "Comments", "df", "t", "p"],
    ["x1", "x2", "x3"], 5, 4)) == """
───────────────────────────────────────────────
         Estimate  Comments  df       t       p
───────────────────────────────────────────────
x1    1.45666         Good    1  -12.56  0.1200
x2  -23.14            Great  56    0.13  0.3467
x3    1.56734e-13     Bad     2    0.00  <1e-15
───────────────────────────────────────────────"""

Random.seed!(10)
m = rand(3,4)
@test sprint(show, CoefTable(m, ["Estimate", "Stderror", "df", "p"], [], 4)) == """
──────────────────────────────────────────
     Estimate   Stderror        df       p
──────────────────────────────────────────
[1]  0.112582  0.0566454  0.381813  0.8198
[2]  0.368314  0.120781   0.815104  0.6699
[3]  0.344454  0.179574   0.242208  0.4531
──────────────────────────────────────────"""

@test sprint(show, StatsBase.PValue(1.0)) == "1.0000"
@test sprint(show, StatsBase.PValue(1e-1)) == "0.1000"
@test sprint(show, StatsBase.PValue(1e-5)) == "<1e-4"
@test sprint(show, StatsBase.PValue(NaN)) == "NaN"
@test_throws ErrorException StatsBase.PValue(-0.1)
@test_throws ErrorException StatsBase.PValue(1.1)

@test sprint(showerror, ConvergenceException(10)) == "failure to converge after 10 iterations."

@test sprint(showerror, ConvergenceException(10, 0.2, 0.1)) ==
    "failure to converge after 10 iterations. Last change (0.2) was greater than tolerance (0.1)."

@test sprint(showerror, ConvergenceException(10, 0.2, 0.1, "Try changing maxIter.")) ==
    "failure to converge after 10 iterations. Last change (0.2) was greater than tolerance (0.1). Try changing maxIter."

err = @test_throws ArgumentError ConvergenceException(10,.1,.2)
@test err.value.msg == "Change must be greater than tol."
