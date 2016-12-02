using StatsBase
using Base.Test

srand(10)
v1 = rand(3)
v2 = ["Good", "Great", "Bad"]
v3 = rand(Int8, 3)
v4 = [StatsBase.PValue(rand()./10000) for i in 1:3]
m = rand(3,4)
@test sprint(show, CoefTable(Any[v1, v2, v3, v4],
    ["Estimate", "Comments", "df", "p"],
    ["x1", "x2", "x3"])) == """\
     Estimate Comments  df     p
x1   0.112582     Good  88 <1e-4
x2   0.368314    Great -90 <1e-4
x3   0.344454      Bad -80 <1e-4
"""

@test sprint(show, CoefTable(m, ["Estimate", "Stderr", "df", "p"],
    ["x1", "x2", "x3"], 4)) == """\
     Estimate   Stderr       df      p
x1   0.819778 0.844007 0.923676 0.1717
x2   0.669931  0.67919 0.066098 0.4204
x3   0.453058  0.72525 0.999172 0.5567
"""

@test sprint(show, StatsBase.PValue(1.0)) == "1.0000"
@test sprint(show, StatsBase.PValue(1e-1)) == "0.1000"
@test sprint(show, StatsBase.PValue(1e-5)) == "<1e-4"
@test sprint(show, StatsBase.PValue(NaN)) == "NaN"
@test_throws ErrorException StatsBase.PValue(-0.1)
@test_throws ErrorException StatsBase.PValue(1.1)

try
    throw(ConvergenceException(10))
catch ex
    @test sprint(showerror, ex) == "failure to converge after 10 iterations"
end
