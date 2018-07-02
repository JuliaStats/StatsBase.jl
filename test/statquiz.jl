# Test taken from http://www.stanford.edu/~clint/bench/wilk.txt

using Test
using DataFrames
using StatsBase
using GLM

using Printf

testeps = sqrt(eps())

nasty = DataFrame(    label = ["One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"],
                    x = collect(1.:9),
                    zero = fill(0.0,9),
                    miss = fill(NA, 9),
                    big = 99999990.0 + collect(1:9),
                    little = (99999990.0 + collect(1:9))/10^8,
                    huge = collect(1.:9)*1e12,
                    tiny = collect(1.:9)*1e-12,
                    round = collect(0.5:8.5))

println(nasty)
println("\nII Real Numbers:\nII A")
print("Test rounding: ")
@test [@sprintf("%1.0f", x) for x in nasty[:round]] == ["1","2","3","4","5","6","7","8","9"]
println("OK")
print("Test math: ")
@test round(Int, 2.6*7 - 0.2) == 18
@test 2 - round(Int, exp(log(sqrt(2)*sqrt(2)))) == 0
@test round(Int, 3 - exp(log(sqrt(2)*sqrt(2)))) == 1
println("OK")

print("Test means: ")
for vars in names(nasty)[2:end]
    if vars == :miss
    @test isna(mean(nasty[vars]))
    else
    @test mean(nasty[vars]) ≈ nasty[vars][5]
    end
end
println("OK")

print("Test standard deviation: ")
for vars in names(nasty)[[2;5:9]]
#    @test (@sprintf("%.9e", std(nasty[vars])))[1:10] == "2.73861278"
    @test repr(std(nasty[vars]))[1:10] == "2.73861278"
end
println("OK")
# Note: Original had one more digit but little fails then. R gives the same as Julia. Stata and SAS less precise returning 2.738612714381 and 2.738612780003404")

println("\nII D")
print("Test correlation: ")
cn = names(nasty)[[2;5:9]]
for i in 1:5
    for j = i+1:6
    @test cor(nasty[cn[i]], nasty[cn[j]]) ≈ 1
    end
end
println("OK")

print("Test spearman correlation: ")
cn = names(nasty)[[2;5:9]]
for i in 1:5
    for j = i+1:6
    @test corspearman(nasty[cn[i]], nasty[cn[j]]) ≈ 1
    end
end
println("OK")

println("\nII F")
print("Testing regression: ")
ctable = coeftable(lm(big ~ x, nasty))
@test typeof(ctable) == CoefTable
# conversion to Vector{Float64} can be removed when 0.4 is no longer supported
@test Vector{Float64}(ctable.cols[1]) ≈ [99999990, 1]

@test sprint(show, ctable) == """\
             Estimate Std.Error t value Pr(>|t|)
(Intercept)     1.0e8       0.0     Inf   <1e-99
x                 1.0       0.0     Inf   <1e-99
"""

println("OK")

println("\nIV Regression:\nIV A")
nasty[:x1] = nasty[:x]
nasty[:x2] = nasty[:x].^2
nasty[:x3] = nasty[:x].^3
nasty[:x4] = nasty[:x].^4
nasty[:x5] = nasty[:x].^5
nasty[:x6] = nasty[:x].^6
nasty[:x7] = nasty[:x].^7
nasty[:x8] = nasty[:x].^8
nasty[:x9] = nasty[:x].^9

## Is it intended that the least squares problem be overdetermined in the lm fit?
## n = 9 and p = 10 because of the implicit intercept.
lm(x1~x2+x3+x4+x5+x6+x7+x8+x9, nasty)
@test coef(lm(x~x, nasty)) ≈ [0,1]
println("OK")
