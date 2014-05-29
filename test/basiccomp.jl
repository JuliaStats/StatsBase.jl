# test basic computation

using StatsBase
using Base.Test

x0 = [1.0, 3.0, 4.0, 6.0, 9.0]
y0 = [5.0, 7.0, 2.5, 3.5, 4.0]

x = copy(x0)
@test is(negate!(x), x)
@test x == -x0

x = copy(x0)
y = copy(y0)
@test is(add!(y, x), y)
@test y == y0 + x0

y = copy(y0)
@test is(add!(y, 1), y)
@test y == y0 .+ 1

x = copy(x0)
y = copy(y0)
@test is(subtract!(y, x), y)
@test y == y0 - x0

y = copy(y0)
@test is(subtract!(y, 1), y)
@test y == y0 .- 1

y = copy(y0)
x = copy(x0)
@test is(addscale!(y, x, 2.0), y)
@test y == y0 + 2.0 * x

y = y0
x = x0
@test addscale(y, x, 2.0) == y0 + x0 * 2.0
@test addscale(y, 1:5, 2.0) == y0 + [2, 4, 6, 8, 10]

ly = rand(10^4)
lx = rand(10^4)
@test addscale(ly, lx, 2.0) == ly + lx * 2.0


z = Float64[]
x = [1.0, -3.0, 4.0, -6.0, 9.0]
y = [-5.0, 7.0, 2.5, -3.5, 4.0]

@test sumabs(z) == 0
@test sumabs(x) == sum(abs(x))
@test maxabs(z) == 0
@test maxabs(x) == maximum(abs(x))
@test sumabs2(z) == 0
@test sumabs2(x) == sum(abs2(x))

@test sumabsdiff(x, y) == sum(abs(x - y))
@test sumabsdiff(x, 2.0) == sum(abs(x .- 2.0))
@test sumabsdiff(2.0, x) == sum(abs(x .- 2.0))
@test maxabsdiff(x, y) == maximum(abs(x - y))
@test maxabsdiff(x, 2.0) == maximum(abs(x .- 2.0))
@test maxabsdiff(2.0, x) == maximum(abs(x .- 2.0))
@test sumabs2diff(x, y) == sum(abs2(x - y))
@test sumabs2diff(x, 2.0) == sum(abs2(x .- 2.0))
@test sumabs2diff(2.0, x) == sum(abs2(x .- 2.0))


@test xlogx(0) === 0.0
@test_approx_eq xlogx(2) 2.0 * log(2.0)

@test xlogy(0, 1) === 0.0
@test_approx_eq xlogy(2, 3) 2.0 * log(3.0)

@test_approx_eq logistic(2) 1.0 / (1.0 + exp(-2.0))
@test_approx_eq logit(0.5) 0.0
@test_approx_eq logit(logistic(2)) 2.0

@test_approx_eq softplus(2.0) log(1.0 + exp(2.0))
@test_approx_eq softplus(-2.0) log(1.0 + exp(-2.0))
@test_approx_eq softplus(10000) 10000.0
@test_approx_eq softplus(-10000) 0.0

@test_approx_eq invsoftplus(softplus(2.0)) 2.0
@test_approx_eq invsoftplus(softplus(-2.0)) -2.0

@test_approx_eq logsumexp(2.0, 3.0) log(exp(2.0) + exp(3.0))
@test_approx_eq logsumexp(10002, 10003) 10000 + logsumexp(2.0, 3.0)

@test_approx_eq logsumexp([1.0, 2.0, 3.0]) 3.40760596444438
@test_approx_eq logsumexp([1.0, 2.0, 3.0] .+ 1000.) 1003.40760596444438

x = [1.0, 2.0, 3.0]
r = exp(x) ./ sum(exp(x))
@test_approx_eq softmax([1.0, 2.0, 3.0]) r
softmax!(x)
@test_approx_eq x r

