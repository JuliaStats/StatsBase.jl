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

