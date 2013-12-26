# Unit testing of intstats

using Stats
using Base.Test

n = 5000

# 1D integer counts

x = rand(1:5, n)
c = counts(x, 1:5)
@test size(c) == (5,)
c0 = Int[nnz(x .== i) for i in 1 : 5]
@test c == c0
@test counts(x+1, 2:6) == c0


x = rand(1:5, n)
w = rand(n)
c = wcounts(x, w, 1:5)
@test size(c) == (5,)
c0 = Float64[sum(w[x .== i]) for i in 1 : 5]
@test_approx_eq c c0
@test_approx_eq wcounts(x+1, w, 2:6) c0

# x = rand(1:4, n)
# y = rand(1:5, n)
# c = counts(x, y)
# @test size(c) == (4, 5)
# c0 = Int[nnz((x .== i) & (y .== j)) for i in 1 : 4, j in 1 : 5]
# @test c == c0

# x = rand(1:4, n)
# y = rand(1:5, n)
# w = rand(n)
# c = wcounts(x, y, w)
# @test size(c) == (4, 5)
# c0 = Float64[sum(w[(x .== i) & (y .== j)]) for i in 1 : 4, j in 1 : 5]
# @test_approx_eq c c0
