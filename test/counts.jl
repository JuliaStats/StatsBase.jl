# Unit testing of intstats

using Stats
using Base.Test

n = 5000

# 1D integer counts

x = rand(1:5, n)
w = weights(rand(n))

c = counts(x, 1:5)
@test size(c) == (5,)
c0 = Int[nnz(x .== i) for i in 1 : 5]
@test c == c0
@test counts(x+1, 2:6) == c0
@test_approx_eq proportions(x, 1:5) (c0 ./ n)

c = counts(x, w, 1:5)
@test size(c) == (5,)
c0 = Float64[sum(w.values[x .== i]) for i in 1 : 5]
@test_approx_eq c c0
@test_approx_eq counts(x+1, w, 2:6) c0
@test_approx_eq proportions(x, w, 1:5) (c0 ./ sum(w))

# 2D integer counts

x = rand(1:4, n)
y = rand(1:5, n)
w = weights(rand(n))

c = counts(x, y, 1:4, 1:5)
@test size(c) == (4, 5)
c0 = Int[nnz((x .== i) & (y .== j)) for i in 1 : 4, j in 1 : 5]
@test c == c0
@test counts(x+2, y+3, 3:6, 4:8) == c0
@test_approx_eq proportions(x, y, 1:4, 1:5) (c0 ./ n)

c = counts(x, y, w, 1:4, 1:5)
@test size(c) == (4, 5)
c0 = Float64[sum(w.values[(x .== i) & (y .== j)]) for i in 1 : 4, j in 1 : 5]
@test_approx_eq c c0
@test_approx_eq counts(x+2, y+3, w, 3:6, 4:8) c0
@test_approx_eq proportions(x, y, w, 1:4, 1:5) (c0 ./ sum(w))
