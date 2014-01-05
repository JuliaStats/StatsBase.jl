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

c = counts(x, 1:5, w)
@test size(c) == (5,)
c0 = Float64[sum(w.values[x .== i]) for i in 1 : 5]
@test_approx_eq c c0
@test_approx_eq counts(x+1, 2:6, w) c0
@test_approx_eq proportions(x, 1:5, w) (c0 ./ sum(w))

# 2D integer counts

x = rand(1:4, n)
y = rand(1:5, n)
w = weights(rand(n))

c = counts(x, y, (1:4, 1:5))
@test size(c) == (4, 5)
c0 = Int[nnz((x .== i) & (y .== j)) for i in 1 : 4, j in 1 : 5]
@test c == c0
@test counts(x+2, y+3, (3:6, 4:8)) == c0
@test_approx_eq proportions(x, y, (1:4, 1:5)) (c0 ./ n)

c = counts(x, y, (1:4, 1:5), w)
@test size(c) == (4, 5)
c0 = Float64[sum(w.values[(x .== i) & (y .== j)]) for i in 1 : 4, j in 1 : 5]
@test_approx_eq c c0
@test_approx_eq counts(x+2, y+3, (3:6, 4:8), w) c0
@test_approx_eq proportions(x, y, (1:4, 1:5), w) (c0 ./ sum(w))

# count map

x = ["a", "b", "a", "a", "b", "c"]
w = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

cm = counts(x)
@test cm["a"] == 3
@test cm["b"] == 2
@test cm["c"] == 1
pm = proportions(x)
@test_approx_eq pm["a"] (1/2)
@test_approx_eq pm["b"] (1/3)
@test_approx_eq pm["c"] (1/6)

cm = counts(x, weights(w))
@test cm["a"] == 5.5
@test cm["b"] == 4.5
@test cm["c"] == 3.5
pm = proportions(x, weights(w))
@test_approx_eq pm["a"] (5.5 / 13.5)
@test_approx_eq pm["b"] (4.5 / 13.5)
@test_approx_eq pm["c"] (3.5 / 13.5)

x = ["a", "b", "a", "a", "b", "c", "b"]
y = [ 1,   2,   2,   2,   1,   2,   1]
w = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

cm = counts(x, y)
@test cm[("a", 1)] == 1
@test cm[("a", 2)] == 2
@test cm[("b", 1)] == 2
@test cm[("b", 2)] == 1
@test cm[("c", 2)] == 1

pm = proportions(x, y)
@test_approx_eq pm[("a", 1)] 1/7
@test_approx_eq pm[("a", 2)] 2/7
@test_approx_eq pm[("b", 1)] 2/7
@test_approx_eq pm[("b", 2)] 1/7
@test_approx_eq pm[("c", 2)] 1/7







