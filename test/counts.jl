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

x = ["a", "b", "c", "a", "c", "c"]
w = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
levels = ["a", "b", "c"]

@test counts(x, levels) == [2, 1, 3] 
@test proportions(x, levels) == [2, 1, 3] / 6
@test counts(x, levels, weights(w)) == [5.0, 2.0, 14.0]
@test proportions(x, levels, weights(w)) == [5.0, 2.0, 14.0] / 21.0

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

x = ["a", "b", "a", "a", "b", "c", "b"]
y = [ 1,   2,   2,   2,   1,   2,   1]
w = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

levels = (["a", "b", "c"], 1:2)
@test counts(x, y, levels) == [1 2; 2 1; 0 1]
@test proportions(x, y, levels) == [1 2; 2 1; 0 1] / 7
@test counts(x, y, levels, weights(w)) == [1.0 7.0; 12.0 2.0; 0.0 6.0]
@test proportions(x, y, levels, weights(w)) == [1.0 7.0; 12.0 2.0; 0.0 6.0] / 28.0

# count map

x = ["a", "b", "a", "a", "b", "c"]
w = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

cm = countmap(x)
@test cm["a"] == 3
@test cm["b"] == 2
@test cm["c"] == 1
pm = proportionmap(x)
@test_approx_eq pm["a"] (1/2)
@test_approx_eq pm["b"] (1/3)
@test_approx_eq pm["c"] (1/6)

cm = countmap(x, weights(w))
@test cm["a"] == 5.5
@test cm["b"] == 4.5
@test cm["c"] == 3.5
pm = proportionmap(x, weights(w))
@test_approx_eq pm["a"] (5.5 / 13.5)
@test_approx_eq pm["b"] (4.5 / 13.5)
@test_approx_eq pm["c"] (3.5 / 13.5)

x = ["a", "b", "a", "a", "b", "c", "b"]
y = [ 1,   2,   2,   2,   1,   2,   1]
w = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

cm = countmap(x, y)
@test cm[("a", 1)] == 1
@test cm[("a", 2)] == 2
@test cm[("b", 1)] == 2
@test cm[("b", 2)] == 1
@test cm[("c", 2)] == 1

pm = proportionmap(x, y)
@test_approx_eq pm[("a", 1)] 1/7
@test_approx_eq pm[("a", 2)] 2/7
@test_approx_eq pm[("b", 1)] 2/7
@test_approx_eq pm[("b", 2)] 1/7
@test_approx_eq pm[("c", 2)] 1/7

cm = countmap(x, y, weights(w))
@test cm[("a", 1)] == 1.0
@test cm[("a", 2)] == 7.0
@test cm[("b", 1)] == 12.0
@test cm[("b", 2)] == 2.0
@test cm[("c", 2)] == 6.0

pm = proportionmap(x, y, weights(w))
@test pm[("a", 1)] == 1.0 / 28.0
@test pm[("a", 2)] == 7.0 / 28.0
@test pm[("b", 1)] == 12.0 / 28.0
@test pm[("b", 2)] == 2.0 / 28.0
@test pm[("c", 2)] == 6.0 / 28.0

