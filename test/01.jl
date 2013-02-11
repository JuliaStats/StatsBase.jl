require("test")
require("test/runtests.jl")

using Stats

@test abs(autocor([1, 2, 3, 4, 5]) - 1.0) < 10e-8

@test iqr([1, 2, 3, 4, 5]) == [2.0, 4.0]

z = [true, true, false, false, true, false, true, true, true]
values, lengths = rle(z)
@test values == [true, false, true, false, true]
@test lengths == [2, 2, 1, 1, 3]
@test inverse_rle(values, lengths) == z

z = [true, true, false, false, true, false, true, true, true, false]
values, lengths = rle(z)
@test values == [true, false, true, false, true, false]
@test lengths == [2, 2, 1, 1, 3, 1]
@test inverse_rle(values, lengths) == z

m = [1 0; 0 1]
d = [0.0 sqrt(2); sqrt(2) 0.0]
@test norm(distances(m) - d) < 10e-8

m = [3.0 1.0; 5.0 1.0]
d = [0.0 2.0; 2.0 0.0]
@test norm(distances(m) - d) < 10e-8

m = [1 0 0; 0 1 0 ; 1 0 1]
d = [0.0 sqrt(2) 1.0; sqrt(2) 0.0 sqrt(3); 1.0 sqrt(3) 0.0]
@test norm(distances(m) - d) < 10e-8

# X = [1 0; 2 1; 3 0; 4 1; 5 10]
# y = [5, 3, 4, 2, 5]
# @assert_approx_eq cov_spearman(X, y)[1] cov_spearman(X[:,1],y)
# @assert_approx_eq cov_spearman(X) cov_spearman(X, X)
# @assert_approx_eq cov_spearman(X, y) [-0.25, -0.1875]
