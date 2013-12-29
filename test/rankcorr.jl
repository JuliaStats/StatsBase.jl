# Test rank correlation

using Stats
using Base.Test

X = [1 0; 2 1; 3 0; 4 1; 5 10]

x1 = X[:,1]
x2 = X[:,2]
y = [5, 3, 4, 2, 5]

# cor_spearman

@test_approx_eq cor_spearman(x1, y) -0.102597835208515
@test_approx_eq cor_spearman(x2, y) -0.081110710565381

@test_approx_eq cor_spearman(X, y) [-0.102597835208515, -0.081110710565381]
@test_approx_eq cor_spearman(y, X) [-0.102597835208515 -0.081110710565381]

c11 = cor_spearman(x1, x1)
c12 = cor_spearman(x1, x2)
c22 = cor_spearman(x2, x2)
@test_approx_eq cor_spearman(X, X) [c11 c12; c12 c22]
@test_approx_eq cor_spearman(X)    [c11 c12; c12 c22]

# @test_approx_eq cor_kendall(X, y) [-0.105409255338946, -0.117851130197758]

