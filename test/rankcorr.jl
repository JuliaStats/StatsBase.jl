using StatsBase
using Test

X = Float64[1 0; 2 1; 3 0; 4 1; 5 10]

x1 = X[:,1]
x2 = X[:,2]
y = [5, 3, 4, 2, 5]

# corspearman

@test corspearman(x1, y) ≈ -0.102597835208515
@test corspearman(x2, y) ≈ -0.081110710565381

@test corspearman(X, y) ≈ [-0.102597835208515, -0.081110710565381]
@test corspearman(y, X) ≈ [-0.102597835208515 -0.081110710565381]

c11 = corspearman(x1, x1)
c12 = corspearman(x1, x2)
c22 = corspearman(x2, x2)
@test c11 ≈ 1.0
@test c22 ≈ 1.0
@test corspearman(X, X) ≈ [c11 c12; c12 c22]
@test corspearman(X)    ≈ [c11 c12; c12 c22]


# corkendall

@test corkendall(x1, y) ≈ -0.105409255338946
@test corkendall(x2, y) ≈ -0.117851130197758

@test corkendall(X, y) ≈ [-0.105409255338946, -0.117851130197758]
@test corkendall(y, X) ≈ [-0.105409255338946 -0.117851130197758]

c11 = corkendall(x1, x1)
c12 = corkendall(x1, x2)
c22 = corkendall(x2, x2)

@test c11 ≈ 1.0
@test c22 ≈ 1.0
@test corkendall(X, X) ≈ [c11 c12; c12 c22]
@test corkendall(X)    ≈ [c11 c12; c12 c22]
