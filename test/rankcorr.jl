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

# corkendall

@test corkendall(x1, y) == -1 / sqrt(90)
@test corkendall(x2, y) == -1 / sqrt(72)
@test corkendall(X, y) == [-1 / sqrt(90), -1 / sqrt(72)]
@test corkendall(y, X) == [-1 / sqrt(90) -1 / sqrt(72)]

n = 100_000 # so vectors tested are 5 times this length
@test corkendall(repeat(x1, n), repeat(y, n)) ≈ -1 / sqrt(90)
@test corkendall(repeat(x2, n), repeat(y, n)) ≈ -1 / sqrt(72)
@test corkendall(repeat(X, n), repeat(y, n)) ≈ [-1 / sqrt(90), -1 / sqrt(72)]
@test corkendall(repeat(y, n), repeat(X, n)) ≈ [-1 / sqrt(90) -1 / sqrt(72)]

c11 = corkendall(x1, x1)
c12 = corkendall(x1, x2)
c22 = corkendall(x2, x2)

@test c11 == 1.0
@test c22 == 1.0
@test c12 == 6 / sqrt(80)

@test corkendall(X, X) ≈ [c11 c12; c12 c22]
@test corkendall(X)    ≈ [c11 c12; c12 c22]

@test corkendall(repeat(X, n), repeat(X, n)) ≈ [c11 c12; c12 c22]
@test corkendall(repeat(X, n))                  ≈ [c11 c12; c12 c22]


@test corkendall(collect(1:100_000), collect(1:100_000)) == 1.0
@test corkendall(collect(1:100_000), reverse(collect(1:100_000))) == -1.0

@test corkendall(repeat([0,1,1,0], n), repeat([1,0,1,0], n)) == 0.0

z = [1  1  1;
    1  1  2;
    1  2  2;
    1  2  2;
    1  2  1;
    2  1  2;
    1  1  2;
    2  2  2]

@test corkendall(z)   == [1 0 1 / 3; 0 1 0;1 / 3 0 1]
@test corkendall(z, z)   == [1 0 1 / 3; 0 1 0;1 / 3 0 1]
@test corkendall(z[:,1], z) == [1 0 1 / 3]
@test corkendall(z, z[:,1]) == [1;0;1 / 3]

z = float(z)
@test corkendall(z)   == [1 0 1 / 3; 0 1 0;1 / 3 0 1]
@test corkendall(z, z)   == [1 0 1 / 3; 0 1 0;1 / 3 0 1]
@test corkendall(z[:,1], z) == [1 0 1 / 3]
@test corkendall(z, z[:,1]) == [1;0;1 / 3]

w = repeat(z, n)
@test corkendall(w)   == [1 0 1 / 3; 0 1 0;1 / 3 0 1]
@test corkendall(w, w)   == [1 0 1 / 3; 0 1 0;1 / 3 0 1]
@test corkendall(w[:,1], w) == [1 0 1 / 3]
@test corkendall(w, w[:,1]) == [1;0;1 / 3]

