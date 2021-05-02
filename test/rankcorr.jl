using StatsBase
using Test

X = Float64[1 0; 2 1; 3 0; 4 1; 5 10]
Y = Float64[5 5 6; 3 4 1; 4 0 4; 2 6 1; 5 7 10]

x1 = X[:,1]
x2 = X[:,2]
y = Y[:,1]

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

@test corspearman(X, Y) ==
     [corspearman(X[:,i], Y[:,j]) for i in axes(X, 2), j in axes(Y, 2)]

# corkendall

# Check error, handling of NaN, Inf etc
@test_throws ErrorException("Vectors must have same length") corkendall([1,2,3,4], [1,2,3])
@test isnan(corkendall([1,2], [3,NaN]))
@test isnan(corkendall([1,1,1], [1,2,3]))
@test corkendall([-Inf,-0.0,Inf],[1,2,3]) == 1.0

# Test, with exact equality, some known results. 
# RealVector, RealVector
@test corkendall(x1, y) == -1/sqrt(90)
@test corkendall(x2, y) == -1/sqrt(72)
# RealMatrix, RealVector
@test corkendall(X, y)  == [-1/sqrt(90), -1/sqrt(72)]
# RealVector, RealMatrix
@test corkendall(y, X)  == [-1/sqrt(90) -1/sqrt(72)]

# n = 78_000 tests for overflow errors on 32 bit
# Testing for overflow errors on 64bit would require n be too large for practicality
# This also tests merge_sort! since n is (much) greater than SMALL_THRESHOLD.
n = 78_000
# Test with many repeats
@test corkendall(repeat(x1, n), repeat(y, n)) ≈ -1/sqrt(90)
@test corkendall(repeat(x2, n), repeat(y, n)) ≈ -1/sqrt(72)
@test corkendall(repeat(X, n), repeat(y, n))  ≈ [-1/sqrt(90), -1/sqrt(72)]
@test corkendall(repeat(y, n), repeat(X, n))  ≈ [-1/sqrt(90) -1/sqrt(72)]
@test corkendall(repeat([0,1,1,0], n), repeat([1,0,1,0], n)) == 0.0

# Test with no repeats, note testing for exact equality
@test corkendall(collect(1:n), collect(1:n))          == 1.0
@test corkendall(collect(1:n), reverse(collect(1:n))) == -1.0

# All elements identical should yield NaN
@test isnan(corkendall(repeat([1], n), collect(1:n)))

c11 = corkendall(x1, x1)
c12 = corkendall(x1, x2)
c22 = corkendall(x2, x2)

# RealMatrix, RealMatrix
@test corkendall(X, X) ≈ [c11 c12; c12 c22]
# RealMatrix
@test corkendall(X)    ≈ [c11 c12; c12 c22]

@test c11 == 1.0
@test c22 == 1.0
@test c12 == 3/sqrt(20)

# Finished testing for overflow, so redefine n for speedier tests
n = 100

@test corkendall(repeat(X, n), repeat(X, n)) ≈ [c11 c12; c12 c22]
@test corkendall(repeat(X, n))               ≈ [c11 c12; c12 c22]

# All eight three-element permutations
z = [1  1  1;
     1  1  2;
     1  2  2;
     1  2  2;
     1  2  1;
     2  1  2;
     1  1  2;
     2  2  2]

@test corkendall(z)         == [1 0 1/3; 0 1 0; 1/3 0 1]
@test corkendall(z, z)      == [1 0 1/3; 0 1 0; 1/3 0 1]
@test corkendall(z[:,1], z) == [1 0 1/3]
@test corkendall(z, z[:,1]) == [1; 0; 1/3]

z = float(z)
@test corkendall(z)         == [1 0 1/3; 0 1 0; 1/3 0 1]
@test corkendall(z, z)      == [1 0 1/3; 0 1 0; 1/3 0 1]
@test corkendall(z[:,1], z) == [1 0 1/3]
@test corkendall(z, z[:,1]) == [1; 0; 1/3]

w = repeat(z, n)
@test corkendall(w)         == [1 0 1/3; 0 1 0; 1/3 0 1]
@test corkendall(w, w)      == [1 0 1/3; 0 1 0; 1/3 0 1]
@test corkendall(w[:,1], w) == [1 0 1/3]
@test corkendall(w, w[:,1]) == [1; 0; 1/3]

StatsBase.midpoint(1,10)        == 5
StatsBase.midpoint(1,widen(10)) == 5


# NaN handling

Xnan = copy(X)
Xnan[1,1] = NaN
Ynan = copy(Y)
Ynan[2,1] = NaN

for f in (corspearman, corkendall)
     @test isnan(f([1.0, NaN, 2.0], [2.0, 1.0, 3.4]))
     @test all(isnan, f([1.0, NaN], [1 2; 3 4]))
     @test all(isnan, f([1 2; 3 4], [1.0, NaN]))
     @test isequal(f([1 NaN; NaN 4]), [1 NaN; NaN 1])
     @test all(isnan, f([1 NaN; NaN 4], [1 NaN; NaN 4]))
     @test all(isnan, f([1 NaN; NaN 4], [NaN 1; NaN 4]))

     @test isequal(f(Xnan, Ynan),
                   [f(Xnan[:,i], Ynan[:,j]) for i in axes(Xnan, 2), j in axes(Ynan, 2)])
     @test isequal(f(Xnan),
                   [i == j ? 1.0 : f(Xnan[:,i], Xnan[:,j])
                    for i in axes(Xnan, 2), j in axes(Xnan, 2)])
     for k in 1:2
          @test isequal(f(Xnan[:,k], Ynan),
                        [f(Xnan[:,k], Ynan[:,j]) for i in 1:1, j in axes(Ynan, 2)])
          # TODO: fix corkendall (PR#659)
          if f === corspearman
               @test isequal(f(Xnan, Ynan[:,k]),
                             [f(Xnan[:,i], Ynan[:,k]) for i in axes(Xnan, 2), j in 1:1])
          else
               @test isequal(f(Xnan, Ynan[:,k]),
                             [f(Xnan[:,i], Ynan[:,k]) for i in axes(Xnan, 2)])
          end
     end
end


# Wrong dimensions

@test_throws DimensionMismatch corspearman([1], [1, 2])
@test_throws DimensionMismatch corspearman([1], [1 2; 3 4])
@test_throws DimensionMismatch corspearman([1 2; 3 4], [1])
@test_throws ArgumentError corspearman([1 2; 3 4: 4 6], [1 2; 3 4])

# TODO: fix corkendall to match corspearman (PR#659)
@test_throws ErrorException corkendall([1], [1, 2])
@test_throws ErrorException corkendall([1], [1 2; 3 4])
@test_throws ErrorException corkendall([1 2; 3 4], [1])
@test_throws ArgumentError corkendall([1 2; 3 4: 4 6], [1 2; 3 4])