using StatsBase
using Base.Test

a = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0]
x = [3.0, 1.0, 2.0, 4.0, 4.0, 2.0, 5.0, 4.0]  # x is a permutated version of a

@test ordinalrank(a) == [1, 2, 3, 4, 5, 6, 7, 8]
srand(123);
@test ordinalrank(a, rand=true) == [1, 3, 2, 4, 5, 6, 7, 8]
@test ordinalrank(a, rand=true) == [1, 2, 3, 4, 6, 7, 5, 8]
@test ordinalrank(x) == [4, 1, 2, 5, 6, 3, 8, 7]
@test ordinalrank(x, rand=true) == [4, 1, 2, 7, 6, 3, 8, 5]
@test ordinalrank(x, rand=true) == [4, 1, 2, 5, 7, 3, 8, 6]

@test competerank(a) == [1, 2, 2, 4, 5, 5, 5, 8]
@test competerank(x) == [4, 1, 2, 5, 5, 2, 8, 5]

@test denserank(a) == [1, 2, 2, 3, 4, 4, 4, 5]
@test denserank(x) == [3, 1, 2, 4, 4, 2, 5, 4]

@test tiedrank(a) == [1.0, 2.5, 2.5, 4.0, 6.0, 6.0, 6.0, 8.0]
@test tiedrank(x) == [4.0, 1.0, 2.5, 6.0, 6.0, 2.5, 8.0, 6.0]
