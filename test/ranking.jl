using StatsBase
using Compat
using Compat.Test

a = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0]
x = [3.0, 1.0, 2.0, 4.0, 4.0, 2.0, 5.0, 4.0]  # x is a permutated version of a
s = ["c", "a", "b", "d", "d", "b", "e", "d"] # s is a vector of strings ordered like x

@test ordinalrank(a) == [1, 2, 3, 4, 5, 6, 7, 8]
@test ordinalrank(x) == [4, 1, 2, 5, 6, 3, 8, 7]
@test ordinalrank(s) == ordinalrank(x)
@test ordinalrank(x, rev = true) == ordinalrank(-x)
@test ordinalrank(x, lt = (x, y) -> isless(y, x)) == ordinalrank(-x)

@test competerank(a) == [1, 2, 2, 4, 5, 5, 5, 8]
@test competerank(x) == [4, 1, 2, 5, 5, 2, 8, 5]
@test competerank(s) == competerank(x)
@test competerank(x, rev = true) == competerank(-x)
@test competerank(x, lt = (x, y) -> isless(y, x)) == competerank(-x)

@test denserank(a) == [1, 2, 2, 3, 4, 4, 4, 5]
@test denserank(x) == [3, 1, 2, 4, 4, 2, 5, 4]
@test denserank(s) == denserank(x)
@test denserank(x, rev = true) == denserank(-x)
@test denserank(x, lt = (x, y) -> isless(y, x)) == denserank(-x)

@test tiedrank(a) == [1.0, 2.5, 2.5, 4.0, 6.0, 6.0, 6.0, 8.0]
@test tiedrank(x) == [4.0, 1.0, 2.5, 6.0, 6.0, 2.5, 8.0, 6.0]
@test tiedrank(s) == tiedrank(x)
@test tiedrank(x, rev = true) == tiedrank(-x)
@test tiedrank(x, lt = (x, y) -> isless(y, x)) == tiedrank(-x)
