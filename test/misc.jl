using StatsBase
using SparseArrays, Test

# rle & inverse_rle

z = [1, 1, 2, 2, 2, 3, 1, 2, 2, 3, 3, 3, 3]
(vals, lens) = rle(z)
@test vals == [1, 2, 3, 1, 2, 3]
@test lens == [2, 3, 1, 1, 2, 4]
@test inverse_rle(vals, lens) == z

z = [true, true, false, false, true, false, true, true, true]
vals, lens = rle(z)
@test vals == [true, false, true, false, true]
@test lens == [2, 2, 1, 1, 3]
@test inverse_rle(vals, lens) == z

# levelsmap
a = [1, 1, 2, 2, 2, 3, 1, 2, 2, 3, 3, 3, 3, 2]
b = [true, false, false, true, false, true, true, false]

@test levelsmap(a) == Dict(2=>2, 3=>3, 1=>1)
@test levelsmap(b) == Dict(false=>2, true=>1)

# indicatormat

II = [false true  false false false;
      true  false false false true;
      false false true  true  false]

x = [2, 1, 3, 3, 2]
@test indicatormat(x, 3) == II
@test Matrix(indicatormat(x, 3; sparse=true)) == II

x = ["b", "a", "c", "c", "b"]
@test indicatormat(x) == II
@test Matrix(indicatormat(x; sparse=true)) == II

io = IOBuffer()
describe(io, collect(1:10))
@test String(take!(io)) == """
                           Summary Stats:
                           Mean:           5.500000
                           Minimum:        1.000000
                           1st Quartile:   3.250000
                           Median:         5.500000
                           3rd Quartile:   7.750000
                           Maximum:        10.000000
                           Length:         10
                           Type:           $Int
                           """

describe(io, fill("s", 3))
@test String(take!(io)) == """
                           Summary Stats:
                           Length:         3
                           Type:           String
                           Number Unique:  1
                           """
