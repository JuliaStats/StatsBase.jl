using StatsBase
using Base.Test

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

# findat

findat([2, 4, 6, 8, 10], [8, 6, 4]) == [4, 3, 2]
findat([2, 4, 6, 8, 10], [2, 3, 4]) == [1, 0, 2]

# indicatormat

I = map(Bool, 
	[0 1 0 0 0; 
     1 0 0 0 1;
     0 0 1 1 0])

x = [2, 1, 3, 3, 2]
@test indicatormat(x, 3) == I
@test full(indicatormat(x, 3; sparse=true)) == I

x = ["b", "a", "c", "c", "b"]
@test indicatormat(x) == I
@test full(indicatormat(x; sparse=true)) == I
