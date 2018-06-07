using StatsBase
using Compat
using Compat.Test

wechsler = Float32[
     7  5  9  8
     8  8  5  6
    16 18 11  9
     8  3  7  9
     6  3 13  9
    11  8 10 10
    12  7  9  8
     8 11  9  3
    14 12 11  4
    13 13 13  6
    13  9  9  9
    13 10 15  7
    14 11 12  8
    15 11 11 10
    13 10 15  9
    10  5  8  6
    10  3  7  7
    17 13 13  7
    10  6 10  7
    10 10 15  8
    14  7 11  5
    16 11 12 11
    10  7 14  6
    10 10  9  6
    10  7 10 10
     7  6  5  9
    15 12 10  6
    17 15 15  8
    16 13 16  9
    13 10 17  8
    13 10 17 10
    19 12 16 10
    19 15 17 11
    13 10  7  8
    15 11 12  8
    16  9 11 11
    14 13 14  9
]

@test @inferred(partialcor(wechsler[:,1], wechsler[:,2], wechsler[:,3:4])) ≈ 0.7118787 rtol=1e-6

X = [ 2 1 0
      4 2 0
     15 3 1
     20 4 1]
@test @inferred(partialcor(view(X,:,1), view(X,:,2), view(X,:,3))) ≈ 0.919145 rtol=1e-6
