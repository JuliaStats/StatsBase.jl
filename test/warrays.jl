# Tests of weighted arrays

using Stats
using Base.Test

x = rand(10)
w = rand(10)

wv = WeightedVector(x, w)
@test values(wv) === x
@test weights(wv) === w
@test weightsum(wv) == sum(w)

