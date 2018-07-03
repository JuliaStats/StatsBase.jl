using StatsBase
using Test

a = [1, 2, 3, 4, 5, 6, 7]
b = [1, 3, 3, 4, 6, 7, 8]

@test counteq(a, b) == 3
@test countne(a, b) == 4

a = rand(5, 6)
b = rand(5, 6)

@test sqL2dist(a, b) ≈ sum(abs2.(a - b))
@test L2dist(a, b)   ≈ sqrt(sqL2dist(a, b))
@test L1dist(a, b)   ≈ sum(abs.(a - b))
@test Linfdist(a, b) ≈ maximum(abs.(a - b))

@test gkldiv(a, b) ≈ sum(a .* log.(a ./ b) - a + b)

@test meanad(a, b)               ≈ mean(abs.(a - b))
@test maxad(a, b)                ≈ maximum(abs.(a - b))
@test msd(a, b)                  ≈ mean(abs2.(a - b))
@test rmsd(a, b)                 ≈ sqrt(msd(a, b))
@test rmsd(a, b; normalize=true) ≈ rmsd(a, b) / (maximum(a) - minimum(a))
@test psnr(a, b, 2)              ≈ 10 * log10(4 / msd(a, b))

