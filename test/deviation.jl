using StatsBase
using Base.Test

a = [1, 2, 3, 4, 5, 6, 7]
b = [1, 3, 3, 4, 6, 7, 8]

@test counteq(a, b) == 3
@test countne(a, b) == 4

a = rand(5, 6)
b = rand(5, 6)

@test_approx_eq sqL2dist(a, b) sum(abs2(a - b))
@test_approx_eq L2dist(a, b)   sqrt(sqL2dist(a, b))
@test_approx_eq L1dist(a, b)   sum(abs(a - b))
@test_approx_eq Linfdist(a, b) maximum(abs(a - b))

@test_approx_eq gkldiv(a, b) sum(a .* log(a ./ b) - a + b)

@test_approx_eq meanad(a, b) mean(abs(a - b))
@test_approx_eq maxad(a, b)  maximum(abs(a - b))
@test_approx_eq msd(a, b)    mean(abs2(a - b))
@test_approx_eq rmsd(a, b)   sqrt(msd(a, b))
@test_approx_eq rmsd(a, b; normalize=true)  rmsd(a, b) / (maximum(a) - minimum(a))
@test_approx_eq psnr(a, b, 2)  10 * log10(4 / msd(a, b))   

