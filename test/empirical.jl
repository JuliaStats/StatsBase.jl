using StatsBase
using Compat
using Compat.Test

X = randn(10000000)
fnecdf = ecdf(X)
@test isapprox(fnecdf([-1.96, -1.644854, -1.281552, -0.6744898, 0, 0.6744898, 1.281552, 1.644854, 1.96]),
    [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975], atol=1e-3)
@test isapprox(fnecdf(1.96), 0.975, atol=1e-3)
@test fnecdf([-1.96, -1.644854, -1.281552, -0.6744898, 0, 0.6744898, 1.281552, 1.644854, 1.96]) ≈
    map(fnecdf, [-1.96, -1.644854, -1.281552, -0.6744898, 0, 0.6744898, 1.281552, 1.644854, 1.96])

fnecdf = ecdf([0.5])
@test fnecdf([zeros(5000); ones(5000)]) == [zeros(5000); ones(5000)]

fnepdf = epdf(X)
@test isapprox(fnepdf.([-10.0, 0.0]), [0.0, 0.383]; atol = 2e-3)

fnepdf = epdf(X; closed = :right)
@test isapprox(fnepdf.([-10.0, 0.0]), [0.0, 0.383]; atol = 2e-3)


fnepdf = epdf(X; nbins = 100)
@test isapprox(fnepdf.([-10.0, 0.0]), [0.0, 0.398]; atol = 2e-3)


fnepdf = epdf(X; nbins = 100, closed = :right)
@test isapprox(fnepdf.([-10.0, 0.0]), [0.0, 0.398]; atol = 2e-3)
