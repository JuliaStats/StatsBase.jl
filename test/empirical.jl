using StatsBase
using Test

@testset "ECDF" begin
    x = randn(10000000)
    fnecdf = ecdf(x)
    y = [-1.96, -1.644854, -1.281552, -0.6744898, 0, 0.6744898, 1.281552, 1.644854, 1.96]
    @test isapprox(fnecdf(y), [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975], atol=1e-3)
    @test isapprox(fnecdf(1.96), 0.975, atol=1e-3)
    @test fnecdf(y) ≈ map(fnecdf, y)
    @test extrema(fnecdf) == (minimum(fnecdf), maximum(fnecdf)) == extrema(x)
    fnecdf = ecdf([0.5])
    @test fnecdf([zeros(5000); ones(5000)]) == [zeros(5000); ones(5000)]
    @test extrema(fnecdf) == (minimum(fnecdf), maximum(fnecdf)) == (0.5, 0.5)
end

@testset "Weighted ECDF" begin
    x = randn(10000000)
    w = weights(rand(10000000))
    fnecdf = ecdf(x, w)
    y = [-1.96, -1.644854, -1.281552, -0.6744898, 0, 0.6744898, 1.281552, 1.644854, 1.96]
    @test fnecdf.weights.values != w.values
    @test fnecdf(y) ≈ map(fnecdf, y)
    @test extrema(fnecdf) == (minimum(fnecdf), maximum(fnecdf)) == extrema(x)
    fnecdf = ecdf([1.0, 0.5], weights([3, 1]))
    @test fnecdf(0.75) == 0.25
    @test extrema(fnecdf) == (minimum(fnecdf), maximum(fnecdf)) == (0.5, 1.0)

    @test_throws MethodError ecdf(rand(8), rand(10))
end
