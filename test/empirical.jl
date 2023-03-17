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
    @test isnan(ecdf([1,2,3])(NaN))
    @test_throws ArgumentError ecdf([1, NaN])
end

@testset "Weighted ECDF" begin
    x = randn(10000000)
    w1 = rand(10000000)
    w2 = weights(w1)
    fnecdf = ecdf(x, weights=w1)
    fnecdfalt = ecdf(x, weights=w2)
    @test fnecdf.sorted_values == fnecdfalt.sorted_values
    @test fnecdf.weights == fnecdfalt.weights
    @test fnecdf.weights != w1  #  check that w wasn't accidentally modified in place
    @test fnecdfalt.weights != w2
    y = [-1.96, -1.644854, -1.281552, -0.6744898, 0, 0.6744898, 1.281552, 1.644854, 1.96]
    @test isapprox(fnecdf(y), [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975], atol=1e-3)
    @test isapprox(fnecdf(1.96), 0.975, atol=1e-3)
    @test fnecdf(y) ≈ map(fnecdf, y)
    @test extrema(fnecdf) == (minimum(fnecdf), maximum(fnecdf)) == extrema(x)
    fnecdf = ecdf([1.0, 0.5], weights=weights([3, 1]))
    @test fnecdf(0.75) == 0.25
    @test extrema(fnecdf) == (minimum(fnecdf), maximum(fnecdf)) == (0.5, 1.0)
    @test_throws ArgumentError ecdf(rand(8), weights=weights(rand(10)))
    #  Check frequency weights
    v = randn(100)
    r = rand(1:100, 100)
    vv = vcat(fill.(v, r)...)  #  repeat elements of v according to r
    fw = fweights(r)
    frecdf1 = ecdf(v, weights=fw)
    frecdf2 = ecdf(vv)
    @test frecdf1(y) ≈ frecdf2(y)
    #  Check probability weights
    a = randn(100)
    b = rand(100)
    b̃ = abs(10randn()) * b
    bw1 = pweights(b)
    bw2 = pweights(b̃)
    precdf1 = ecdf(a, weights=bw1)
    precdf2 = ecdf(a, weights=bw2)
    @test precdf1(y) ≈ precdf2(y)
end
