using StatsBase
using LinearAlgebra, Random, Test

@testset "StatsBase.Reliability" begin
    # basic vanilla test
    cov_X = [10 6 6 6;
             6 11 6 6;
             6 6 12 6;
             6 6 6 13]
    reliability_X = crombach_alpha(cov_X)
    @test reliability_X.alpha ≈ 0.8135593220338981
    @test reliability_X.dropped[1].second ≈ 0.75
    @test reliability_X.dropped[2].second ≈ 0.7605633802816901
    @test reliability_X.dropped[3].second ≈ 0.7714285714285715
    @test reliability_X.dropped[4].second ≈ 0.782608695652174

    # testing Rational
    cov_rational = 1 .// cov_X
    reliability_rational = crombach_alpha(cov_X)
    @test reliability_rational.alpha ≈ 0.8135593220338981
    @test reliability_rational.dropped[1].second ≈ 0.75
    @test reliability_rational.dropped[2].second ≈ 0.7605633802816901
    @test reliability_rational.dropped[3].second ≈ 0.7714285714285715
    @test reliability_rational.dropped[4].second ≈ 0.782608695652174

    # testing BigFloat
    cov_bigfloat = BigFloat.(cov_X)
    reliability_bigfloat = crombach_alpha(cov_X)
    @test reliability_bigfloat.alpha ≈ 0.8135593220338981
    @test reliability_bigfloat.dropped[1].second ≈ 0.75
    @test reliability_bigfloat.dropped[2].second ≈ 0.7605633802816901
    @test reliability_bigfloat.dropped[3].second ≈ 0.7714285714285715
    @test reliability_bigfloat.dropped[4].second ≈ 0.782608695652174

    # testing corner cases
    @test_throws MethodError crombach_alpha([1.0, 2.0])
    cov_k2 = [10 6;
              6 11]
    reliability_k2 = crombach_alpha(cov_k2)
    @test reliability_k2.alpha ≈ 0.7272727272727273
    @test isempty(reliability_k2.dropped)
end # @testset "StatsBase.Reliability"
