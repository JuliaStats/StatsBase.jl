using StatsBase
using LinearAlgebra, Random, Test

@testset "StatsBase.Reliability" begin
    cov_X = [10 6 6 6;
             6 11 6 6;
             6 6 12 6;
             6 6 6 13]
    @test StatsBase._crombach_alpha(cov_X) ≈ 0.8135593220338984
    reliabity_x = crombach_alpha(cov_X)
    @test reliabity_x.dropped[1].second ≈ 0.75
    @test reliabity_x.dropped[2].second ≈ 0.7605633802816901
    @test reliabity_x.dropped[3].second ≈ 0.7714285714285715
    @test reliabity_x.dropped[4].second ≈ 0.782608695652174
end # @testset "StatsBase.Reliability"
