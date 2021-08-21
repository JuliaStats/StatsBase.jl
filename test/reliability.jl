using StatsBase
using LinearAlgebra, Random, Test

@testset "StatsBase.Reliability" begin
    # basic vanilla test
    cov_X = [10 6 6 6;
             6 11 6 6;
             6 6 12 6;
             6 6 6 13]
    reliability_X = crombachalpha(cov_X)
    @test reliability_X isa Reliability{Float64}
    @test reliability_X.alpha ≈ 0.8135593220338981
    @test reliability_X.dropped[1] ≈ 0.75
    @test reliability_X.dropped[2] ≈ 0.7605633802816901
    @test reliability_X.dropped[3] ≈ 0.7714285714285715
    @test reliability_X.dropped[4] ≈ 0.782608695652174

    # testing Rational
    cov_rational = cov_X .// 1
    reliability_rational = crombachalpha(cov_rational)
    @test reliability_rational isa Reliability{Rational{Int}}
    @test reliability_rational.alpha == 48 // 59
    @test reliability_rational.dropped[1] == 3 // 4
    @test reliability_rational.dropped[2] == 54 // 71
    @test reliability_rational.dropped[3] == 27 // 35
    @test reliability_rational.dropped[4] == 18 // 23

    # testing BigFloat
    cov_bigfloat = BigFloat.(cov_X)
    reliability_bigfloat = crombachalpha(cov_bigfloat)
    @test reliability_bigfloat isa Reliability{BigFloat}
    @test reliability_bigfloat.alpha ≈ 0.8135593220338981
    @test reliability_bigfloat.dropped[1] ≈ 0.75
    @test reliability_bigfloat.dropped[2] ≈ 0.7605633802816901
    @test reliability_bigfloat.dropped[3] ≈ 0.7714285714285715
    @test reliability_bigfloat.dropped[4] ≈ 0.782608695652174

    # testing corner cases
    @test_throws MethodError crombachalpha([1.0, 2.0])
    cov_k2 = [10 6;
              6 11]
    reliability_k2 = crombachalpha(cov_k2)
    @test reliability_k2.alpha ≈ 0.7272727272727273
    @test isempty(reliability_k2.dropped)

    # testing when Matrix is not positive-definite
    cov_not_pos = [-1 1;
                   -1 1]
    @test_throws ArgumentError crombachalpha(cov_not_pos)

    # testing with a zero
    cov_zero = [1 2;
                0 1]
    @test_throws ArgumentError crombachalpha(cov_not_pos)

    # testing with one column
    cov_k1 = reshape([1, 2], 2, 1)
    @test_throws ArgumentError crombachalpha(cov_k1)

    # testing with Missing
    cov_missing = [1 2;
                   missing 1]
    @test_throws MethodError crombachalpha(cov_missing)


    # testing Base.show
    reliability_X = crombachalpha(cov_X)
    io = IOBuffer()
    show(io, reliability_X)
    str = String(take!(io))
    @test str == """
        Reliability for all items: 0.8136
        Reliability if an item is dropped:
        item 1: 0.7500
        item 2: 0.7606
        item 3: 0.7714
        item 4: 0.7826
        """
end # @testset "StatsBase.Reliability"
