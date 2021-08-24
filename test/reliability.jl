using StatsBase
using LinearAlgebra, Random, Test

@testset "Cronbach's Alpha" begin
    # basic vanilla test
    cov_X = [10 6 6 6;
             6 11 6 6;
             6 6 12 6;
             6 6 6 13]
    cronbach_X = cronbachalpha(cov_X)
    @test cronbach_X isa CronbachAlpha{Float64}
    @test cronbach_X.alpha ≈ 0.8135593220338981
    @test cronbach_X.dropped ≈
        [0.75, 0.7605633802816901, 0.7714285714285715, 0.782608695652174]

    # testing Rational
    cov_rational = cov_X .// 1
    cronbach_rational = cronbachalpha(cov_rational)
    @test cronbach_rational isa CronbachAlpha{Rational{Int}}
    @test cronbach_rational.alpha == 48 // 59
    @test cronbach_rational.dropped ==
        [3 // 4, 54 // 71, 27 // 35, 18 // 23]

    # testing BigFloat
    cov_bigfloat = BigFloat.(cov_X)
    cronbach_bigfloat = cronbachalpha(cov_bigfloat)
    @test cronbach_bigfloat isa CronbachAlpha{BigFloat}
    @test cronbach_bigfloat.alpha ≈ 0.8135593220338981
    @test cronbach_bigfloat.dropped ≈
        [0.75, 0.7605633802816901, 0.7714285714285715, 0.782608695652174]

    # testing corner cases
    @test_throws MethodError cronbachalpha([1.0, 2.0])
    cov_k2 = [10 6;
              6 11]
    cronbach_k2 = cronbachalpha(cov_k2)
    @test cronbach_k2.alpha ≈ 0.7272727272727273
    @test isempty(cronbach_k2.dropped)

    # testing when Matrix is not positive-definite
    cov_not_pos = [-1 1;
                   -1 1]
    @test_throws ArgumentError cronbachalpha(cov_not_pos)

    # testing with a zero
    cov_zero = [1 2;
                0 1]
    @test_throws ArgumentError cronbachalpha(cov_not_pos)

    # testing with one column
    cov_k1 = reshape([1, 2], 2, 1)
    @test_throws ArgumentError cronbachalpha(cov_k1)

    # testing with Missing
    cov_missing = [1 2;
                   missing 1]
    @test_throws MethodError cronbachalpha(cov_missing)


    # testing Base.show
    cronbach_X = cronbachalpha(cov_X)
    io = IOBuffer()
    show(io, cronbach_X)
    str = String(take!(io))
    @test str == """
        Cronbach's alpha for all items: 0.8136

        Cronbach's alpha if an item is dropped:
        item 1: 0.7500
        item 2: 0.7606
        item 3: 0.7714
        item 4: 0.7826
        """
    # for two columns
    io = IOBuffer()
    show(io, cronbach_k2)
    str = String(take!(io))
    @test str == "Cronbach's alpha for all items: 0.7273\n"

end # @testset "Cronbach's Alpha"
