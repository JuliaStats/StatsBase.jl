using StatsBase
using Test

a = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0]
x = [3.0, 1.0, 2.0, 4.0, 4.0, 2.0, 5.0, 4.0]  # x is a permutated version of a
xm = [3.0, 1.0, 2.0, 4.0, 4.0, 2.0, 5.0, 4.0, missing]
s = ["c", "a", "b", "d", "d", "b", "e", "d"] # s is a vector of strings ordered like x

@test ordinalrank(a) == [1, 2, 3, 4, 5, 6, 7, 8]
@test ordinalrank(x) == [4, 1, 2, 5, 6, 3, 8, 7]
@test isequal(ordinalrank(xm), [4, 1, 2, 5, 6, 3, 8, 7, missing])
@test isequal(ordinalrank([missing, missing]), [missing, missing])
@test ordinalrank(s) == ordinalrank(x)
@test ordinalrank(x, rev = true) == ordinalrank(-x)
@test ordinalrank(x, lt = (x, y) -> isless(y, x)) == ordinalrank(-x)

@test competerank(a) == [1, 2, 2, 4, 5, 5, 5, 8]
@test competerank(x) == [4, 1, 2, 5, 5, 2, 8, 5]
@test isequal(competerank(xm), [4, 1, 2, 5, 5, 2, 8, 5, missing])
@test isequal(competerank([missing, missing]), [missing, missing])
@test competerank(s) == competerank(x)
@test competerank(x, rev = true) == competerank(-x)
@test competerank(x, lt = (x, y) -> isless(y, x)) == competerank(-x)

@test denserank(a) == [1, 2, 2, 3, 4, 4, 4, 5]
@test denserank(x) == [3, 1, 2, 4, 4, 2, 5, 4]
@test isequal(denserank(xm), [3, 1, 2, 4, 4, 2, 5, 4, missing])
@test isequal(denserank([missing, missing]), [missing, missing])
@test denserank(s) == denserank(x)
@test denserank(x, rev = true) == denserank(-x)
@test denserank(x, lt = (x, y) -> isless(y, x)) == denserank(-x)

@test tiedrank(a) == [1.0, 2.5, 2.5, 4.0, 6.0, 6.0, 6.0, 8.0]
@test tiedrank(x) == [4.0, 1.0, 2.5, 6.0, 6.0, 2.5, 8.0, 6.0]
@test isequal(tiedrank(xm), [4.0, 1.0, 2.5, 6.0, 6.0, 2.5, 8.0, 6.0, missing])
@test isequal(tiedrank([missing, missing]), [missing, missing])
@test tiedrank(s) == tiedrank(x)
@test tiedrank(x, rev = true) == tiedrank(-x)
@test tiedrank(x, lt = (x, y) -> isless(y, x)) == tiedrank(-x)

@testset "quantilerank and percentilerank" begin

    pcrk = percentilerank
    qtrk = quantilerank
 
    @testset "score as number and array" begin        
        for (method, res1, res2) in [(:rank, 40.0, [40.0, 85.0]), 
                                    (:strict, 30.0, [30.0, 80.0]),
                                    (:weak, 50.0, [50.0, 90.0]), 
                                    (:mean, 40.0, [40.0, 85.0])]
            @test pcrk([1,2,3,4,4,5,6,7,8,9], 4, method=method) == res1
            @test qtrk([1,2,3,4,4,5,6,7,8,9], 4, method=method) == res1 / 100
            @test pcrk([1,2,3,4,4,5,6,7,8,9], [4, 8], method=method) == res2
            @test qtrk([1,2,3,4,4,5,6,7,8,9], [4, 8], method=method) == res2 / 100
        end
    end
 
    @testset "missing and error" begin
        for method in (:rank, :mean, :strict, :weak)
            @test pcrk([1,2,3,5,6,7,8,9,10,11], 4, method=method) == 30
            @test qtrk([1,2,3,5,6,7,8,9,10,11], 4, method=method) == 0.3
        end

        @test_throws ArgumentError pcrk([1, 2, 3, 3, 4], 3, method=:wrongargument)
        @test_throws ArgumentError qtrk([1, 2, 3, 3, 4], 3, method=:wrongargument)
    end
end
