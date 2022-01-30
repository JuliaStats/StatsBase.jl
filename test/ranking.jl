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

@testset "quantilerank" begin
    @testset "value as number and array" begin
        @testset ":inc and :exc" begin
            v1 = [1, 1, 1, 2, 3, 4, 8, 11, 12, 13]
            v2 = [1, 2, 3, 6, 6, 6, 7, 8, 9]
            v3 = [1, 2, 4, 3, 4]
            v4 = [1, 2, 1, 3, 4]
            @test quantilerank(v1, 2, method=:inc)    == 1/3
            @test quantilerank(v1, 4, method=:inc)    == 5/9
            @test quantilerank(v1, 8, method=:inc)    == 2/3
            @test quantilerank(v1, 5, method=:inc)    == 7/12        
            @test quantilerank(v2, 7, method=:exc)    == 0.7
            @test quantilerank(v2, 5.43, method=:exc) == 0.381
            @test quantilerank(v3, 4, method=:exc)    == 6/9
            @test quantilerank(v3, 4, method=:inc)    == 3/4
            @test quantilerank(v4, 1, method=:exc)    == 1/6
            @test quantilerank(v4, -100, method=:inc) == 0.0
            @test quantilerank(v4,  100, method=:inc) == 1.0
            @test quantilerank(v4, -100, method=:exc) == 0.0
            @test quantilerank(v4,  100, method=:exc) == 1.0
        end
        @testset ":compete" begin
            v = [0, 0, 1, 1, 2, 2, 2, 2, 4, 4]
            @test quantilerank(v, 1, method=:compete) == 2/9
            @test quantilerank(v, 2, method=:compete) == 4/9
            @test quantilerank(v, 4, method=:compete) == 8/9
            @test quantilerank(v, -100, method=:compete) == 0.0
            @test quantilerank(v,  100, method=:compete) == 1.0
        end
        @testset ":strict, :weak and :tied" begin
            v = [7, 8, 2, 1, 3, 4, 5, 4, 6, 9]
            for (method, res1, res2) in [(:tied, .4, [.4, .85]),
                                         (:strict, .3, [.3, .8]),
                                         (:weak, .5, [.5, .9])]
                @test quantilerank(v, 4, method=method) == res1
            end
        end
    end
    @testset "errors" begin
        v1 = [1, 2, 3, 5, 6, missing, 8]
        v2 = [missing, missing]
        v3 = [1, 2, 3, 5, 6, NaN, 8]
        v4 = [1, 2, 3, 3, 4]
        for method in (:tied, :strict, :weak)
            @test_throws ArgumentError quantilerank(v1, 4, method=method)
            @test_throws ArgumentError quantilerank(v2, 4, method=method)
            @test_throws ArgumentError quantilerank(v3,  4, method=method)
        end
        @test_throws ArgumentError quantilerank(v4, 3, method=:wrongargument)
        @test_throws ArgumentError quantilerank(v4,  NaN)
        @test_throws ArgumentError quantilerank(v4,  missing)
        @test_throws ArgumentError quantilerank([], 3)
        @test_throws ArgumentError quantilerank([1], 3)
    end
end

@testset "percentilerank" begin
    v1 = [1, 1, 1, 2, 3, 4, 8, 11, 12, 13]
    v2 = [1, 2, 3, 6, 6, 6, 7, 8, 9]
    v3 = [0, 0, 1, 1, 2, 2, 2, 2, 4, 4]
    # :inc and :exc
    @test percentilerank(v1,  2, method=:inc) == 100 * quantilerank(v1, 2, method=:inc)
    @test percentilerank(v2, 7, method=:exc)     == 100 * quantilerank(v2, 7, method=:exc)
    @test percentilerank(v2, 5.43, method=:exc)  == 100 * quantilerank(v2, 5.43, method=:exc)
    
    # :compete
    @test percentilerank(v3,  1, method=:compete) == 100 * quantilerank(v3, 1, method=:compete)
    
    @testset ":strict, :weak and :tied" begin
        v = [7, 8, 2, 1, 3, 4, 5, 4, 6, 9]
        for (method, res1) in [(:tied, .4),
                               (:strict, .3),
                               (:weak, .5)]
            @test percentilerank(v,  4, method=method) == res1 * 100
        end
    end

    @testset "errors" begin
        v1 = [1, 2, 3, 5, 6, missing, 8]
        v2 = [missing, missing]
        v3 = [1, 2, 3, 3, 4]
        for method in (:tied, :strict, :weak)
            @test_throws ArgumentError percentilerank(v1,  4, method=method)
            @test_throws ArgumentError percentilerank(v2,  4, method=method)
        end
        @test_throws ArgumentError percentilerank(v3,  3, method=:wrongargument)
    end
end
