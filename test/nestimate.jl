using StatsBase 
using Test 

@testset "nestimate (MVUE and MOM estimators)" begin 
    @testset "mvue: Minimum-variance unbiased istimator (of N)" begin 
        vect = [19, 40, 42, 60]
        result = nestimate(vect)

        @test result == nestimate(vect, method = "mvue")
        @test result isa Float64
        @test result == 74 
    end 

    @testset "Method of moments (mom) estimate of N" begin 
        vect = [19, 40, 42, 60]
        result = nestimate(vect, method = "mom")

        @test result isa Float64 
        @test result == 79.5
    end  
end 