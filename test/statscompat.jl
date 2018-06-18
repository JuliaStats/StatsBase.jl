using StatsBase, Test

@testset "StatsCompat for replacing some Compat functions" begin
    @test StatsCompat.varm([1 2; 3 4], -1) == 18
    @test StatsCompat.varm([1 2; 3 4], [-1 -2], dims=1) == [20 52]
    @test StatsCompat.varm([1 2; 3 4], [-1, -2], dims=2) == hcat([13, 61])
    @test StatsCompat.var([1 2; 3 4]) == 5/3
    @test StatsCompat.var([1 2; 3 4], dims=1) == [2 2]
    @test StatsCompat.var([1 2; 3 4], dims=2) == hcat([0.5, 0.5])
    @test StatsCompat.var([1 2; 3 4], corrected=false) == 1.25
    @test StatsCompat.var([1 2; 3 4], corrected=false, dims=1) == [1 1]
    @test StatsCompat.var([1 2; 3 4], corrected=false, dims=2) == hcat([0.25, 0.25])
    @test StatsCompat.std([1 2; 3 4]) == sqrt(5/3)
    @test StatsCompat.std([1 2; 3 4], dims=1) == [sqrt(2) sqrt(2)]
    @test StatsCompat.std([1 2; 3 4], dims=2) == hcat([sqrt(0.5), sqrt(0.5)])
    @test StatsCompat.std([1 2; 3 4], corrected=false) == sqrt(1.25)
    @test StatsCompat.std([1 2; 3 4], corrected=false, dims=1) == [sqrt(1) sqrt(1)]
    @test StatsCompat.std([1 2; 3 4], corrected=false, dims=2) == hcat([sqrt(0.25), sqrt(0.25)])
    @test StatsCompat.cov([1 2; 3 4]) == [2 2; 2 2]
    @test StatsCompat.cov([1 2; 3 4], dims=1) == [2 2; 2 2]
    @test StatsCompat.cov([1 2; 3 4], dims=2) == [0.5 0.5; 0.5 0.5]
    @test StatsCompat.cov([1 2; 3 4], [4; 5]) == hcat([1; 1])
    @test StatsCompat.cov([1 2; 3 4], [4; 5], dims=1) == hcat([1; 1])
    @test StatsCompat.cov([1 2; 3 4], [4; 5], dims=2) == hcat([0.5; 0.5])
    @test StatsCompat.cov([1 2; 3 4], [4; 5], corrected=false) == hcat([0.5; 0.5])
    @test StatsCompat.cov([1 2; 3 4], [4; 5], corrected=false, dims=1) == hcat([0.5; 0.5])
    @test StatsCompat.cov([1 2; 3 4], [4; 5], corrected=false, dims=2) == hcat([0.25; 0.25])
    @test StatsCompat.cor([1 2; 3 4]) ≈ [1 1; 1 1]
    @test StatsCompat.cor([1 2; 3 4], dims=1) ≈ [1 1; 1 1]
    @test StatsCompat.cor([1 2; 3 4], dims=2) ≈ [1 1; 1 1]
    @test StatsCompat.cor([1 2; 3 4], [4; 5]) ≈ [1; 1]
    @test StatsCompat.cor([1 2; 3 4], [4; 5], dims=1) ≈ [1; 1]
    @test StatsCompat.cor([1 2; 3 4], [4; 5], dims=2) ≈ [1; 1]
end
