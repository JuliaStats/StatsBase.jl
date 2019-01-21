using StatsBase
import StatsBase: transform, reconstruct
using Statistics
using Test

@testset "Transformations" begin
    X = rand(5, 8)

    t = fit(ZScoreTransform, X; center=false, scale=false)
    Y = transform(t, X)
    @test isa(t, AbstractDataTransform)
    @test isempty(t.mean)
    @test isempty(t.scale)
    @test isequal(X, Y)
    @test transform(t, X[:,1]) ≈ Y[:,1]
    @test reconstruct(t, Y[:,1]) ≈ X[:,1]
    @test reconstruct(t, Y) ≈ X

    t = fit(ZScoreTransform, X; center=false)
    Y = transform(t, X)
    @test isempty(t.mean)
    @test length(t.scale) == 5
    @test Y ≈ X ./ std(X, dims=2)
    @test transform(t, X[:,1]) ≈ Y[:,1]
    @test reconstruct(t, Y[:,1]) ≈ X[:,1]
    @test reconstruct(t, Y) ≈ X

    t = fit(ZScoreTransform, X; scale=false)
    Y = transform(t, X)
    @test length(t.mean) == 5
    @test isempty(t.scale)
    @test Y ≈ X .- mean(X, dims=2)
    @test transform(t, X[:,1]) ≈ Y[:,1]
    @test reconstruct(t, Y[:,1]) ≈ X[:,1]
    @test reconstruct(t, Y) ≈ X

    t = fit(ZScoreTransform, X)
    Y = transform(t, X)
    @test length(t.mean) == 5
    @test length(t.scale) == 5
    @test Y ≈ (X .- mean(X, dims=2)) ./ std(X, dims=2)
    @test transform(t, X[:,1]) ≈ Y[:,1]
    @test reconstruct(t, Y[:,1]) ≈ X[:,1]
    @test reconstruct(t, Y) ≈ X

    t = fit(UnitRangeTransform, X)
    Y = transform(t, X)
    @test isa(t, AbstractDataTransform)
    @test length(t.min) == 5
    @test length(t.scale) == 5
    @test Y ≈ (X .- minimum(X, dims=2)) ./ (maximum(X, dims=2) .- minimum(X, dims=2))
    @test transform(t, X[:,1]) ≈ Y[:,1]
    @test reconstruct(t, Y[:,1]) ≈ X[:,1]
    @test reconstruct(t, Y) ≈ X

    t = fit(UnitRangeTransform, X; unit=false)
    Y = transform(t, X)
    @test length(t.min) == 5
    @test length(t.scale) == 5
    @test Y ≈ X ./ (maximum(X, dims=2) .- minimum(X, dims=2))
    @test transform(t, X[:,1]) ≈ Y[:,1]
    @test reconstruct(t, Y[:,1]) ≈ X[:,1]
    @test reconstruct(t, Y) ≈ X

    @test Y == standardize(UnitRangeTransform, X; unit=false)
end
