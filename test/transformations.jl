using StatsBase
import StatsBase: transform, reconstruct, transform!, reconstruct!
using Statistics
using CUDA
using Test

CUDA.allowscalar(false)

@warn "CUDA related method definitions included only in test script"
include("cuda.jl")

@testset "Transformations on $device" for device in (:cpu, :gpu)
    if device==:gpu
        if !CUDA.functional(true)
            @test_broken "GPU not functional - transformations not tested on GPU"
           continue
        end
    end

    to_device = device == :cpu ? identity : cu

    # matrix
    X = rand(5, 8) |> to_device
    X_ = copy(X)

    t = fit(ZScoreTransform, X, dims=1, center=false, scale=false)
    Y = transform(t, X)
    @test isa(t, AbstractDataTransform)
    @test isempty(t.mean)
    @test isempty(t.scale)
    @test X == Y
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(ZScoreTransform, X, dims=1, center=false)
    Y = transform(t, X)
    @test isempty(t.mean)
    @test length(t.scale) == 8
    @test Y ≈ X ./ std(X, dims=1)
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(ZScoreTransform, X, dims=1, scale=false)
    Y = transform(t, X)
    @test length(t.mean) == 8
    @test isempty(t.scale)
    @test Y ≈ X .- mean(X, dims=1)
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(ZScoreTransform, X, dims=1)
    Y = transform(t, X)
    @test length(t.mean) == 8
    @test length(t.scale) == 8
    @test Y ≈ (X .- mean(X, dims=1)) ./ std(X, dims=1)
    @test reconstruct(t, Y) ≈ X
    @test Y ≈ standardize(ZScoreTransform, X, dims=1)
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(ZScoreTransform, X, dims=2)
    Y = transform(t, X)
    @test length(t.mean) == 5
    @test length(t.scale) == 5
    @test Y ≈ (X .- mean(X, dims=2)) ./ std(X, dims=2)
    @test reconstruct(t, Y) ≈ X
    @test Y ≈ standardize(ZScoreTransform, X, dims=2)
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(UnitRangeTransform, X, dims=1, unit=false)
    Y = transform(t, X)
    @test length(t.min) == 8
    @test length(t.scale) == 8
    @test Y ≈ X ./ (maximum(X, dims=1) .- minimum(X, dims=1))
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(UnitRangeTransform, X, dims=1)
    Y = transform(t, X)
    @test isa(t, AbstractDataTransform)
    @test length(t.min) == 8
    @test length(t.scale) == 8
    @test Y ≈ (X .- minimum(X, dims=1)) ./ (maximum(X, dims=1) .- minimum(X, dims=1))
    @test reconstruct(t, Y) ≈ X
    @test Y ≈ standardize(UnitRangeTransform, X, dims=1)
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(UnitRangeTransform, X, dims=2)
    Y = transform(t, X)
    @test isa(t, AbstractDataTransform)
    @test length(t.min) == 5
    @test length(t.scale) == 5
    @test Y ≈ (X .- minimum(X, dims=2)) ./ (maximum(X, dims=2) .- minimum(X, dims=2))
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    # vector
    if device == :gpu
        @test_broken "vector-valued transformations on GPU"
        continue
    end
    X = rand(10) |> to_device
    X_ = copy(X)

    t = fit(ZScoreTransform, X, dims=1, center=false, scale=false)
    Y = transform(t, X)
    @test transform(t, X) ≈ Y
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(ZScoreTransform, X, dims=1, center=false)
    Y = transform(t, X)
    @test Y ≈ X ./ std(X, dims=1)
    @test transform(t, X) ≈ Y
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(ZScoreTransform, X, dims=1, scale=false)
    Y = transform(t, X)
    @test Y ≈ X .- mean(X, dims=1)
    @test transform(t, X) ≈ Y
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(ZScoreTransform, X, dims=1)
    Y = transform(t, X)
    @test Y ≈ (X .- mean(X, dims=1)) ./ std(X, dims=1)
    @test transform(t, X) ≈ Y
    @test reconstruct(t, Y) ≈ X
    @test Y ≈ standardize(ZScoreTransform, X, dims=1)
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(UnitRangeTransform, X, dims=1)
    Y = transform(t, X)
    @test Y ≈ (X .- minimum(X, dims=1)) ./ (maximum(X, dims=1) .- minimum(X, dims=1))
    @test transform(t, X) ≈ Y
    @test reconstruct(t, Y) ≈ X
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

    X = copy(X_)
    t = fit(UnitRangeTransform, X, dims=1, unit=false)
    Y = transform(t, X)
    @test Y ≈ X ./ (maximum(X, dims=1) .- minimum(X, dims=1))
    @test transform(t, X) ≈ Y
    @test reconstruct(t, Y) ≈ X
    @test Y ≈ standardize(UnitRangeTransform, X, dims=1, unit=false)
    @test transform!(t, X) === X
    @test X == Y
    @test reconstruct!(t, Y) === Y
    @test Y ≈ X_

end
