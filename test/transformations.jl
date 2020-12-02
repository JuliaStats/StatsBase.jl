using StatsBase
import StatsBase: transform, reconstruct, transform!, reconstruct!
using Statistics
using Test
cuda_available = try
    using StatsBase.CUDA
    true
catch
    @warn "CUDA not available - transformations not tested on GPU"
    false
end

devices = [:cpu]
if cuda_available
    if CUDA.functional(true)
        push!(devices, :gpu)
    else
        @warn "GPU not functional - transformations not tested on GPU"
    end
end

@testset "Transformations on $device" for device in devices
    # matrix
    X = rand(5, 8)
    if device == :gpu
        X = cu(X)
    end
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
    X = rand(10)
    if device == :gpu
        X = cu(X)
    end

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
