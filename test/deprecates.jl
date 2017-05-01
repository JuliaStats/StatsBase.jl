using StatsBase
using Base.Test
using Compat
import Compat: view

@testset "StatsBase.Deprecates" begin

@testset "Deprecates WeightVec and weights" begin
    @test isa(weights([1, 2, 3]), WeightVec{Int})
    @test isa(weights([1., 2., 3.]), WeightVec{Float64})
    @test isa(weights([1 2 3; 4 5 6]), WeightVec{Int})

    @test isa(WeightVec([1, 2, 3], 6), WeightVec{Int})

    @test isempty(weights(Float64[]))
    @test size(weights([1, 2, 3])) == (3,)

    w  = [1., 2., 3.]
    wv = weights(w)
    @test eltype(wv) === Float64
    @test length(wv) === 3
    @test values(wv) === w
    @test sum(wv) === 6.0
    @test !isempty(wv)

    b  = trues(3)
    bv = weights(b)
    @test eltype(bv) === Bool
    @test length(bv) === 3
    @test values(bv) === b
    @test sum(bv)    === 3
    @test !isempty(bv)

    ba = BitArray([true, false, true])
    sa = sparsevec([1., 0., 2.])

    @test sum(ba, wv) === 4.0
    @test sum(sa, wv) === 7.0
end

@testset "Moments" begin
    @testset "Vectors" begin
        x = rand(10)
        wv = weights(rand(10))
        m = mean(x, wv)

        @testset "var" begin
            @test var(x, wv)           ≈ sum(abs2.(x .- m), wv) ./ sum(wv)
            @test var(x, wv; mean=0)   ≈ sum(abs2.(x), wv) ./ sum(wv)
            @test var(x, wv; mean=1.0) ≈ sum(abs2.(x .- 1.0), wv) ./ sum(wv)
        end

        @testset "std" begin
            @test std(x, wv)           ≈ sqrt(var(x, wv))
            @test std(x, wv; mean=0)   ≈ sqrt(var(x, wv; mean=0))
            @test std(x, wv; mean=1.0) ≈ sqrt(var(x, wv; mean=1.0))
        end

        @testset "mean_and_var" begin
            (m, v) = mean_and_var(x)
            @test m == mean(x)
            @test v == var(x)

            (m, v) = mean_and_var(x, wv)
            @test m == mean(x, wv)
            @test v == var(x, wv)
        end

        @testset "mean_and_std" begin
            (m, s) = mean_and_std(x)
            @test m == mean(x)
            @test s == std(x)

            (m, s) = mean_and_std(x, wv)
            @test m == mean(x, wv)
            @test s == std(x, wv)
        end
    end

    @testset "Matrices" begin
        x = rand(5, 6)
        w1 = rand(5)
        w2 = rand(6)
        wv1 = weights(w1)
        wv2 = weights(w2)
        m1 = mean(x, wv1, 1)
        m2 = mean(x, wv2, 2)

        @testset "var" begin
            @test var(x, wv1, 1; mean=0) ≈ sum(abs2.(x) .* w1, 1) ./ sum(wv1)
            @test var(x, wv2, 2; mean=0) ≈ sum(abs2.(x) .* w2', 2) ./ sum(wv2)
            @test var(x, wv1, 1; mean=m1) ≈ sum(abs2.(x .- m1) .* w1, 1) ./ sum(wv1)
            @test var(x, wv2, 2; mean=m2) ≈ sum(abs2.(x .- m2) .* w2', 2) ./ sum(wv2)
            @test var(x, wv1, 1) ≈ sum(abs2.(x .- m1) .* w1, 1) ./ sum(wv1)
            @test var(x, wv2, 2) ≈ sum(abs2.(x .- m2) .* w2', 2) ./ sum(wv2)
        end

        @testset "std" begin
            @test std(x, wv1, 1)          ≈ sqrt.(var(x, wv1, 1))
            @test std(x, wv2, 2)          ≈ sqrt.(var(x, wv2, 2))
            @test std(x, wv1, 1; mean=0)  ≈ sqrt.(var(x, wv1, 1; mean=0))
            @test std(x, wv2, 2; mean=0)  ≈ sqrt.(var(x, wv2, 2; mean=0))
            @test std(x, wv1, 1; mean=m1) ≈ sqrt.(var(x, wv1, 1; mean=m1))
            @test std(x, wv2, 2; mean=m2) ≈ sqrt.(var(x, wv2, 2; mean=m2))
        end

        @testset "mean_and_var" begin
            for d in 1:2
                (m, v) = mean_and_var(x, d)
                @test m == mean(x, d)
                @test v == var(x, d)
            end

            (m, v) = mean_and_var(x, wv1, 1)
            @test m == mean(x, wv1, 1)
            @test v == var(x, wv1, 1)

            (m, v) = mean_and_var(x, wv2, 2)
            @test m == mean(x, wv2, 2)
            @test v == var(x, wv2, 2)
        end

        @testset "mean_and_std" begin
            for d in 1:2
                (m, s) = mean_and_std(x, d)
                @test m == mean(x, d)
                @test s == std(x, d)
            end

            (m, s) = mean_and_std(x, wv1, 1)
            @test m == mean(x, wv1, 1)
            @test s == std(x, wv1, 1)

            (m, s) = mean_and_std(x, wv2, 2)
            @test m == mean(x, wv2, 2)
            @test s == std(x, wv2, 2)
        end
    end
end

@testset "Covariance" begin
    X = randn(3, 8)

    Z1 = X .- mean(X, 1)
    Z2 = X .- mean(X, 2)

    w1 = rand(3)
    w2 = rand(8)

    wv1 = weights(w1)
    wv2 = weights(w2)

    Z1w = X .- mean(X, wv1, 1)
    Z2w = X .- mean(X, wv2, 2)

    S1 = Z1'Z1
    S2 = Z2 * Z2'

    Sz1 = X'X
    Sz2 = X * X'

    S1w = Z1w' * diagm(w1) * Z1w
    S2w = Z2w * diagm(w2) * Z2w'

    Sz1w = X' * diagm(w1) * X
    Sz2w = X * diagm(w2) * X'

    @testset "cov" begin
        @test cov(X, wv1)    ≈ S1w ./ sum(wv1)
        @test cov(X, wv2, 2) ≈ S2w ./ sum(wv2)

        @test Base.covm(X, 0, wv1)    ≈ Sz1w ./ sum(wv1)
        @test Base.covm(X, 0, wv2, 2) ≈ Sz2w ./ sum(wv2)

        @test Base.covm(X, mean(X, wv1, 1), wv1)    ≈ S1w ./ sum(wv1)
        @test Base.covm(X, mean(X, wv2, 2), wv2, 2) ≈ S2w ./ sum(wv2)

        @test Base.covm(X, zeros(1,8), wv1)  ≈ Sz1w ./ sum(wv1)
        @test Base.covm(X, zeros(3), wv2, 2) ≈ Sz2w ./ sum(wv2)
    end

    @testset "mean_and_cov" begin
        (m, C) = mean_and_cov(X, 1)
        @test m == mean(X, 1)
        @test C == cov(X, 1)

        (m, C) = mean_and_cov(X, 2)
        @test m == mean(X, 2)
        @test C == cov(X, 2)

        (m, C) = mean_and_cov(X, wv1)
        @test m == mean(X, wv1, 1)
        @test C == cov(X, wv1, 1)

        (m, C) = mean_and_cov(X, wv1, 1)
        @test m == mean(X, wv1, 1)
        @test C == cov(X, wv1, 1)

        (m, C) = mean_and_cov(X, wv2, 2)
        @test m == mean(X, wv2, 2)
        @test C == cov(X, wv2, 2)
    end
end

end # @testset "StatsBase.Deprecates"
