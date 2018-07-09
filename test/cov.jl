using StatsBase
using LinearAlgebra, Random, Test

@testset "StatsBase.Covariance" begin
weight_funcs = (weights, aweights, fweights, pweights)

@testset "$f" for f in weight_funcs
    X = randn(3, 8)

    Z1 = X .- mean(X, dims = 1)
    Z2 = X .- mean(X, dims = 2)

    w1 = rand(3)
    w2 = rand(8)

    # varcorrection is negative if sum of weights is smaller than 1
    if f === fweights
        w1[1] += 1
        w2[1] += 1
    end

    wv1 = f(w1)
    wv2 = f(w2)

    Z1w = X .- mean(X, wv1, 1)
    Z2w = X .- mean(X, wv2, 2)

    ## reference results

    S1 = Z1'Z1
    S2 = Z2 * Z2'

    Sz1 = X'X
    Sz2 = X * X'

    S1w = Z1w' * Matrix(Diagonal(w1)) * Z1w
    S2w = Z2w * Matrix(Diagonal(w2)) * Z2w'

    Sz1w = X' * Matrix(Diagonal(w1)) * X
    Sz2w = X * Matrix(Diagonal(w2)) * X'

    @testset "Scattermat" begin
        @test scattermat(X)    ≈ S1
        @test scattermat(X, 2) ≈ S2

        @test StatsBase.scattermatm(X, 0)    ≈ Sz1
        @test StatsBase.scattermatm(X, 0, 2) ≈ Sz2

        @test StatsBase.scattermatm(X, mean(X, dims = 1))    ≈ S1
        @test StatsBase.scattermatm(X, mean(X, dims = 2), 2) ≈ S2

        @test StatsBase.scattermatm(X, zeros(1,8))  ≈ Sz1
        @test StatsBase.scattermatm(X, zeros(3), 2) ≈ Sz2

        @testset "Weighted" begin
            @test scattermat(X, wv1)    ≈ S1w
            @test scattermat(X, wv2, 2) ≈ S2w

            @test StatsBase.scattermatm(X, 0, wv1)    ≈ Sz1w
            @test StatsBase.scattermatm(X, 0, wv2, 2) ≈ Sz2w

            @test StatsBase.scattermatm(X, mean(X, wv1, 1), wv1)    ≈ S1w
            @test StatsBase.scattermatm(X, mean(X, wv2, 2), wv2, 2) ≈ S2w

            @test StatsBase.scattermatm(X, zeros(1,8), wv1)  ≈ Sz1w
            @test StatsBase.scattermatm(X, zeros(3), wv2, 2) ≈ Sz2w
        end
    end

    @testset "Uncorrected" begin
        @testset "Weighted Covariance" begin
            @test cov(X, wv1; corrected=false)    ≈ S1w ./ sum(wv1)
            @test cov(X, wv2, 2; corrected=false) ≈ S2w ./ sum(wv2)

            @test StatsBase.covm(X, 0, wv1, 1; corrected=false) ≈ Sz1w ./ sum(wv1)
            @test StatsBase.covm(X, 0, wv2, 2; corrected=false) ≈ Sz2w ./ sum(wv2)

            @test StatsBase.covm(X, mean(X, wv1, 1), wv1, 1; corrected=false) ≈ S1w ./ sum(wv1)
            @test StatsBase.covm(X, mean(X, wv2, 2), wv2, 2; corrected=false) ≈ S2w ./ sum(wv2)

            @test StatsBase.covm(X, zeros(1,8), wv1, 1; corrected=false) ≈ Sz1w ./ sum(wv1)
            @test StatsBase.covm(X, zeros(3), wv2, 2; corrected=false)   ≈ Sz2w ./ sum(wv2)
        end

        @testset "Mean and covariance" begin
            (m, C) = mean_and_cov(X; corrected=false)
            @test m == mean(X, dims = 1)
            @test C == cov(X, dims = 1, corrected=false)

            (m, C) = mean_and_cov(X, 1; corrected=false)
            @test m == mean(X, dims = 1)
            @test C == cov(X, dims = 1, corrected = false)

            (m, C) = mean_and_cov(X, 2; corrected=false)
            @test m == mean(X, dims = 2)
            @test C == cov(X, dims = 2, corrected = false)

            (m, C) = mean_and_cov(X, wv1; corrected=false)
            @test m == mean(X, wv1, 1)
            @test C == cov(X, wv1, 1, corrected=false)

            (m, C) = mean_and_cov(X, wv1, 1; corrected=false)
            @test m == mean(X, wv1, 1)
            @test C == cov(X, wv1, 1, corrected=false)

            (m, C) = mean_and_cov(X, wv2, 2; corrected=false)
            @test m == mean(X, wv2, 2)
            @test C == cov(X, wv2, 2, corrected=false)
        end
        @testset "Conversions" begin
            std1 = std(X, wv1, 1; corrected=false)
            std2 = std(X, wv2, 2; corrected=false)

            cov1 = cov(X, wv1, 1; corrected=false)
            cov2 = cov(X, wv2, 2; corrected=false)

            cor1 = cor(X, wv1, 1)
            cor2 = cor(X, wv2, 2)

            @testset "cov2cor" begin
                @test cov2cor(cov(X, dims = 1), std(X, dims = 1)) ≈ cor(X, dims = 1)
                @test cov2cor(cov(X, dims = 2), std(X, dims = 2)) ≈ cor(X, dims = 2)
                @test cov2cor(cov1, std1) ≈ cor1
                @test cov2cor(cov2, std2) ≈ cor2
            end
            @testset "cor2cov" begin
                @test cor2cov(cor(X, dims = 1), std(X, dims = 1)) ≈ cov(X, dims = 1)
                @test cor2cov(cor(X, dims = 2), std(X, dims = 2)) ≈ cov(X, dims = 2)
                @test cor2cov(cor1, std1) ≈ cov1
                @test cor2cov(cor2, std2) ≈ cov2
            end
        end
    end

    @testset "Corrected" begin
        @testset "Weighted Covariance" begin
            if isa(wv1, Weights)
                @test_throws ArgumentError cov(X, wv1; corrected=true)
            else
                var_corr1 = StatsBase.varcorrection(wv1, true)
                var_corr2 = StatsBase.varcorrection(wv2, true)

                @test cov(X, wv1; corrected=true)    ≈ S1w .* var_corr1
                @test cov(X, wv2, 2; corrected=true) ≈ S2w .* var_corr2

                @test StatsBase.covm(X, 0, wv1, 1; corrected=true) ≈ Sz1w .* var_corr1
                @test StatsBase.covm(X, 0, wv2, 2; corrected=true) ≈ Sz2w .* var_corr2

                @test StatsBase.covm(X, mean(X, wv1, 1), wv1, 1; corrected=true) ≈ S1w .* var_corr1
                @test StatsBase.covm(X, mean(X, wv2, 2), wv2, 2; corrected=true) ≈ S2w .* var_corr2

                @test StatsBase.covm(X, zeros(1,8), wv1, 1; corrected=true) ≈ Sz1w .* var_corr1
                @test StatsBase.covm(X, zeros(3), wv2, 2; corrected=true)   ≈ Sz2w .* var_corr2
            end
        end
        @testset "Mean and covariance" begin
            (m, C) = mean_and_cov(X; corrected=true)
            @test m == mean(X, dims =1)
            @test C == cov(X, dims = 1, corrected = true)

            (m, C) = mean_and_cov(X, 1; corrected=true)
            @test m == mean(X, dims = 1)
            @test C == cov(X, dims = 1, corrected = true)

            (m, C) = mean_and_cov(X, 2; corrected=true)
            @test m == mean(X, dims = 2)
            @test C == cov(X, dims = 2, corrected = true)

            if isa(wv1, Weights)
                @test_throws ArgumentError mean_and_cov(X, wv1; corrected=true)
            else
                (m, C) = mean_and_cov(X, wv1; corrected=true)
                @test m == mean(X, wv1, 1)
                @test C == cov(X, wv1, 1; corrected=true)

                (m, C) = mean_and_cov(X, wv1, 1; corrected=true)
                @test m == mean(X, wv1, 1)
                @test C == cov(X, wv1, 1; corrected=true)

                (m, C) = mean_and_cov(X, wv2, 2; corrected=true)
                @test m == mean(X, wv2, 2)
                @test C == cov(X, wv2, 2; corrected=true)
            end
        end
        @testset "Conversions" begin
            if !isa(wv1, Weights)
                std1 = std(X, wv1, 1; corrected=true)
                std2 = std(X, wv2, 2; corrected=true)

                cov1 = cov(X, wv1, 1; corrected=true)
                cov2 = cov(X, wv2, 2; corrected=true)

                cor1 = cor(X, wv1, 1)
                cor2 = cor(X, wv2, 2)

                @testset "cov2cor" begin
                    @test cov2cor(cov(X, dims = 1), std(X, dims = 1)) ≈ cor(X, dims = 1)
                    @test cov2cor(cov(X, dims = 2), std(X, dims = 2)) ≈ cor(X, dims = 2)
                    @test cov2cor(cov1, std1) ≈ cor1
                    @test cov2cor(cov2, std2) ≈ cor2
                end

                @testset "cov2cor!" begin
                    tmp_cov1 = copy(cov1)
                    @test !(tmp_cov1 ≈ cor1)
                    StatsBase.cov2cor!(tmp_cov1, std1)
                    @test tmp_cov1 ≈ cor1

                    tmp_cov2 = copy(cov2)
                    @test !(tmp_cov2 ≈ cor2)
                    StatsBase.cov2cor!(tmp_cov2, std2)
                    @test tmp_cov2 ≈ cor2
                end

                @testset "cor2cov" begin
                    @test cor2cov(cor(X, dims = 1), std(X, dims = 1)) ≈ cov(X, dims = 1)
                    @test cor2cov(cor(X, dims = 2), std(X, dims = 2)) ≈ cov(X, dims = 2)
                    @test cor2cov(cor1, std1) ≈ cov1
                    @test cor2cov(cor2, std2) ≈ cov2
                end

                @testset "cor2cov!" begin
                    tmp_cor1 = copy(cor1)
                    @test !(tmp_cor1 ≈ cov1)
                    StatsBase.cor2cov!(tmp_cor1, std1)
                    @test tmp_cor1 ≈ cov1

                    tmp_cor2 = copy(cor2)
                    @test !(tmp_cor2 ≈ cov2)
                    StatsBase.cor2cov!(tmp_cor2, std2)
                    @test tmp_cor2 ≈ cov2
                end
            end
        end
    end

    @testset "Correlation" begin
        @test cor(X, f(ones(3)), 1) ≈ cor(X, dims = 1)
        @test cor(X, f(ones(8)), 2) ≈ cor(X, dims = 2)

        cov1 = cov(X, wv1, 1; corrected=false)
        std1 = std(X, wv1, 1; corrected=false)
        cov2 = cov(X, wv2, 2; corrected=false)
        std2 = std(X, wv2, 2; corrected=false)
        expected_cor1 = StatsBase.cov2cor!(cov1, std1)
        expected_cor2 = StatsBase.cov2cor!(cov2, std2)

        @test cor(X, wv1, 1) ≈ expected_cor1
        @test cor(X, wv2, 2) ≈ expected_cor2
    end
end
end # @testset "StatsBase.Covariance"
