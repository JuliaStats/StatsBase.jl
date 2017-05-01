using StatsBase
using Base.Test

@testset "StatsBase.Moments" begin
weight_funcs = (aweights, fweights, pweights)

@testset "Variance and Standard Deviation" begin
    @testset "Vectors" begin
        x = [0.57, 0.10, 0.91, 0.72, 0.46]
        w = [3.84, 2.70, 8.29, 8.91, 9.71]

        @testset "Uncorrected" begin
            @testset "Variance with $f" for f in weight_funcs
                wv = f(w)
                m = mean(x, wv)
                @test var(x, wv, false)           ≈ sum(abs2.(x .- m), wv) ./ sum(wv)
                @test var(x, wv, false; mean=0)   ≈ sum(abs2.(x), wv) ./ sum(wv)
                @test var(x, wv, false; mean=1.0) ≈ sum(abs2.(x .- 1.0), wv) ./ sum(wv)
            end

            @testset "Standard Deviation with $f" for f in weight_funcs
                wv = f(w)
                m = mean(x, wv)
                @test std(x, wv, false)           ≈ sqrt(var(x, wv, false))
                @test std(x, wv, false; mean=0)   ≈ sqrt(var(x, wv, false; mean=0))
                @test std(x, wv, false; mean=1.0) ≈ sqrt(var(x, wv, false; mean=1.0))
            end

            @testset "Mean and Variance with $f" for f in weight_funcs
                wv = f(w)
                (m, v) = mean_and_var(x, false)
                @test m == mean(x)
                @test v == var(x, corrected=false)

                (m, v) = mean_and_var(x, wv, false)
                @test m == mean(x, wv)
                @test v == var(x, wv, false)
            end

            @testset "Mean and Standard Deviation with $f" for f in weight_funcs
                wv = f(w)
                (m, s) = mean_and_std(x, false)
                @test m == mean(x)
                @test s == std(x, corrected=false)

                (m, s) = mean_and_std(x, wv, false)
                @test m == mean(x, wv)
                @test s == std(x, wv, false)
            end
        end

        @testset "Corrected" begin
            @testset "Variance" begin
                # expected `var` output for (aweights, fweights, pweights)
                expected = (0.0694434191182236, 0.05466601256158146, 0.06628969012045285)
                expected_0 = (0.5798908707332937, 0.45649137134052387, 0.5535554932735426)
                expected_1 = (0.25422659392845115, 0.20012773497688754, 0.24268105381165922)

                @testset "$(weight_funcs[i])" for i in 1:3
                    wv = weight_funcs[i](w)
                    m = mean(x, wv)

                    @test var(x, wv, true)           ≈ expected[i]
                    @test var(x, wv, true; mean=0)   ≈ expected_0[i]
                    @test var(x, wv, true; mean=1.0) ≈ expected_1[i]
                end
            end

            @testset "Standard Deviation with $f" for f in weight_funcs
                wv = f(w)
                m = mean(x, wv)
                @test std(x, wv, true)           ≈ sqrt(var(x, wv, true))
                @test std(x, wv, true; mean=0)   ≈ sqrt(var(x, wv, true; mean=0))
                @test std(x, wv, true; mean=1.0) ≈ sqrt(var(x, wv, true; mean=1.0))
            end

            @testset "Mean and Variance with $f" for f in weight_funcs
                wv = f(w)

                (m, v) = mean_and_var(x, true)
                @test m == mean(x)
                @test v == var(x, corrected=true)

                (m, v) = mean_and_var(x, wv, true)
                @test m == mean(x, wv)
                @test v == var(x, wv, true)
            end

            @testset "Mean and Standard Deviation with $f" for f in weight_funcs
                wv = f(w)

                (m, s) = mean_and_std(x, true)
                @test m == mean(x)
                @test s == std(x, corrected=true)

                (m, s) = mean_and_std(x, wv, true)
                @test m == mean(x, wv)
                @test s == std(x, wv, true)
            end
        end
    end

    @testset "Matrices" begin
        x = rand(5, 6)
        w1 = rand(5)
        w2 = rand(6)
        wv1 = fweights(w1)
        wv2 = fweights(w2)
        m1 = mean(x, wv1, 1)
        m2 = mean(x, wv2, 2)

        @test var(x, wv1, 1, false; mean=0) ≈ sum(abs2.(x) .* w1, 1) ./ sum(wv1)
        @test var(x, wv2, 2, false; mean=0) ≈ sum(abs2.(x) .* w2', 2) ./ sum(wv2)

        @test var(x, wv1, 1, false; mean=m1) ≈ sum(abs2.(x .- m1) .* w1, 1) ./ sum(wv1)
        @test var(x, wv2, 2, false; mean=m2) ≈ sum(abs2.(x .- m2) .* w2', 2) ./ sum(wv2)

        @test var(x, wv1, 1, false) ≈ sum(abs2.(x .- m1) .* w1, 1) ./ sum(wv1)
        @test var(x, wv2, 2, false) ≈ sum(abs2.(x .- m2) .* w2', 2) ./ sum(wv2)

        @test std(x, wv1, 1, false)          ≈ sqrt.(var(x, wv1, 1, false))
        @test std(x, wv2, 2, false)          ≈ sqrt.(var(x, wv2, 2, false))
        @test std(x, wv1, 1, false; mean=0)  ≈ sqrt.(var(x, wv1, 1, false; mean=0))
        @test std(x, wv2, 2, false; mean=0)  ≈ sqrt.(var(x, wv2, 2, false; mean=0))
        @test std(x, wv1, 1, false; mean=m1) ≈ sqrt.(var(x, wv1, 1, false; mean=m1))
        @test std(x, wv2, 2, false; mean=m2) ≈ sqrt.(var(x, wv2, 2, false; mean=m2))

        for d in 1:2
            (m, v) = mean_and_var(x, d, false)
            @test m == mean(x, d)
            @test v == var(x, d; corrected=false)

            (m, s) = mean_and_std(x, d, false)
            @test m == mean(x, d)
            @test s == std(x, d; corrected=false)
        end

        (m, v) = mean_and_var(x, wv1, 1, true)
        @test m == mean(x, wv1, 1)
        @test v == var(x, wv1, 1, true)

        (m, v) = mean_and_var(x, wv2, 2, false)
        @test m == mean(x, wv2, 2)
        @test v == var(x, wv2, 2, false)

        (m, s) = mean_and_std(x, wv1, 1, false)
        @test m == mean(x, wv1, 1)
        @test s == std(x, wv1, 1, false)

        (m, s) = mean_and_std(x, wv2, 2, false)
        @test m == mean(x, wv2, 2)
        @test s == std(x, wv2, 2, false)
    end
end

@testset "Skewness and Kurtosis" begin
    wv = fweights(ones(5) * 2.0)

    @test skewness(1:5)             ≈  0.0
    @test skewness([1, 2, 3, 4, 5]) ≈  0.0
    @test skewness([1, 2, 2, 2, 5]) ≈  1.1731251294063556
    @test skewness([1, 4, 4, 4, 5]) ≈ -1.1731251294063556

    @test skewness([1, 2, 2, 2, 5], wv) ≈ 1.1731251294063556

    @test kurtosis(1:5)             ≈ -1.3
    @test kurtosis([1, 2, 3, 4, 5]) ≈ -1.3
    @test kurtosis([1, 2, 3, 3, 2]) ≈ -1.1530612244897953

    @test kurtosis([1, 2, 3, 4, 5], wv) ≈ -1.3
end

@testset "General Moments" begin
    x = collect(2.0:8.0)
    @test moment(x, 2) ≈ sum((x .- 5).^2) / length(x)
    @test moment(x, 3) ≈ sum((x .- 5).^3) / length(x)
    @test moment(x, 4) ≈ sum((x .- 5).^4) / length(x)
    @test moment(x, 5) ≈ sum((x .- 5).^5) / length(x)

    @test moment(x, 2, 4.0) ≈ sum((x .- 4).^2) / length(x)
    @test moment(x, 3, 4.0) ≈ sum((x .- 4).^3) / length(x)
    @test moment(x, 4, 4.0) ≈ sum((x .- 4).^4) / length(x)
    @test moment(x, 5, 4.0) ≈ sum((x .- 4).^5) / length(x)

    w = fweights([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
    x2 = collect(2.0:6.0)
    @test moment(x, 2, w) ≈ sum((x2 .- 4).^2) / 5
    @test moment(x, 3, w) ≈ sum((x2 .- 4).^3) / 5
    @test moment(x, 4, w) ≈ sum((x2 .- 4).^4) / 5
    @test moment(x, 5, w) ≈ sum((x2 .- 4).^5) / 5
end

end # @testset "StatsBase.Moments"
