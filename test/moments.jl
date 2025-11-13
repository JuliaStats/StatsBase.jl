using StatsBase
using Statistics
using Test

@testset "StatsBase.Moments" begin
weight_funcs = (weights, aweights, fweights, pweights)

##### weighted var & std

x = [0.57, 0.10, 0.91, 0.72, 0.46, 0.0]
xf0 = Float32.(x)
w = [3.84, 2.70, 8.29, 8.91, 9.71, 0.0]
wf0 = Float32.(w)

@testset "Uncorrected with $f (values of type $(eltype(x)), weights of type $(eltype(w)))" for f in weight_funcs, x in (x, xf0), w in (w, wf0)
    TX = eltype(x)
    T = promote_type(TX, eltype(w))
    wv = f(w)
    m = @inferred(mean(x, wv))::T

    # expected uncorrected output
    expected_var = sum(abs2.(x .- m), wv) / sum(wv)
    expected_std = sqrt.(expected_var)

    @testset "Variance" begin
        @test @inferred(var(x, wv; corrected=false))::T           ≈ expected_var
        @test @inferred(var(x, wv; mean=m, corrected=false))::T   ≈ expected_var
        @test @inferred(varm(x, wv, m; corrected=false))::T       ≈ expected_var
    end

    @testset "Standard Deviation" begin
        @test @inferred(std(x, wv; corrected=false))::T           ≈ expected_std
        @test @inferred(std(x, wv; mean=m, corrected=false))::T   ≈ expected_std
        @test @inferred(stdm(x, wv, m; corrected=false))::T       ≈ expected_std
    end

    @testset "Mean and Variance" begin
        (m, v) = @inferred(mean_and_var(x; corrected=false))::Tuple{TX,TX}
        @test m == mean(x)
        @test v == var(x; corrected=corrected=false)

        (m, v) = @inferred(mean_and_var(x, wv; corrected=false))::Tuple{T,T}
        @test m == mean(x, wv)
        @test v == var(x, wv; corrected=false)
    end

    @testset "Mean and Standard Deviation" begin
        (m, s) = @inferred(mean_and_std(x; corrected=false))::Tuple{TX,TX}
        @test m == mean(x)
        @test s == std(x; corrected=false)

        (m, s) = @inferred(mean_and_std(x, wv; corrected=false))::Tuple{T,T}
        @test m == mean(x, wv)
        @test s == std(x, wv; corrected=false)
    end
end

# expected corrected output for (weights, aweights, fweights, pweights)
expected_var = [NaN, 0.0694434191182236, 0.05466601256158146, 0.06628969012045285]
expected_std = sqrt.(expected_var)

@testset "Corrected with $(weight_funcs[i]) (values of type $(eltype(x)), weights of type $(eltype(w)))" for i in eachindex(weight_funcs), x in (x, xf0), w in (w, wf0)
    TX = eltype(x)
    TW = eltype(w)
    T = promote_type(TX, TW)
    TR = TX === Float32 || TW === Float32 ? Float32 : Float64
    wv = weight_funcs[i](w)
    m = @inferred(mean(x, wv))::T

    @testset "Variance" begin
        if isa(wv, Weights)
            @test_throws ArgumentError var(x, wv; corrected=true)
        else
            @test @inferred(var(x, wv; corrected=true))::T           ≈ TR(expected_var[i])
            @test @inferred(var(x, wv; mean=m, corrected=true))::T   ≈ TR(expected_var[i])
            @test @inferred(varm(x, wv, m; corrected=true))::T       ≈ TR(expected_var[i])
        end
    end

    @testset "Standard Deviation" begin
        if isa(wv, Weights)
            @test_throws ArgumentError std(x, wv; corrected=true)
        else
            @test @inferred(std(x, wv; corrected=true))::T           ≈ TR(expected_std[i])
            @test @inferred(std(x, wv; mean=m, corrected=true))::T   ≈ TR(expected_std[i])
            @test @inferred(stdm(x, wv, m; corrected=true))::T       ≈ TR(expected_std[i])
        end
    end

    @testset "Mean and Variance" begin
        (m, v) = @inferred(mean_and_var(x; corrected=true))::Tuple{TX,TX}
        @test m == mean(x)
        @test v == var(x; corrected=true)

        if isa(wv, Weights)
            @test_throws ArgumentError mean_and_var(x, wv; corrected=true)
        else
            (m, v) = @inferred(mean_and_var(x, wv; corrected=true))::Tuple{T,T}
            @test m == mean(x, wv)
            @test v == var(x, wv; corrected=true)
        end
    end

    @testset "Mean and Standard Deviation" begin
        (m, s) = @inferred(mean_and_std(x; corrected=true))::Tuple{TX,TX}
        @test m == mean(x)
        @test s == std(x; corrected=true)

        if isa(wv, Weights)
            @test_throws ArgumentError mean_and_std(x, wv; corrected=true)
        else
            (m, s) = @inferred(mean_and_std(x, wv; corrected=true))::Tuple{T,T}
            @test m == mean(x, wv)
            @test s == std(x, wv; corrected=true)
        end
    end
end

x = rand(5, 6)
xf0 = Float32.(x)
w1 = [0.57, 5.10, 0.91, 1.72, 0.0]
w1f0 = Float32.(w1)
w2 = [3.84, 2.70, 8.29, 8.91, 9.71, 0.0]
w2f0 = Float32.(w2)

@testset "Uncorrected with $f (values of type $(eltype(x)), 1st weights of type $(eltype(w1)), 2nd weights of type $(eltype(w2)))" for f in weight_funcs, x in (x, xf0), w1 in (w1, w1f0), w2 in (w2, w2f0)
    TX = eltype(x)
    TW1 = eltype(w1)
    TW2 = eltype(w2)
    T1 = promote_type(TX, TW1)
    T2 = promote_type(TX, TW2)
    wv1 = f(w1)
    wv2 = f(w2)
    m1 = @inferred(mean(x, wv1, dims=1))::Matrix{T1}
    m2 = @inferred(mean(x, wv2, dims=2))::Matrix{T2}

    expected_var1 = sum(abs2.(x .- m1) .* w1, dims = 1) ./ sum(wv1)
    expected_var2 = sum(abs2.(x .- m2) .* w2', dims = 2) ./ sum(wv2)
    expected_std1 = sqrt.(expected_var1)
    expected_std2 = sqrt.(expected_var2)

    @testset "Variance" begin
        @test @inferred(var(x, wv1, 1; corrected=false))::Matrix{T1}          ≈ expected_var1
        @test @inferred(var(x, wv2, 2; corrected=false))::Matrix{T2}          ≈ expected_var2
        @test @inferred(var(x, wv1, 1; mean=m1, corrected=false))::Matrix{T1} ≈ expected_var1
        @test @inferred(var(x, wv2, 2; mean=m2, corrected=false))::Matrix{T2} ≈ expected_var2
        @test @inferred(varm(x, wv1, m1, 1; corrected=false))::Matrix{T1}     ≈ expected_var1
        @test @inferred(varm(x, wv2, m2, 2; corrected=false))::Matrix{T2}     ≈ expected_var2
    end

    @testset "Standard Deviation" begin
        @test @inferred(std(x, wv1, 1; corrected=false))::Matrix{T1}          ≈ expected_std1
        @test @inferred(std(x, wv2, 2; corrected=false))::Matrix{T2}         ≈ expected_std2
        @test @inferred(std(x, wv1, 1; mean=m1, corrected=false))::Matrix{T1} ≈ expected_std1
        @test @inferred(std(x, wv2, 2; mean=m2, corrected=false))::Matrix{T2} ≈ expected_std2
        @test @inferred(stdm(x, wv1, m1, 1; corrected=false))::Matrix{T1}     ≈ expected_std1
        @test @inferred(stdm(x, wv2, m2, 2; corrected=false))::Matrix{T2}     ≈ expected_std2
    end

    @testset "Mean and Variance" begin
        for d in 1:2
            (m, v) = @inferred(mean_and_var(x, d; corrected=false))::Tuple{Matrix{TX},Matrix{TX}}
            @test m == mean(x, dims=d)
            @test v == var(x, dims=d, corrected=false)
        end

        (m, v) = @inferred(mean_and_var(x, wv1, 1; corrected=false))::Tuple{Matrix{T1},Matrix{T1}}
        @test m == mean(x, wv1, dims=1)
        @test v == var(x, wv1, 1; corrected=false)

        (m, v) = @inferred(mean_and_var(x, wv2, 2; corrected=false))::Tuple{Matrix{T2},Matrix{T2}}
        @test m == mean(x, wv2, dims=2)
        @test v == var(x, wv2, 2; corrected=false)
    end

    @testset "Mean and Standard Deviation" begin
        for d in 1:2
            (m, s) = @inferred(mean_and_std(x, d; corrected=false))::Tuple{Matrix{TX},Matrix{TX}}
            @test m == mean(x, dims=d)
            @test s == std(x, dims=d; corrected=false)
        end

        (m, s) = @inferred(mean_and_std(x, wv1, 1; corrected=false))::Tuple{Matrix{T1},Matrix{T1}}
        @test m == mean(x, wv1, dims=1)
        @test s == std(x, wv1, 1; corrected=false)

        (m, s) = @inferred(mean_and_std(x, wv2, 2; corrected=false))::Tuple{Matrix{T2},Matrix{T2}}
        @test m == mean(x, wv2, dims=2)
        @test s == std(x, wv2, 2; corrected=false)
    end
end

@testset "Corrected with $f (values of type $(eltype(x)), weights of type $(eltype(w1)))" for f in weight_funcs, x in (Float32.(x), Float64.(x)), (w1, w2) in ((Float32.(w1), Float32.(w2)), (Float64.(w1), Float64.(w2)))
    TX = eltype(x)
    TW1 = eltype(w1)
    TW2 = eltype(w2)
    T1 = promote_type(TX, TW1)
    T2 = promote_type(TX, TW2)
    wv1 = f(w1)
    wv2 = f(w2)
    m1 = @inferred(mean(x, wv1, dims=1))::Matrix{T1}
    m2 = @inferred(mean(x, wv2, dims=2))::Matrix{T2}

    if !isa(wv1, Weights)
        expected_var1 = sum(abs2.(x .- m1) .* w1, dims = 1) .* StatsBase.varcorrection(wv1, true)
        expected_var2 = sum(abs2.(x .- m2) .* w2', dims = 2) .* StatsBase.varcorrection(wv2, true)
        expected_std1 = sqrt.(expected_var1)
        expected_std2 = sqrt.(expected_var2)
    end

    @testset "Variance" begin
        if isa(wv1, Weights)
            @test_throws ArgumentError var(x, wv1, 1; corrected=true)
        else
            @test @inferred(var(x, wv1, 1; corrected=true))::Matrix{T1}          ≈ expected_var1
            @test @inferred(var(x, wv2, 2; corrected=true))::Matrix{T2}          ≈ expected_var2
            @test @inferred(var(x, wv1, 1; mean=m1, corrected=true))::Matrix{T1} ≈ expected_var1
            @test @inferred(var(x, wv2, 2; mean=m2, corrected=true))::Matrix{T2} ≈ expected_var2
            @test @inferred(varm(x, wv1, m1, 1; corrected=true))::Matrix{T1}     ≈ expected_var1
            @test @inferred(varm(x, wv2, m2, 2; corrected=true))::Matrix{T2}     ≈ expected_var2
        end
    end

    @testset "Standard Deviation" begin
        if isa(wv1, Weights)
            @test_throws ArgumentError std(x, wv1, 1; corrected=true)
        else
            @test @inferred(std(x, wv1, 1; corrected=true))::Matrix{T1}          ≈ expected_std1
            @test @inferred(std(x, wv2, 2; corrected=true))::Matrix{T2}          ≈ expected_std2
            @test @inferred(std(x, wv1, 1; mean=m1, corrected=true))::Matrix{T1} ≈ expected_std1
            @test @inferred(std(x, wv2, 2; mean=m2, corrected=true))::Matrix{T2} ≈ expected_std2
            @test @inferred(stdm(x, wv1, m1, 1; corrected=true))::Matrix{T1}     ≈ expected_std1
            @test @inferred(stdm(x, wv2, m2, 2; corrected=true))::Matrix{T2}     ≈ expected_std2
        end
    end

    @testset "Mean and Variance" begin
        for d in 1:2
            (m, v) = @inferred(mean_and_var(x, d; corrected=true))::Tuple{Matrix{TX},Matrix{TX}}
            @test m == mean(x, dims=d)
            @test v == var(x, dims=d, corrected=true)
        end

        if isa(wv1, Weights)
            @test_throws ArgumentError mean_and_var(x, wv1, 1; corrected=true)
        else
            (m, v) = @inferred(mean_and_var(x, wv1, 1; corrected=true))::Tuple{Matrix{T1},Matrix{T1}}
            @test m == mean(x, wv1, dims=1)
            @test v == var(x, wv1, 1; corrected=true)

            (m, v) = @inferred(mean_and_var(x, wv2, 2; corrected=true))::Tuple{Matrix{T2},Matrix{T2}}
            @test m == mean(x, wv2, dims=2)
            @test v == var(x, wv2, 2; corrected=true)
        end
    end

    @testset "Mean and Standard Deviation" begin
        for d in 1:2
            (m, s) = @inferred(mean_and_std(x, d; corrected=true))::Tuple{Matrix{TX},Matrix{TX}}
            @test m == mean(x, dims=d)
            @test s == std(x, dims=d, corrected=true)
        end

        if isa(wv1, Weights)
            @test_throws ArgumentError mean_and_std(x, wv1, 1; corrected=true)
        else
            (m, s) = @inferred(mean_and_std(x, wv1, 1; corrected=true))::Tuple{Matrix{T1},Matrix{T1}}
            @test m == mean(x, wv1, dims=1)
            @test s == std(x, wv1, 1; corrected=true)

            (m, s) = @inferred(mean_and_std(x, wv2, 2; corrected=true))::Tuple{Matrix{T2},Matrix{T2}}
            @test m == mean(x, wv2, dims=2)
            @test s == std(x, wv2, 2; corrected=true)
        end
    end
end

@testset "Skewness and Kurtosis with $f" for f in weight_funcs
    for T in (Int, Float32, Float64)
        for v in (T(1):T(5), collect(T, 1:5))
            s = @inferred(skewness(v))
            @test s isa float(T)
            @test iszero(s)

            k = @inferred(kurtosis(v))
            @test k isa float(T)
            @test k ≈ oftype(k, -1.3)
        end

        v = T[1, 2, 2, 2, 5]
        s = @inferred(skewness(v))
        @test s isa float(T)
        @test s ≈ oftype(s, 1.1731251294063556)

        v = T[1, 4, 4, 4, 5]
        s = @inferred(skewness(v))
        @test s isa float(T)
        @test s ≈ oftype(s, -1.1731251294063556)

        v = T[1, 2, 3, 3, 2]
        k = @inferred(kurtosis(v))
        @test k isa float(T)
        @test k ≈ oftype(k, -1.1530612244897953)

        # Empty arrays
        s = @inferred(skewness(T[]))
        @test s isa float(T)
        @test isnan(s)
        k = @inferred(kurtosis(T[]))
        @test k isa float(T)
        @test isnan(k)

        for T2 in (Int, Float32, Float64)
            wv = f(fill(T2(2), 5))
            v = T[1, 2, 2, 2, 5]
            s = @inferred(skewness(v, wv))
            @test s isa float(promote_type(T, T2))
            @test s ≈ oftype(s, 1.1731251294063556)

            v = collect(T, 1:5)
            k = @inferred(kurtosis(v, wv))
            @test k isa float(promote_type(T, T2))
            @test k ≈ oftype(k, -1.3)

            # Empty arrays
            wv = f(T2[])
            s = @inferred(skewness(T[], wv))
            @test s isa float(promote_type(T, T2))
            @test isnan(s)
            k = @inferred(kurtosis(T[], wv))
            @test k isa float(promote_type(T, T2))
            @test isnan(k)
        end

        # Invalid arguments
        v = collect(T, 1:5)
        for n in (length(x) - 1, length(x) + 1)
            @test_throws DimensionMismatch("Inconsistent array lengths.") kurtosis(v, f(ones(T, n)))
            @test_throws DimensionMismatch("Inconsistent array lengths.") skewness(v, f(ones(T, n)))
        end
    end
end

@testset "General Moments with $f" for f in weight_funcs
    for T in (Int, Float32, Float64)
        x = collect(T, 2:8)
        for k in 2:5
            momk = @inferred(moment(x, k))
            @test momk isa float(T)
            @test momk ≈ sum((x .- 5).^k) / length(x)

            # Empty array
            momk = @inferred(moment(T[], k))
            @test momk isa float(T)
            @test isnan(momk)

            for TM in (Int, Float32, Float64)
                m = TM(4)
                momk = @inferred(moment(x, k, m))
                @test momk isa float(promote_type(T, TM))
                @test momk ≈ sum((x .- 4).^k) / length(x)

                # Empty array
                momk = @inferred(moment(T[], k, zero(TM)))
                @test momk isa float(promote_type(T, TM))
                @test isnan(momk)
            end
        end

        for T2 in (Int, Float32, Float64)
            wv = f(T2[1, 1, 1, 1, 1, 0, 0])
            x2 = collect(T, 2:6)
            for k in 2:5
                momk = @inferred(moment(x, k, wv))
                @test momk isa float(promote_type(T, T2))
                @test momk ≈ sum((x2 .- 4).^k) / 5

                # Empty array
                momk = @inferred(moment(T[], k, f(T2[])))
                @test momk isa float(promote_type(T, T2))
                @test isnan(momk)

                for TM in (Int, Float32, Float64)
                    m = TM(3)
                    momk = @inferred(moment(x, k, wv, m))
                    @test momk isa float(promote_type(T, T2, TM))
                    @test momk ≈ sum((x2 .- 3).^k) / 5

                    # Empty array
                    momk = @inferred(moment(T[], k, f(T2[]), zero(TM)))
                    @test momk isa float(promote_type(T, T2, TM))
                    @test isnan(momk)
                end
            end
        end
    end
end

@testset "Cumulants with $f" for f in weight_funcs
    for T in (Int, Float32, Float64)
        x = collect(T, 2:8)
        for k in 1:6
            cumk = @inferred(cumulant(x, k))
            @test cumk isa float(T)
            if k == 1
                @test cumk ≈ mean(x)
            elseif k == 2 || k == 3
                @test cumk ≈ moment(x, k)
            elseif k == 4
                @test cumk ≈ moment(x, 4) - 3*moment(x, 2)^2
            elseif k == 5
                @test cumk ≈ moment(x, 5) - 10*moment(x, 3)*moment(x, 2)
            else
                @assert k == 6
                @test cumk ≈ moment(x, 6) - 15*moment(x, 4)*moment(x, 2) - 10*moment(x, 3)^2 + 30*moment(x, 2)^3
            end
        end
        cumks = @inferred(cumulant(x, 1:6))
        @test cumks isa Vector{float(T)}
        @test cumks == [cumulant(x, i) for i in 1:6]

        for TM in (Int, Float32, Float64)
            m = TM(4)
            for k in 1:6
                cumk = @inferred(cumulant(x, k, m))
                @test cumk isa float(promote_type(T, TM))
                if k == 1
                    @test cumk ≈ m
                elseif k == 2 || k == 3
                    @test cumk ≈ moment(x, k, m)
                elseif k == 4
                    @test cumk ≈ moment(x, 4, m) - 3*moment(x, 2, m)^2
                elseif k == 5
                    @test cumk ≈ moment(x, 5, m) - 10*moment(x, 3, m)*moment(x, 2, m)
                else
                    @assert k == 6
                    @test cumk ≈ moment(x, 6, m) - 15*moment(x, 4, m)*moment(x, 2, m) - 10*moment(x, 3, m)^2 + 30*moment(x, 2, m)^3
                end
            end
            cumks = @inferred(cumulant(x, 1:6, m))
            @test cumks isa Vector{float(promote_type(T, TM))}
            @test cumks == [cumulant(x, i, m) for i in 1:6]
        end

        for T2 in (Int, Float32, Float64)
            wv = f(T2[1, 1, 1, 1, 1, 0, 0])
            x2 = collect(T, 2:6)
            for k in 1:6
                cumk = @inferred(cumulant(x, k, wv))
                @test cumk isa float(promote_type(T, T2))
                @test cumk ≈ cumulant(x2, k) rtol = cbrt(eps(typeof(cumk)))
            end
            cumks = @inferred(cumulant(x, 1:6, wv))
            @test cumks isa Vector{float(promote_type(T, T2))}
            @test cumks == [cumulant(x, i, wv) for i in 1:6]

            for TM in (Int, Float32, Float64)
                m = TM(3)
                for k in 1:6
                    cumk = @inferred(cumulant(x, k, wv, m))
                    @test cumk isa float(promote_type(T, T2, TM))
                    @test cumk ≈ cumulant(x2, k, m) rtol = cbrt(eps(typeof(cumk)))
                end
                cumks = @inferred(cumulant(x, 1:6, wv, m))
                @test cumks isa Vector{float(promote_type(T, T2, TM))}
                @test cumks == [cumulant(x, i, wv, m) for i in 1:6]
            end
        end

        # Invalid arguments
        @test_throws ArgumentError cumulant(x, -1)
        @test_throws ArgumentError cumulant(x, 0)
        @test_throws ArgumentError cumulant(x, 0:3)
        @test_throws ArgumentError cumulant(x, -1:3)
        @test_throws ArgumentError cumulant(x, 1:0)

        for n in (length(x) - 1, length(x) + 1), krange in (1, 1:3)
            @test_throws DimensionMismatch("Inconsistent array lengths.") cumulant(x, krange, f(ones(n)))
            @test_throws DimensionMismatch("Inconsistent array lengths.") cumulant(x, krange, f(ones(n)), 0.0)
        end
    end
end

end # @testset "StatsBase.Moments"
