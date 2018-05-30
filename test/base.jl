# This file contains code formerly part of Julia. License is MIT: https://julialang.org/license

using Test, Random, LinearAlgebra, SparseArrays

@testset "var & std" begin
    # edge case: empty vector
    # iterable; this has to throw for type stability
    @test_throws ArgumentError var(())
    @test_throws ArgumentError var((); corrected=false)
    @test_throws ArgumentError var((); mean=2)
    @test_throws ArgumentError var((); mean=2, corrected=false)
    # reduction
    @test isnan(var(Int[]))
    @test isnan(var(Int[]; corrected=false))
    @test isnan(var(Int[]; mean=2))
    @test isnan(var(Int[]; mean=2, corrected=false))
    # reduction across dimensions
    @test isequal(var(Int[], dims=1), [NaN])
    @test isequal(var(Int[], dims=1; corrected=false), [NaN])
    @test isequal(var(Int[], dims=1; mean=[2]), [NaN])
    @test isequal(var(Int[], dims=1; mean=[2], corrected=false), [NaN])

    # edge case: one-element vector
    # iterable
    @test isnan(@inferred(var((1,))))
    @test var((1,); corrected=false) === 0.0
    @test var((1,); mean=2) === Inf
    @test var((1,); mean=2, corrected=false) === 1.0
    # reduction
    @test isnan(@inferred(var([1])))
    @test var([1]; corrected=false) === 0.0
    @test var([1]; mean=2) === Inf
    @test var([1]; mean=2, corrected=false) === 1.0
    # reduction across dimensions
    @test isequal(@inferred(var([1], dims=1)), [NaN])
    @test var([1], dims=1; corrected=false) ≈ [0.0]
    @test var([1], dims=1; mean=[2]) ≈ [Inf]
    @test var([1], dims=1; mean=[2], corrected=false) ≈ [1.0]

    @test var(1:8) == 6.
    @test varm(1:8,1) == varm(Vector(1:8),1)
    @test isnan(varm(1:1,1))
    @test isnan(var(1:1))
    @test isnan(var(1:-1))

    @test @inferred(var(1.0:8.0)) == 6.
    @test varm(1.0:8.0,1.0) == varm(Vector(1.0:8.0),1)
    @test isnan(varm(1.0:1.0,1.0))
    @test isnan(var(1.0:1.0))
    @test isnan(var(1.0:-1.0))

    @test @inferred(var(1.0f0:8.0f0)) === 6.f0
    @test varm(1.0f0:8.0f0,1.0f0) == varm(Vector(1.0f0:8.0f0),1)
    @test isnan(varm(1.0f0:1.0f0,1.0f0))
    @test isnan(var(1.0f0:1.0f0))
    @test isnan(var(1.0f0:-1.0f0))

    @test varm([1,2,3], 2) ≈ 1.
    @test var([1,2,3]) ≈ 1.
    @test var([1,2,3]; corrected=false) ≈ 2.0/3
    @test var([1,2,3]; mean=0) ≈ 7.
    @test var([1,2,3]; mean=0, corrected=false) ≈ 14.0/3

    @test varm((1,2,3), 2) ≈ 1.
    @test var((1,2,3)) ≈ 1.
    @test var((1,2,3); corrected=false) ≈ 2.0/3
    @test var((1,2,3); mean=0) ≈ 7.
    @test var((1,2,3); mean=0, corrected=false) ≈ 14.0/3
    @test_throws ArgumentError var((1,2,3); mean=())

    @test var([1 2 3 4 5; 6 7 8 9 10], dims=2) ≈ [2.5 2.5]'
    @test var([1 2 3 4 5; 6 7 8 9 10], dims=2; corrected=false) ≈ [2.0 2.0]'

    @test var(collect(1:99), dims=1) ≈ [825]
    @test var(Matrix(transpose(collect(1:99))), dims=2) ≈ [825]

    @test stdm([1,2,3], 2) ≈ 1.
    @test std([1,2,3]) ≈ 1.
    @test std([1,2,3]; corrected=false) ≈ sqrt(2.0/3)
    @test std([1,2,3]; mean=0) ≈ sqrt(7.0)
    @test std([1,2,3]; mean=0, corrected=false) ≈ sqrt(14.0/3)

    @test stdm([1.0,2,3], 2) ≈ 1.
    @test std([1.0,2,3]) ≈ 1.
    @test std([1.0,2,3]; corrected=false) ≈ sqrt(2.0/3)
    @test std([1.0,2,3]; mean=0) ≈ sqrt(7.0)
    @test std([1.0,2,3]; mean=0, corrected=false) ≈ sqrt(14.0/3)

    @test std([1.0,2,3]; dims=1)[] ≈ 1.
    @test std([1.0,2,3]; dims=1, corrected=false)[] ≈ sqrt(2.0/3)
    @test std([1.0,2,3]; dims=1, mean=[0])[] ≈ sqrt(7.0)
    @test std([1.0,2,3]; dims=1, mean=[0], corrected=false)[] ≈ sqrt(14.0/3)

    @test stdm((1,2,3), 2) ≈ 1.
    @test std((1,2,3)) ≈ 1.
    @test std((1,2,3); corrected=false) ≈ sqrt(2.0/3)
    @test std((1,2,3); mean=0) ≈ sqrt(7.0)
    @test std((1,2,3); mean=0, corrected=false) ≈ sqrt(14.0/3)

    @test std([1 2 3 4 5; 6 7 8 9 10], dims=2) ≈ sqrt.([2.5 2.5]')
    @test std([1 2 3 4 5; 6 7 8 9 10], dims=2; corrected=false) ≈ sqrt.([2.0 2.0]')

    let A = ComplexF64[exp(i*im) for i in 1:10^4]
        @test varm(A, 0.) ≈ sum(map(abs2, A)) / (length(A) - 1)
        @test varm(A, mean(A)) ≈ var(A)
    end

    @test var([1//1, 2//1]) isa Rational{Int}
    @test var([1//1, 2//1], dims=1) isa Vector{Rational{Int}}

    @test std([1//1, 2//1]) isa Float64
    @test std([1//1, 2//1], dims=1) isa Vector{Float64}
end

function safe_cov(x, y, zm::Bool, cr::Bool)
    n = length(x)
    if !zm
        x = x .- mean(x)
        y = y .- mean(y)
    end
    dot(vec(x), vec(y)) / (n - Int(cr))
end
X = [1.0  5.0;
     2.0  4.0;
     3.0  6.0;
     4.0  2.0;
     5.0  1.0]
Y = [6.0  2.0;
     1.0  7.0;
     5.0  8.0;
     3.0  4.0;
     2.0  3.0]

@testset "covariance" begin
    for vd in [1, 2], zm in [true, false], cr in [true, false]
        # println("vd = $vd: zm = $zm, cr = $cr")
        if vd == 1
            k = size(X, 2)
            Cxx = zeros(k, k)
            Cxy = zeros(k, k)
            for i = 1:k, j = 1:k
                Cxx[i,j] = safe_cov(X[:,i], X[:,j], zm, cr)
                Cxy[i,j] = safe_cov(X[:,i], Y[:,j], zm, cr)
            end
            x1 = vec(X[:,1])
            y1 = vec(Y[:,1])
        else
            k = size(X, 1)
            Cxx = zeros(k, k)
            Cxy = zeros(k, k)
            for i = 1:k, j = 1:k
                Cxx[i,j] = safe_cov(X[i,:], X[j,:], zm, cr)
                Cxy[i,j] = safe_cov(X[i,:], Y[j,:], zm, cr)
            end
            x1 = vec(X[1,:])
            y1 = vec(Y[1,:])
        end

        c = zm ? StatsBase.covm(x1, 0, corrected=cr) :
                 cov(x1, corrected=cr)
        @test isa(c, Float64)
        @test c ≈ Cxx[1,1]
        @inferred cov(x1, corrected=cr)

        @test cov(X) == StatsBase.covm(X, mean(X, dims=1))
        C = zm ? StatsBase.covm(X, 0, vd, corrected=cr) :
                 cov(X, dims=vd, corrected=cr)
        @test size(C) == (k, k)
        @test C ≈ Cxx
        @inferred cov(X, dims=vd, corrected=cr)

        @test cov(x1, y1) == StatsBase.covm(x1, mean(x1), y1, mean(y1))
        c = zm ? StatsBase.covm(x1, 0, y1, 0, corrected=cr) :
                 cov(x1, y1, corrected=cr)
        @test isa(c, Float64)
        @test c ≈ Cxy[1,1]
        @inferred cov(x1, y1, corrected=cr)

        if vd == 1
            @test cov(x1, Y) == StatsBase.covm(x1, mean(x1), Y, mean(Y, dims=1))
        end
        C = zm ? StatsBase.covm(x1, 0, Y, 0, vd, corrected=cr) :
                 cov(x1, Y, dims=vd, corrected=cr)
        @test size(C) == (1, k)
        @test vec(C) ≈ Cxy[1,:]
        @inferred cov(x1, Y, dims=vd, corrected=cr)

        if vd == 1
            @test cov(X, y1) == StatsBase.covm(X, mean(X, dims=1), y1, mean(y1))
        end
        C = zm ? StatsBase.covm(X, 0, y1, 0, vd, corrected=cr) :
                 cov(X, y1, dims=vd, corrected=cr)
        @test size(C) == (k, 1)
        @test vec(C) ≈ Cxy[:,1]
        @inferred cov(X, y1, dims=vd, corrected=cr)

        @test cov(X, Y) == StatsBase.covm(X, mean(X, dims=1), Y, mean(Y, dims=1))
        C = zm ? StatsBase.covm(X, 0, Y, 0, vd, corrected=cr) :
                 cov(X, Y, dims=vd, corrected=cr)
        @test size(C) == (k, k)
        @test C ≈ Cxy
        @inferred cov(X, Y, dims=vd, corrected=cr)
    end
end

function safe_cor(x, y, zm::Bool)
    if !zm
        x = x .- mean(x)
        y = y .- mean(y)
    end
    x = vec(x)
    y = vec(y)
    dot(x, y) / (sqrt(dot(x, x)) * sqrt(dot(y, y)))
end
@testset "correlation" begin
    for vd in [1, 2], zm in [true, false]
        # println("vd = $vd: zm = $zm")
        if vd == 1
            k = size(X, 2)
            Cxx = zeros(k, k)
            Cxy = zeros(k, k)
            for i = 1:k, j = 1:k
                Cxx[i,j] = safe_cor(X[:,i], X[:,j], zm)
                Cxy[i,j] = safe_cor(X[:,i], Y[:,j], zm)
            end
            x1 = vec(X[:,1])
            y1 = vec(Y[:,1])
        else
            k = size(X, 1)
            Cxx = zeros(k, k)
            Cxy = zeros(k, k)
            for i = 1:k, j = 1:k
                Cxx[i,j] = safe_cor(X[i,:], X[j,:], zm)
                Cxy[i,j] = safe_cor(X[i,:], Y[j,:], zm)
            end
            x1 = vec(X[1,:])
            y1 = vec(Y[1,:])
        end

        c = zm ? StatsBase.corm(x1, 0) : cor(x1)
        @test isa(c, Float64)
        @test c ≈ Cxx[1,1]
        @inferred cor(x1)

        @test cor(X) == StatsBase.corm(X, mean(X, dims=1))
        C = zm ? StatsBase.corm(X, 0, vd) : cor(X, dims=vd)
        @test size(C) == (k, k)
        @test C ≈ Cxx
        @inferred cor(X, dims=vd)

        @test cor(x1, y1) == StatsBase.corm(x1, mean(x1), y1, mean(y1))
        c = zm ? StatsBase.corm(x1, 0, y1, 0) : cor(x1, y1)
        @test isa(c, Float64)
        @test c ≈ Cxy[1,1]
        @inferred cor(x1, y1)

        if vd == 1
            @test cor(x1, Y) == StatsBase.corm(x1, mean(x1), Y, mean(Y, dims=1))
        end
        C = zm ? StatsBase.corm(x1, 0, Y, 0, vd) : cor(x1, Y, dims=vd)
        @test size(C) == (1, k)
        @test vec(C) ≈ Cxy[1,:]
        @inferred cor(x1, Y, dims=vd)

        if vd == 1
            @test cor(X, y1) == StatsBase.corm(X, mean(X, dims=1), y1, mean(y1))
        end
        C = zm ? StatsBase.corm(X, 0, y1, 0, vd) : cor(X, y1, dims=vd)
        @test size(C) == (k, 1)
        @test vec(C) ≈ Cxy[:,1]
        @inferred cor(X, y1, dims=vd)

        @test cor(X, Y) == StatsBase.corm(X, mean(X, dims=1), Y, mean(Y, dims=1))
        C = zm ? StatsBase.corm(X, 0, Y, 0, vd) : cor(X, Y, dims=vd)
        @test size(C) == (k, k)
        @test C ≈ Cxy
        @inferred cor(X, Y, dims=vd)
    end

    @test cor(repeat(1:17, 1, 17))[2] <= 1.0
    @test cor(1:17, 1:17) <= 1.0
    @test cor(1:17, 18:34) <= 1.0
    let tmp = range(1, stop=85, length=100)
        tmp2 = Vector(tmp)
        @test cor(tmp, tmp) <= 1.0
        @test cor(tmp, tmp2) <= 1.0
    end
end

@testset "variance of complex arrays (#13309)" begin
    z = rand(ComplexF64, 10)
    @test var(z) ≈ invoke(var, Tuple{Any}, z) ≈ cov(z) ≈ var(z,dims=1)[1] ≈ sum(abs2, z .- mean(z))/9
    @test isa(var(z), Float64)
    @test isa(invoke(var, Tuple{Any}, z), Float64)
    @test isa(cov(z), Float64)
    @test isa(var(z,dims=1), Vector{Float64})
    @test varm(z, 0.0) ≈ invoke(varm, Tuple{Any,Float64}, z, 0.0) ≈ sum(abs2, z)/9
    @test isa(varm(z, 0.0), Float64)
    @test isa(invoke(varm, Tuple{Any,Float64}, z, 0.0), Float64)
    @test cor(z) === 1.0
    v = varm([1.0+2.0im], 0; corrected = false)
    @test v ≈ 5
    @test isa(v, Float64)
end

@testset "cov and cor of complex arrays (issue #21093)" begin
    x = [2.7 - 3.3im, 0.9 + 5.4im, 0.1 + 0.2im, -1.7 - 5.8im, 1.1 + 1.9im]
    y = [-1.7 - 1.6im, -0.2 + 6.5im, 0.8 - 10.0im, 9.1 - 3.4im, 2.7 - 5.5im]
    @test cov(x, y) ≈ 4.8365 - 12.119im
    @test cov(y, x) ≈ 4.8365 + 12.119im
    @test cov(x, reshape(y, :, 1)) ≈ reshape([4.8365 - 12.119im], 1, 1)
    @test cov(reshape(x, :, 1), y) ≈ reshape([4.8365 - 12.119im], 1, 1)
    @test cov(reshape(x, :, 1), reshape(y, :, 1)) ≈ reshape([4.8365 - 12.119im], 1, 1)
    @test cov([x y]) ≈ [21.779 4.8365-12.119im;
                        4.8365+12.119im 54.548]
    @test cor(x, y) ≈ 0.14032104449218274 - 0.35160772008699703im
    @test cor(y, x) ≈ 0.14032104449218274 + 0.35160772008699703im
    @test cor(x, reshape(y, :, 1)) ≈ reshape([0.14032104449218274 - 0.35160772008699703im], 1, 1)
    @test cor(reshape(x, :, 1), y) ≈ reshape([0.14032104449218274 - 0.35160772008699703im], 1, 1)
    @test cor(reshape(x, :, 1), reshape(y, :, 1)) ≈ reshape([0.14032104449218274 - 0.35160772008699703im], 1, 1)
    @test cor([x y]) ≈ [1.0                                          0.14032104449218274-0.35160772008699703im
                        0.14032104449218274+0.35160772008699703im  1.0]
end

@testset "Issue #17153 and PR #17154" begin
    a = rand(10,10)
    b = copy(a)
    x = var(a, dims=1)
    @test b == a
    x = var(a, dims=2)
    @test b == a
    x = std(a, dims=1)
    @test b == a
    x = std(a, dims=2)
    @test b == a
end

@testset "Unitful elements" begin
    r = Furlong(1):Furlong(1):Furlong(2)
    a = Vector(r)
    @test var(r) == var(a) == Furlong{2}(0.5)
    @test std(r) == std(a) == Furlong{1}(sqrt(0.5))

    # Issue #21786
    A = [Furlong{1}(rand(-5:5)) for i in 1:2, j in 1:2]
    @test var(A, dims=1)[1] === var(A[:, 1])
    @test std(A, dims=1)[1] === std(A[:, 1])
end

# Issue #22901
@testset "var and std of Any arrays" begin
    x = Any[1, 2, 4, 10]
    y = Any[1, 2, 4, 10//1]
    @test var(x) === 16.25
    @test var(y) === 65//4
    @test std(x) === sqrt(16.25)
end

@testset "Promotion in covzm. Issue #8080" begin
    A = [1 -1 -1; -1 1 1; -1 1 -1; 1 -1 -1; 1 -1 1]
    @test StatsBase.covzm(A) - mean(A, dims=1)'*mean(A, dims=1)*size(A, 1)/(size(A, 1) - 1) ≈ cov(A)
    A = [1//1 -1 -1; -1 1 1; -1 1 -1; 1 -1 -1; 1 -1 1]
    @test (A'A - size(A, 1)*Base.mean(A, dims=1)'*Base.mean(A, dims=1))/4 == cov(A)
end

@testset "cov/var/std of Vector{Vector}" begin
    x = [[2,4,6],[4,6,8]]
    @test var(x) ≈ vec(var([x[1] x[2]], dims=2))
    @test std(x) ≈ vec(std([x[1] x[2]], dims=2))
    @test cov(x) ≈ cov([x[1] x[2]], dims=2)
end

# Faster covariance function for sparse matrices
# Prevents densifying the input matrix when subtracting the mean
# Test against dense implementation
# PR https://github.com/JuliaLang/julia/pull/22735
# Part of this test needed to be hacked due to the treatment
# of Inf in sparse matrix algebra
# https://github.com/JuliaLang/julia/issues/22921
# The issue will be resolved in
# https://github.com/JuliaLang/julia/issues/22733
@testset "optimizing sparse $elty covariance" for elty in (Float64, Complex{Float64})
    n = 10
    p = 5
    np2 = div(n*p, 2)
    nzvals, x_sparse = guardsrand(1) do
        if elty <: Real
            nzvals = randn(np2)
        else
            nzvals = complex.(randn(np2), randn(np2))
        end
        nzvals, sparse(rand(1:n, np2), rand(1:p, np2), nzvals, n, p)
    end
    x_dense  = convert(Matrix{elty}, x_sparse)
    @testset "Test with no Infs and NaNs, vardim=$vardim, corrected=$corrected" for vardim in (1, 2),
                                                                                 corrected in (true, false)
        @test cov(x_sparse, dims=vardim, corrected=corrected) ≈
              cov(x_dense , dims=vardim, corrected=corrected)
    end

    @testset "Test with $x11, vardim=$vardim, corrected=$corrected" for x11 in (NaN, Inf),
                                                                     vardim in (1, 2),
                                                                  corrected in (true, false)
        x_sparse[1,1] = x11
        x_dense[1 ,1] = x11

        cov_sparse = cov(x_sparse, dims=vardim, corrected=corrected)
        cov_dense  = cov(x_dense , dims=vardim, corrected=corrected)
        @test cov_sparse[2:end, 2:end] ≈ cov_dense[2:end, 2:end]
        @test isfinite.(cov_sparse) == isfinite.(cov_dense)
        @test isfinite.(cov_sparse) == isfinite.(cov_dense)
    end

    @testset "Test with NaN and Inf, vardim=$vardim, corrected=$corrected" for vardim in (1, 2),
                                                                            corrected in (true, false)
        x_sparse[1,1] = Inf
        x_dense[1 ,1] = Inf
        x_sparse[2,1] = NaN
        x_dense[2 ,1] = NaN

        cov_sparse = cov(x_sparse, dims=vardim, corrected=corrected)
        cov_dense  = cov(x_dense , dims=vardim, corrected=corrected)
        @test cov_sparse[(1 + vardim):end, (1 + vardim):end] ≈
              cov_dense[ (1 + vardim):end, (1 + vardim):end]
        @test isfinite.(cov_sparse) == isfinite.(cov_dense)
        @test isfinite.(cov_sparse) == isfinite.(cov_dense)
    end
end

# issue #6672
@test std(AbstractFloat[1,2,3], dims=1) == [1.0]

@testset "var of sparse array" begin
    se33 = SparseMatrixCSC{Float64}(I, 3, 3)
    sA = sprandn(3, 7, 0.5)
    pA = sparse(rand(3, 7))

    for arr in (se33, sA, pA)
        farr = Array(arr)
        @test var(arr) ≈ var(farr)
        @test var(arr, dims=1) ≈ var(farr, dims=1)
        @test var(arr, dims=2) ≈ var(farr, dims=2)
        @test var(arr, dims=(1, 2)) ≈ [var(farr)]
        @test isequal(var(arr, dims=3), var(farr, dims=3))
    end

    @testset "empty cases" begin
        @test var(sparse(Int[])) === NaN
        @test isequal(var(spzeros(0, 1), dims=1), var(Matrix{Int}(I, 0, 1), dims=1))
        @test isequal(var(spzeros(0, 1), dims=2), var(Matrix{Int}(I, 0, 1), dims=2))
        @test isequal(var(spzeros(0, 1), dims=(1, 2)), var(Matrix{Int}(I, 0, 1), dims=(1, 2)))
        @test isequal(var(spzeros(0, 1), dims=3), var(Matrix{Int}(I, 0, 1), dims=3))
    end
end

@testset "var: empty cases" begin
    A = Matrix{Int}(undef, 0,1)
    @test var(A) === NaN

    @test isequal(var(A, dims=1), fill(NaN, 1, 1))
    @test isequal(var(A, dims=2), fill(NaN, 0, 1))
    @test isequal(var(A, dims=(1, 2)), fill(NaN, 1, 1))
    @test isequal(var(A, dims=3), fill(NaN, 0, 1))
end

@testset "linreg" begin
    # make sure unequal input arrays throw an error
    x = [2; 5; 6]
    y = [3; 7; 10; 10]
    @test_throws DimensionMismatch linreg(x, y)
    x = [2 5 6]
    y = [3; 7; 10]
    @test_throws MethodError linreg(x, y)

    # check (UnitRange, Array)
    x = 1:12
    y = [5.5; 6.3; 7.6; 8.8; 10.9; 11.79; 13.48; 15.02; 17.77; 20.81; 22.0; 22.99]
    @test [linreg(x,y)...] ≈ [2.5559090909090867, 1.6960139860139862]
    @test [linreg(view(x,1:6),view(y,1:6))...] ≈ [3.8366666666666642,1.3271428571428574]

    # check (LinRange, UnitRange)
    x = range(1.0, stop=12.0, length=100)
    y = -100:-1
    @test [linreg(x, y)...] ≈ [-109.0, 9.0]

    # check (UnitRange, UnitRange)
    x = 1:12
    y = 12:-1:1
    @test [linreg(x, y)...] ≈ [13.0, -1.0]

    # check (LinRange, LinRange)
    x = range(-5, stop=10, length=100)
    y = range(50, stop=200, length=100)
    @test [linreg(x, y)...] ≈ [100.0, 10.0]

    # check (Array, Array)
    # Anscombe's quartet (https://en.wikipedia.org/wiki/Anscombe%27s_quartet)
    x123 = [10.0; 8.0; 13.0; 9.0; 11.0; 14.0; 6.0; 4.0; 12.0; 7.0; 5.0]
    y1 = [8.04; 6.95; 7.58; 8.81; 8.33; 9.96; 7.24; 4.26; 10.84; 4.82; 5.68]
    @test [linreg(x123,y1)...] ≈ [3.0,0.5] atol=15e-5

    y2 = [9.14; 8.14; 8.74; 8.77; 9.26; 8.10; 6.12; 3.10; 9.13; 7.26; 4.74]
    @test [linreg(x123,y2)...] ≈ [3.0,0.5] atol=10e-3

    y3 = [7.46; 6.77; 12.74; 7.11; 7.81; 8.84; 6.08; 5.39; 8.15; 6.42; 5.73]
    @test [linreg(x123,y3)...] ≈ [3.0,0.5] atol=10e-3

    x4 = [8.0; 8.0; 8.0; 8.0; 8.0; 8.0; 8.0; 19.0; 8.0; 8.0; 8.0]
    y4 = [6.58; 5.76; 7.71; 8.84; 8.47; 7.04; 5.25; 12.50; 5.56; 7.91; 6.89]
    @test [linreg(x4,y4)...] ≈ [3.0,0.5] atol=10e-3
end
