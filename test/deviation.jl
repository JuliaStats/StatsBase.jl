using StatsBase, OffsetArrays
using Test

@testset "counting (arrays with element types $T1 and $T2)" for T1 in (Int, Float32, Float64), T2 in (Int, Float32, Float64)
    a = T1[1, 2, 3, 4, 5, 6, 7]
    b = T2[1, 3, 3, 4, 6, 7, 8]
    a_offset = OffsetArray(a, -5:1)
    b_offset = OffsetArray(b, -5:1)
    for (a, b) in ((a, b), (a_offset, b_offset))
        @test @inferred(counteq(a, b))::Int == 3
        @test @inferred(countne(a, b))::Int == 4
    end

    # Empty arrays
    empty_a = T1[]
    empty_b = T2[]
    empty_a_offset = OffsetArray(empty_a, -5)
    empty_b_offset = OffsetArray(empty_b, -5)
    for (a, b) in ((empty_a, empty_b), (empty_a_offset, empty_b_offset))
        @test @inferred(counteq(a, b))::Int == 0
        @test @inferred(countne(a, b))::Int == 0
    end

    # Inconsistent lengths
    err = DimensionMismatch("Inconsistent array lengths.")
    for (a, b) in ((a, empty_b), (empty_a, b), (a_offset, empty_b_offset), (empty_a_offset, b_offset))
        @test_throws err counteq(a, b)
        @test_throws err countne(a, b)
    end
end

@testset "deviation (arrays with element types $T1 and $T2)" for T1 in (Float32, Float64), T2 in (Float32, Float64)
    T = promote_type(T1, T2)
    a = rand(T1, 5, 6)
    b = rand(T2, 5, 6)
    a_offset = OffsetArray(a, 5, -10)
    b_offset = OffsetArray(b, 5, -10)
    for (a, b) in ((a, b), (a_offset, b_offset))
        @test @inferred(sqL2dist(a, b))::T ≈ sum(abs2.(a - b))
        @test @inferred(L2dist(a, b))::T ≈ sqrt(sqL2dist(a, b))
        @test @inferred(L1dist(a, b))::T ≈ sum(abs.(a - b))
        @test @inferred(Linfdist(a, b))::T ≈ maximum(abs.(a - b))
        @test @inferred(gkldiv(a, b))::T ≈ sum(a .* log.(a ./ b) - a + b)
        @test @inferred(meanad(a, b))::T ≈ mean(abs.(a - b))
        @test @inferred(maxad(a, b))::T ≈ maximum(abs.(a - b))
        @test @inferred(msd(a, b))::T ≈ mean(abs2.(a - b))
        @test @inferred(rmsd(a, b))::T ≈ sqrt(msd(a, b))
        @test @inferred(rmsd(a, b; normalize = true))::T ≈ rmsd(a, b) / (maximum(a) - minimum(a))
        for T2 in (Int, Float32, Float64)
            S = promote_type(T, T2)
            @test @inferred(psnr(a, b, T2(2)))::S ≈ 10 * log10(4 / msd(a, b))
        end
    end

    # Empty arrays
    empty_a = T1[]
    empty_b = T2[]
    empty_a_offset = OffsetArray(empty_a, 5)
    empty_b_offset = OffsetArray(empty_b, 5)
    for (a, b) in ((empty_a, empty_b), (empty_a_offset, empty_b_offset))
        @test iszero(@inferred(sqL2dist(a, b))::T)
        @test iszero(@inferred(L2dist(a, b))::T)
        @test iszero(@inferred(L1dist(a, b))::T)
        @test iszero(@inferred(Linfdist(a, b))::T)
        @test iszero(@inferred(gkldiv(a, b))::T)
        @test isnan(@inferred(meanad(a, b))::T)
        @test iszero(@inferred(maxad(a, b))::T)
        @test isnan(@inferred(msd(a, b))::T)
        @test isnan(@inferred(rmsd(a, b))::T)
        @test isnan(@inferred(rmsd(a, b; normalize = true))::T)
        for T2 in (Int, Float32, Float64)
            S = promote_type(T, T2)
            @test isnan(@inferred(psnr(a, b, T2(2)))::S)
        end
    end

    err = DimensionMismatch("Inconsistent array lengths.")
    for (a, b) in ((a, empty_b), (empty_a, b), (a_offset, empty_b_offset), (empty_a_offset, b_offset))
        @test_throws err sqL2dist(a, b)
        @test_throws err L2dist(a, b)
        @test_throws err L1dist(a, b)
        @test_throws err Linfdist(a, b)
        @test_throws err gkldiv(a, b)
        @test_throws err meanad(a, b)
        @test_throws err maxad(a, b)
        @test_throws err msd(a, b)
        @test_throws err rmsd(a, b)
        @test_throws err rmsd(a, b; normalize = true)
        for T2 in (Int, Float32, Float64)
            @test_throws err psnr(a, b, T2(2))
        end
    end
end
