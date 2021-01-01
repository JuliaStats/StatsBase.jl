using StatsBase
using Test, Random, Statistics, LinearAlgebra
using Missings

const ≅ = isequal

Random.seed!(1)

# to avoid using specialized method
arbitrary_fun(x, y) = cor(x, y)

@testset "pairwise with $f" for f in (arbitrary_fun, cor, cov)
    @testset "basic interface" begin
        x = [rand(10) for _ in 1:4]
        y = [rand(Float32, 10) for _ in 1:5]
        # to test case where inference of returned eltype fails
        z = [Vector{Any}(rand(Float32, 10)) for _ in 1:5]

        res = @inferred pairwise(f, x, y)
        @test res isa Matrix{Float64}
        @test res == [f(xi, yi) for xi in x, yi in y]

        res = pairwise(f, y, z)
        @test res isa Matrix{Float32}
        @test res == [f(yi, zi) for yi in y, zi in z]

        res = pairwise(f, Any[[1.0, 2.0, 3.0], [1.0f0, 3.0f0, 10.5f0]])
        @test res isa Matrix{AbstractFloat}
        @test res == [f(xi, yi) for xi in ([1.0, 2.0, 3.0], [1.0f0, 3.0f0, 10.5f0]),
                                    yi in ([1.0, 2.0, 3.0], [1.0f0, 3.0f0, 10.5f0])]
        @test typeof.(res) == [Float64 Float64
                               Float64 Float32]

        @inferred pairwise(f, x, y)

        @test_throws ArgumentError pairwise(f, [Int[]], [Int[]])
    end

    @testset "missing values handling interface" begin
        xm = [ifelse.(rand(100) .> 0.9, missing, rand(100)) for _ in 1:4]
        ym = [ifelse.(rand(100) .> 0.9, missing, rand(Float32, 100)) for _ in 1:4]
        zm = [ifelse.(rand(100) .> 0.9, missing, rand(Float32, 100)) for _ in 1:4]

        res = pairwise(f, xm, ym)
        @test res isa Matrix{Missing}
        @test res ≅ [missing for xi in xm, yi in ym]

        res = pairwise(f, xm, ym, skipmissing=:pairwise)
        @test res isa Matrix{Float64}
        @test isapprox(res, [f(collect.(skipmissings(xi, yi))...) for xi in xm, yi in ym],
                       rtol=1e-6)

        res = pairwise(f, ym, zm, skipmissing=:pairwise)
        @test res isa Matrix{Float32}
        @test isapprox(res, [f(collect.(skipmissings(yi, zi))...) for yi in ym, zi in zm],
                       rtol=1e-6)

        nminds = mapreduce(x -> .!ismissing.(x),
                        (x, y) -> x .& y,
                        [xm; ym])
        res = pairwise(f, xm, ym, skipmissing=:listwise)
        @test res isa Matrix{Float64}
        @test isapprox(res, [f(view(xi, nminds), view(yi, nminds)) for xi in xm, yi in ym],
                       rtol=1e-6)

        if VERSION >= v"1.6.0-DEV"
            # inference of cor fails so use an inferrable function
            # to check that pairwise itself is inferrable
            for skipmissing in (:none, :pairwise, :listwise)
                g(x, y=x) = pairwise((x, y) -> x[1] * y[1], x, y, skipmissing=skipmissing)
                @test Core.Compiler.return_type(g, Tuple{Vector{Vector{Union{Float64, Missing}}}}) ==
                    Core.Compiler.return_type(g, Tuple{Vector{Vector{Union{Float64, Missing}}},
                                                    Vector{Vector{Union{Float64, Missing}}}}) ==
                    Matrix{<: Union{Float64, Missing}}
                if skipmissing in (:pairwise, :listwise)
                    @test_broken Core.Compiler.return_type(g, Tuple{Vector{Vector{Union{Float64, Missing}}}}) ==
                        Core.Compiler.return_type(g, Tuple{Vector{Vector{Union{Float64, Missing}}},
                                                        Vector{Vector{Union{Float64, Missing}}}}) ==
                        Matrix{Float64}
                end
            end
        end

        @test_throws ArgumentError pairwise(f, xm, ym, skipmissing=:something)

        # variable with only missings
        xm = [fill(missing, 10), rand(10)]
        ym = [rand(10), rand(10)]

        res = pairwise(f, xm, ym)
        @test res isa Matrix{Union{Float64, Missing}}
        @test res ≅ [f(xi, yi) for xi in xm, yi in ym]

        if VERSION >= v"1.5" # Fails with UndefVarError on Julia 1.0
            @test_throws ArgumentError pairwise(f, xm, ym, skipmissing=:pairwise)
            @test_throws ArgumentError pairwise(f, xm, ym, skipmissing=:listwise)
        end
    end

    @testset "iterators" begin
        x = (v for v in [rand(10) for _ in 1:4])
        y = (v for v in [rand(10) for _ in 1:4])

        @test pairwise(f, x, y) == pairwise(f, collect(x), collect(y))
        @test pairwise(f, x) == pairwise(f, collect(x))
    end

    @testset "two-argument method" begin
        x = [rand(10) for _ in 1:4]
        @test pairwise(f, x) == pairwise(f, x, x)
    end

    @testset "symmetric" begin
        x = [rand(10) for _ in 1:4]
        y = [rand(10) for _ in 1:4]
        @test pairwise(f, x, x, symmetric=true) ==
            pairwise(f, x, symmetric=true) ==
            Symmetric(pairwise(f, x, x), :U)
        @test_throws ArgumentError pairwise(f, x, y, symmetric=true)
    end

    @testset "cor corner cases" begin
        # Integer inputs must give a Float64 output
        res = pairwise(cor, [[1, 2, 3], [1, 5, 2]])
        @test res isa Matrix{Float64}
        @test res == [cor(xi, yi) for xi in ([1, 2, 3], [1, 5, 2]),
                                      yi in ([1, 2, 3], [1, 5, 2])]

        # NaNs are ignored for the diagonal
        res = pairwise(cor, [[1, 2, NaN], [1, 5, 2]])
        @test res isa Matrix{Float64}
        @test res ≅ [1.0 NaN
                     NaN 1.0]

        # missings are propagated even for the diagonal
        res = pairwise(cor, [[1, 2, 7], [1, 5, missing]])
        @test res isa Matrix{Union{Float64, Missing}}
        @test res ≅ [1.0 missing
                     missing missing]

        for sm in (:pairwise, :listwise)
            res = pairwise(cor, [[1, 2, NaN, 4], [1, 5, 5, missing]], skipmissing=sm)
            @test res isa Matrix{Float64}
            @test res ≅ [1.0 NaN
                        NaN 1.0]
        end
    end
end