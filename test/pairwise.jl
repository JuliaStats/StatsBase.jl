using StatsBase
using Test, Random, Statistics, LinearAlgebra
using Missings

const ≅ = isequal

Random.seed!(1)

# to avoid using specialized method
arbitrary_fun(x, y) = cor(x, y)

@testset "pairwise and pairwise! with $f" for f in (arbitrary_fun, cor, cov)
    @testset "basic interface" begin
        x = [rand(10) for _ in 1:4]
        y = [rand(Float32, 10) for _ in 1:5]
        # to test case where inference of returned eltype fails
        z = [Vector{Any}(rand(Float32, 10)) for _ in 1:5]

        res = @inferred pairwise(f, x, y)
        @test res isa Matrix{Float64}
        res2 = zeros(Float64, size(res))
        @test pairwise!(f, res2, x, y) === res2
        @test res == res2 == [f(xi, yi) for xi in x, yi in y]

        res = pairwise(f, y, z)
        @test res isa Matrix{Float32}
        res2 = zeros(Float32, size(res))
        @test pairwise!(f, res2, y, z) === res2
        @test res == res2 == [f(yi, zi) for yi in y, zi in z]

        res = pairwise(f, Any[[1.0, 2.0, 3.0], [1.0f0, 3.0f0, 10.5f0]])
        @test res isa Matrix{Float64}
        res2 = zeros(AbstractFloat, size(res))
        @test pairwise!(f, res2, Any[[1.0, 2.0, 3.0], [1.0f0, 3.0f0, 10.5f0]]) === res2
        @test res == res2 ==
            [f(xi, yi) for xi in ([1.0, 2.0, 3.0], [1.0f0, 3.0f0, 10.5f0]),
                           yi in ([1.0, 2.0, 3.0], [1.0f0, 3.0f0, 10.5f0])]
        @test res isa Matrix{Float64}

        @inferred pairwise(f, x, y)

        @test_throws Union{ArgumentError,MethodError} pairwise(f, [Int[]], [Int[]])
        @test_throws Union{ArgumentError,MethodError} pairwise!(f, zeros(1, 1), [Int[]], [Int[]])

        res = pairwise(f, [], [])
        @test size(res) == (0, 0)
        @test res isa Matrix{Any}
        res2 = zeros(0, 0)
        @test pairwise!(f, res2, [], []) === res2

        res = pairwise(f, Vector{Int}[], Vector{Int}[])
        @test size(res) == (0, 0)
        @test res isa Matrix{Float64}
        res2 = zeros(0, 0)
        @test pairwise!(f, res2, Vector{Int}[], Vector{Int}[]) === res2

        res = pairwise(f, [[1, 2]], Vector{Int}[])
        @test size(res) == (1, 0)
        @test res isa Matrix{Float64}
        res2 = zeros(1, 0)
        @test pairwise!(f, res2, [[1, 2]], Vector{Int}[]) === res2

        res = pairwise(f, Vector{Int}[], [[1, 2], [2, 3]])
        @test size(res) == (0, 2)
        @test res isa Matrix{Float64}
        res2 = zeros(0, 2)
        @test pairwise!(f, res2, [], [[1, 2], [2, 3]]) === res2

        @test_throws DimensionMismatch pairwise!(f, zeros(1, 2), x, y)
        @test_throws DimensionMismatch pairwise!(f, zeros(1, 2), [], [])
        @test_throws DimensionMismatch pairwise!(f, zeros(0, 0),
                                                 [], [[1, 2], [2, 3]])
    end

    @testset "missing values handling interface" begin
        xm = [ifelse.(rand(100) .> 0.9, missing, rand(100)) for _ in 1:4]
        ym = [ifelse.(rand(100) .> 0.9, missing, rand(Float32, 100)) for _ in 1:4]
        zm = [ifelse.(rand(100) .> 0.9, missing, rand(Float32, 100)) for _ in 1:4]

        res = pairwise(f, xm, ym)
        @test res isa Matrix{Missing}
        res2 = zeros(Union{Float64, Missing}, size(res))
        @test pairwise!(f, res2, xm, ym) === res2
        @test res ≅ res2 ≅ [missing for xi in xm, yi in ym]

        res = pairwise(f, xm, ym, skipmissing=:pairwise)
        @test res isa Matrix{Float64}
        res2 = zeros(Union{Float64, Missing}, size(res))
        @test pairwise!(f, res2, xm, ym, skipmissing=:pairwise) === res2
        @test res ≅ res2
        @test isapprox(res, [f(collect.(skipmissings(xi, yi))...) for xi in xm, yi in ym],
                       rtol=1e-6)

        res = pairwise(f, ym, zm, skipmissing=:pairwise)
        @test res isa Matrix{Float32}
        res2 = zeros(Union{Float32, Missing}, size(res))
        @test pairwise!(f, res2, ym, zm, skipmissing=:pairwise) === res2
        @test res ≅ res2
        @test isapprox(res, [f(collect.(skipmissings(yi, zi))...) for yi in ym, zi in zm],
                       rtol=1e-6)

        nminds = mapreduce(x -> .!ismissing.(x),
                        (x, y) -> x .& y,
                        [xm; ym])
        res = pairwise(f, xm, ym, skipmissing=:listwise)
        @test res isa Matrix{Float64}
        res2 = zeros(Union{Float64, Missing}, size(res))
        @test pairwise!(f, res2, xm, ym, skipmissing=:listwise) === res2
        @test res ≅ res2
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
        @test_throws ArgumentError pairwise!(f, zeros(Union{Float64, Missing},
                                                      length(xm), length(ym)), xm, ym,
                                             skipmissing=:something)

        # variable with only missings
        xm = [fill(missing, 10), rand(10)]
        ym = [rand(10), rand(10)]

        res = pairwise(f, xm, ym)
        @test res isa Matrix{Union{Float64, Missing}}
        res2 = zeros(Union{Float64, Missing}, size(res))
        @test pairwise!(f, res2, xm, ym) === res2
        @test res ≅ res2 ≅ [f(xi, yi) for xi in xm, yi in ym]

        if VERSION >= v"1.5" # Fails with UndefVarError on Julia 1.0
            @test_throws Union{ArgumentError,MethodError} pairwise(f, xm, ym, skipmissing=:pairwise)
            @test_throws Union{ArgumentError,MethodError} pairwise(f, xm, ym, skipmissing=:listwise)

            res = zeros(Union{Float64, Missing}, length(xm), length(ym))
            @test_throws Union{ArgumentError,MethodError} pairwise!(f, res, xm, ym, skipmissing=:pairwise)
            @test_throws Union{ArgumentError,MethodError} pairwise!(f, res, xm, ym, skipmissing=:listwise)
        end

        for sm in (:pairwise, :listwise)
            @test_throws ArgumentError pairwise(f, [[1, 2]], [1], skipmissing=sm)
            @test_throws ArgumentError pairwise(f, [1], [[1, 2]], skipmissing=sm)
            @test_throws ArgumentError pairwise(f, [1], [1], skipmissing=sm)
        end
    end

    @testset "iterators" begin
        x = (v for v in [rand(10) for _ in 1:4])
        y = (v for v in [rand(10) for _ in 1:4])

        res = @inferred pairwise(f, x, y)
        res2 = zeros(size(res))
        @test pairwise!(f, res2, x, y) === res2
        @test res == res2 == pairwise(f, collect(x), collect(y))

        res = @inferred(pairwise(f, x))
        res2 = zeros(size(res))
        @test pairwise!(f, res2, x) === res2
        @test res == res2 == pairwise(f, collect(x))
    end

    @testset "non-vector entries" begin
        x = (Iterators.drop(v, 1) for v in [rand(10) for _ in 1:4])
        y = (Iterators.drop(v, 1) for v in [rand(10) for _ in 1:4])

        @test pairwise((x, y) -> f(collect(x), collect(y)), x, y) ==
            [f(collect(xi), collect(yi)) for xi in x, yi in y]
        @test pairwise((x, y) -> f(collect(x), collect(y)), x) ==
            [f(collect(xi1), collect(xi2)) for xi1 in x, xi2 in x]
        @test_throws ArgumentError pairwise((x, y) -> f(collect(x), collect(y)), x, y,
                                            skipmissing=:pairwise)
        @test_throws ArgumentError pairwise((x, y) -> f(collect(x), collect(y)), x, y,
                                            skipmissing=:listwise)
    end

    @testset "two-argument method" begin
        x = [rand(10) for _ in 1:4]
        res = pairwise(f, x)
        res2 = zeros(size(res))
        @test pairwise!(f, res2, x) === res2
        @test res == res2 == pairwise(f, x, x)
    end

    @testset "symmetric" begin
        x = [rand(10) for _ in 1:4]
        y = [rand(10) for _ in 1:4]

        @test pairwise(f, x, x, symmetric=true) ==
            pairwise(f, x, symmetric=true) ==
            Symmetric(pairwise(f, x, x), :U)

        res = zeros(4, 4)
        res2 = zeros(4, 4)
        @test pairwise!(f, res, x, x, symmetric=true) === res
        @test pairwise!(f, res2, x, symmetric=true) === res2
        @test res == res2 == Symmetric(pairwise(f, x, x), :U)

        @test_throws ArgumentError pairwise(f, x, y, symmetric=true)
        @test_throws ArgumentError pairwise!(f, res, x, y, symmetric=true)
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

        # missings are ignored for the diagonal
        res = pairwise(cor, [[1, 2, 7], [1, 5, missing]])
        @test res isa Matrix{Union{Float64, Missing}}
        @test res ≅ [1.0 missing
                     missing 1.0]
        res = pairwise(cor, Vector{Union{Int, Missing}}[[missing, missing, missing],
                                                        [missing, missing, missing]])
        @test res isa Matrix{Union{Float64, Missing}}
        @test res ≅ [1.0 missing
                     missing 1.0]
        if VERSION >= v"1.5"
            # except when eltype is Missing
            res = pairwise(cor, [[missing, missing, missing],
                                 [missing, missing, missing]])
            @test res isa Matrix{Missing}
            @test res ≅ [missing missing
                         missing missing]
        end

        for sm in (:pairwise, :listwise)
            res = pairwise(cor, [[1, 2, NaN, 4], [1, 5, 5, missing]], skipmissing=sm)
            @test res isa Matrix{Float64}
            @test res ≅ [1.0 NaN
                         NaN 1.0]
            if VERSION >= v"1.5"
                @test_throws ArgumentError pairwise(cor, [[missing, missing, missing],
                                                          [missing, missing, missing]],
                                                    skipmissing=sm)
            end
        end
    end

    @testset "promote_type_union" begin
        @test StatsBase.promote_type_union(Int) === Int
        @test StatsBase.promote_type_union(Real) === Real
        @test StatsBase.promote_type_union(Union{Int, Float64}) === Float64
        @test StatsBase.promote_type_union(Union{Int, Missing}) === Union{Int, Missing}
        @test StatsBase.promote_type_union(Union{Int, String}) === Any
        @test StatsBase.promote_type_union(Vector) === Any
        @test StatsBase.promote_type_union(Union{}) === Union{}
        if VERSION >= v"1.6.0-DEV"
            @test StatsBase.promote_type_union(Tuple{Union{Int, Float64}}) ===
                Tuple{Real}
        else
            @test StatsBase.promote_type_union(Tuple{Union{Int, Float64}}) ===
                Any
        end
    end

    @testset "type-unstable corner case (#771)" begin
        v = [rand(5) for _=1:10]
        function f(v)
            pairwise(v) do x, y
                (x[1] < 0 ? nothing :
                    x[1] > y[1] ? 1 : 1.5,
                 0)
            end
        end
        res = f(v)
        @test res isa Matrix{Tuple{Real, Int}}
    end
end