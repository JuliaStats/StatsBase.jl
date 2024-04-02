using StatsBase, Random
using Test

@testset "rank_correlation_auxiliary_fns" begin

    #Auxiliary functions for corkendall
    n = 100
    x = [1, 2, 3, missing, 4]
    y = [missing, 1, 2, 3, 4]
    u = [missing, missing, 1, 2]
    v = [3, 4, missing, missing]

    mx = [1 2
        missing 3
        4 missing
        missing missing
        5 6]

    @test StatsBase.handle_pairwise(x, y) == ([2, 3, 4], [1, 2, 4])
    @test StatsBase.handle_pairwise(float.(x), y) == ([2.0, 3.0, 4.0], [1, 2, 4])
    @test StatsBase.handle_pairwise(x, float.(y)) == ([2, 3, 4], [1.0, 2.0, 4.0])
    @test StatsBase.handle_pairwise(u, v) == (Int64[], Int64[])

    v = collect(100:-1:1)
    StatsBase.insertion_sort!(v, 1, n)
    @test v == 1:n

    v = collect(1000:-1:1)
    StatsBase.merge_sort!(v, 1, 1000)
    @test v == 1:1000

    @test StatsBase.midpoint(1, 10) == 5
    @test StatsBase.midpoint(1, widen(10)) == 5

    @test StatsBase.equal_sum_subsets(0, 1) == Vector{Int64}[]
    @test sum.(StatsBase.equal_sum_subsets(100, 5)) == repeat([1010], 5)
    @test sort(vcat(StatsBase.equal_sum_subsets(500, 7)...)) == collect(1:500)

end

#=
Edge cases (to refer to when writing tests)
===========================================
In the symmetric case (x===y) on-diagonal terms should be 1.0 apart from the special case of
eltype(x) === Missing && skipmissing == :none, in which case the on-diagonal terms should be
missing.
Otherwise, if x and y are equal-length vectors with length(x)<2
f(x,y) == NaN
otherwise if x and y are equal-length vectors and either contains missing
f(x,y) == missing
otherwise if x and y are equal-length vectors and either contains NaN
f(x,y) == NaN
otherwise
f(x,y) = the required Kendall or Spearman correlation.

This behaviour is close to but not quite identical to Statistics.cor. e.g.:

julia> cor(fill(missing,3,2))
2×2 Matrix{Missing}:
 missing  missing
 missing  missing

julia> corkendall(fill(missing,3,2)) #SAME behavior as cor
2×2 Matrix{Missing}:
 missing  missing
 missing  missing

julia> cor(Matrix{Union{Missing,Float64}}(missing,5,3))
3×3 Matrix{Missing}:
 missing  missing  missing
 missing  missing  missing
 missing  missing  missing

julia> corkendall(Matrix{Union{Missing,Float64}}(missing,5,3)) #DIFFERENT behaviour from cor
3×3 Matrix{Union{Missing, Float64}}:
 1.0        missing   missing
  missing  1.0        missing
  missing   missing  1.0
=#

@testset "$f" for f in (corkendall, corspearman)

    n = 100
    big_n = 78_000
    X = Float64[1 0; 2 1; 3 0; 4 1; 5 10]
    Y = Float64[5 5 6; 3 4 1; 4 0 4; 2 6 1; 5 7 10]
    Xm = [1 0; missing 1; 2 1; 3 0; 4 1; 5 10]
    Ym = [5 5 6; 1 2 3; 3 4 1; 4 0 4; 2 6 1; 5 7 10]
    xm = [missing, missing, missing, missing, missing]
    xmm = hcat(xm, xm)
    a = [5, 2, 3, 4, 1]
    b = [1, 4, 2, 3, 5]

    x1 = X[:, 1]
    x2 = X[:, 2]
    y1 = Y[:, 1]

    c11 = f(x1, x1)
    c12 = f(x1, x2)
    c22 = f(x2, x2)

    # Test some known results
    if f == corkendall
        # Vector, Vector
        @test f(x1, y1) == -1 / sqrt(90)
        @test f(x2, y1) == -1 / sqrt(72)
        # Matrix, Vector
        @test f(X, y1) == [-1 / sqrt(90), -1 / sqrt(72)]
        # Vector, Matrix
        @test f(y1, X) == [-1 / sqrt(90) -1 / sqrt(72)]

        # big_n = 78_000 tests for overflow errors on 32 bit
        # Testing for overflow errors on 64 bit is not practical
        # This also tests merge_sort! since big_n is (much) greater than SMALL_THRESHOLD.
        # Test with many repeats
        @test f(repeat(x1, big_n), repeat(y1, big_n)) ≈ -1 / sqrt(90)
        @test f(repeat(x2, big_n), repeat(y1, big_n)) ≈ -1 / sqrt(72)
        @test f(repeat(X, big_n), repeat(y1, big_n)) ≈ [-1 / sqrt(90), -1 / sqrt(72)]
        @test f(repeat(y1, big_n), repeat(X, big_n)) ≈ [-1 / sqrt(90) -1 / sqrt(72)]
    elseif f == corspearman
        @test corspearman(x1, y1) ≈ -1 / sqrt(95)
        @test corspearman(repeat(x1, 1000), repeat(y1, 1000)) ≈ -1 / sqrt(95)
        @test corspearman(vcat(missing, x1, missing), vcat(missing, y1, missing), skipmissing=:pairwise) ≈ -1 / sqrt(95)
        @test corspearman(x2, y1) ≈ -3 / sqrt(1368)
        @test corspearman(X, y1) ≈ [-1 / sqrt(95), -3 / sqrt(1368)]
        @test corspearman(y1, X) ≈ [-1 / sqrt(95) -3 / sqrt(1368)]
    end

    # Matrix, Matrix
    @test f(X, X) ≈ [c11 c12; c12 c22]
    # Matrix
    @test f(X) ≈ [c11 c12; c12 c22]

    @test c11 == 1.0
    @test c22 == 1.0
    if f == corkendall
        @test c12 == 3 / sqrt(20)
    elseif f == corspearman
        @test c12 == 7 / sqrt(90)
    end
    @test f(repeat(X, n), repeat(X, n)) ≈ [c11 c12; c12 c22]
    @test f(repeat(X, n)) ≈ [c11 c12; c12 c22]
    @test f(X, Y) ≈
          [f(X[:, i], Y[:, j]) for i in axes(X, 2), j in axes(Y, 2)]

    @test f(vcat(missing, a), vcat(missing, b), skipmissing=:pairwise) == f(a, b)
    @test f(vcat(a, missing), vcat(missing, b), skipmissing=:pairwise) ==
          f(a[2:end], b[1:(end-1)])
    @test f(hcat(vcat(a, missing), vcat(missing, b)), skipmissing=:listwise) ==
          f(hcat(a[2:end], b[1:(end-1)]))
    @test f(Xm, Xm, skipmissing=:pairwise) == f(X, X)
    @test f(Xm, Xm, skipmissing=:listwise) == f(X, X)
    @test f(Xm, Ym, skipmissing=:listwise) == f(X, Y)
    if f == corkendall
        @test f(Xm, Ym, skipmissing=:pairwise) ≈ [-1/√90 0.4 1/√90; -2/√154 7/√165 -1/√154]
    end
    @test isnan(f([1, 2, 3, 4, 5], xm, skipmissing=:pairwise))
    @test isnan(f(xm, [1, 2, 3, 4, 5], skipmissing=:pairwise))
    @test isequal(f(xmm, skipmissing=:pairwise), [1.0 NaN; NaN 1.0])
    @test isequal(f(xmm, skipmissing=:none), [missing missing; missing missing])
    @test isequal(f(xmm, xmm, skipmissing=:none), [missing missing; missing missing])
    @test isequal(f(xmm, copy(xmm), skipmissing=:none),
        [missing missing; missing missing])
    @test isequal(f(xmm, xmm, skipmissing=:listwise), [1.0 NaN; NaN 1.0])
    @test isequal(f(xmm, copy(xmm), skipmissing=:listwise), [NaN NaN; NaN NaN])
    @test isequal(f(xmm, copy(xmm), skipmissing=:pairwise), [NaN NaN; NaN NaN])
    @test ismissing(f([1, 2, 3, 4, 5], xm, skipmissing=:none))
    @test ismissing(f([1, 2, 3, 4, 5], xm, skipmissing=:none))
    @test isequal(f(xmm, skipmissing=:none), [missing missing; missing missing])
    @test isequal(f(xmm, copy(xmm), skipmissing=:none),
        [missing missing; missing missing])
    @test isequal(f(hcat(Y, xm), skipmissing=:none), vcat(hcat(f(Y, skipmissing=:none),
            [missing, missing, missing]), [missing missing missing 1.0]))
    @test_throws ArgumentError f([1, 2, 3, 4], [4, 3, 2, 1], skipmissing=:listwise)

    # All elements identical should yield NaN
    @test isnan(f(repeat([1], n), collect(1:n)))
    @test f(repeat([0, 1, 1, 0], n), repeat([1, 0, 1, 0], n)) == 0.0

    # Test with no repeats
    @test f(collect(1:n), collect(1:n)) == 1.0
    @test f(collect(1:n), reverse(collect(1:n))) == -1.0

    # Inf and -Inf work ok
    @test f([-Inf, -0.0, Inf], [1, 2, 3]) == 1.0

    #Interaction of NaN and missing with skipmissing argument
    nan_and_missing = hcat(fill(NaN, 10), fill(missing, 10),
        vcat(fill(NaN, 5), fill(missing, 5)))
    @test isequal(f(nan_and_missing, skipmissing=:none),
        [1.0 missing missing; missing 1.0 missing; missing missing 1.0])
    @test isequal(f(nan_and_missing, copy(nan_and_missing), skipmissing=:none),
        [NaN missing missing; missing missing missing; missing missing missing])
    @test isequal(f(nan_and_missing, skipmissing=:pairwise),
        [1.0 NaN NaN; NaN 1.0 NaN; NaN NaN 1.0])
    @test isequal(f(nan_and_missing, copy(nan_and_missing), skipmissing=:pairwise),
        fill(NaN, 3, 3))
    @test isequal(f(nan_and_missing, skipmissing=:listwise),
        [1.0 NaN NaN; NaN 1.0 NaN; NaN NaN 1.0])
    @test isequal(f(nan_and_missing, copy(nan_and_missing), skipmissing=:listwise),
        fill(NaN, 3, 3))

    #Reject nonsense skipmissing argument
    @test_throws ArgumentError f(X; skipmissing=:foo)
    @test_throws ArgumentError f(Xm; skipmissing=:foo)

    # Inputs have fewer than 2 rows
    @test isequal(f([], []), NaN)
    @test isequal(f(fill(1, 0, 2), fill(1, 0, 2)), [NaN NaN; NaN NaN])
    @test isequal(f(fill(1, 0, 2)), [1.0 NaN; NaN 1.0])
    @test isequal(f([1;;], [1;;]), [NaN;;])
    @test isequal(f([1;;]), [1.0;;])
    @test isequal(f([missing], [missing]), NaN)
    @test isequal(f([1], [1]), NaN)
    @test isequal(f([NaN], [NaN]), NaN)
    @test isequal(f([1], [1]), NaN)
    @test isequal(f([]), 1.0)
    @test isequal(f([1]), 1.0)
    @test isequal(f([NaN]), 1.0)
    @test isequal(f([missing]), missing)
    @test isequal(f([missing]), missing)
    @test isequal(f([missing], [missing missing]), [NaN NaN])
    @test isequal(f([missing missing]), [missing NaN; NaN missing])
    @test isequal(f([missing missing], [missing missing]), [NaN NaN; NaN NaN])

    # Exercise catch block in method _cor (when f === corspearman)
    @test isequal(f(vcat([1 2 3], fill(missing, 2, 3))),
        [1.0 missing missing; missing 1.0 missing; missing missing 1.0])
    # Exercise "correction" of on-diagonal terms in method
    # _pairwise!(::Val{:none}, f::typeof(corspearman),...
    @test isequal(f([missing missing; 1 1; 1 1]), [1.0 missing; missing 1.0])

    # Works for not-numbers
    @test isequal(f(["a", "b", "c"], ["a", "b", "c"]), 1.0)
    @test isequal(f(["a", "b", "c"], ["c", "b", "a"]), -1.0)
    @test (f(["a" "z"; "b" "y"; "c" "x"]) ≈ [1.0 -1.0; -1.0 1.0])
    @test (f(["a" 3; "b" 2; "c" 1]) ≈ [1.0 -1.0; -1.0 1.0])

    #Works for zero size input ( [;;] not compatible with Julia 1.0.5)
    let nada = Array{Any,2}(undef, 0, 0)
        @test isequal(f(nada), nada)
        @test isequal(f(nada, nada), nada)
        @test isequal(f(nada, nada, skipmissing=:pairwise), nada)
        @test isequal(f(nada, nada, skipmissing=:listwise), nada)
    end

    # Wrong dimensions
    @test_throws DimensionMismatch f([1], [1, 2])
    @test_throws DimensionMismatch f([1], [1 2; 3 4])
    @test_throws DimensionMismatch f([1 2; 3 4], [1])
    @test_throws DimensionMismatch f([1 2; 3 4; 4 6], [1 2; 3 4])

    # All eight three-element permutations, for these cases corspearman and corkendall
    # have the same return values
    z = [1 1 1;
        1 1 2;
        1 2 2;
        1 2 2;
        1 2 1;
        2 1 2;
        1 1 2;
        2 2 2]

    @test f(z) ≈ [1 0 1/3; 0 1 0; 1/3 0 1]
    @test f(z, z) ≈ [1 0 1/3; 0 1 0; 1/3 0 1]
    @test f(z[:, 1], z) ≈ [1 0 1 / 3]
    @test f(z, z[:, 1]) ≈ [1; 0; 1 / 3]

    z = float(z)
    @test f(z) ≈ [1 0 1/3; 0 1 0; 1/3 0 1]
    @test f(z, z) ≈ [1 0 1/3; 0 1 0; 1/3 0 1]
    @test f(z[:, 1], z) ≈ [1 0 1 / 3]
    @test f(z, z[:, 1]) ≈ [1; 0; 1 / 3]

    w = repeat(z, n)
    @test f(w) ≈ [1 0 1/3; 0 1 0; 1/3 0 1]
    @test f(w, w) ≈ [1 0 1/3; 0 1 0; 1/3 0 1]
    @test f(w[:, 1], w) ≈ [1 0 1 / 3]
    @test f(w, w[:, 1]) ≈ [1; 0; 1 / 3]

    # NaN handling
    Xnan = copy(X)
    Xnan[1, 1] = NaN
    Ynan = copy(Y)
    Ynan[2, 1] = NaN

    @test isnan(f([1.0, NaN, 2.0], [2.0, 1.0, 3.4]))
    @test all(isnan, f([1.0, NaN], [1 2; 3 4]))
    @test all(isnan, f([1 2; 3 4], [1.0, NaN]))
    @test isequal(f([1 NaN; NaN 4]), [1 NaN; NaN 1])
    @test all(isnan, f([1 NaN; NaN 4], [1 NaN; NaN 4]))
    @test all(isnan, f([1 NaN; NaN 4], [NaN 1; NaN 4]))

    @test isequal(f(Xnan, Ynan),
        [f(Xnan[:, i], Ynan[:, j]) for i in axes(Xnan, 2), j in axes(Ynan, 2)])
    @test isequal(f(Xnan),
        [i == j ? 1.0 : f(Xnan[:, i], Xnan[:, j])
         for i in axes(Xnan, 2), j in axes(Xnan, 2)])
    for k in 1:2
        @test isequal(f(Xnan[:, k], Ynan),
            [f(Xnan[:, k], Ynan[:, j]) for i in 1:1, j in axes(Ynan, 2)])
        # TODO: fix corkendall (PR#659)
        if f === corspearman
            @test isequal(f(Xnan, Ynan[:, k]),
                [f(Xnan[:, i], Ynan[:, k]) for i in axes(Xnan, 2), j in 1:1])
        else
            @test isequal(f(Xnan, Ynan[:, k]),
                [f(Xnan[:, i], Ynan[:, k]) for i in axes(Xnan, 2)])
        end
    end

end

# Don't artificially boost coverage stats when checking for mutation and allocation size
# COV_EXCL_START
@testset "Check no mutation in $f" for f in (corkendall, corspearman)

    nr = 50
    nc = 5
    cutoff = min(0.1, 10 / (nr * nc))
    X = randn(nr, nc)
    x = randn(nr)
    Xm = ifelse.(rand(nr, nc) .< cutoff, missing, randn(nr, nc))
    xm = ifelse.(rand(nr) .< cutoff, missing, randn(nr))
    Y = randn(nr, nc)
    y = randn(nr)
    Ym = ifelse.(rand(nr, nc) .< cutoff, missing, randn(nr, nc))
    ym = ifelse.(rand(nr) .< cutoff, missing, randn(nr))

    for arg1 in (X, x, Xm, xm)
        for arg2 in (Y, y, Ym, ym)
            for skipmissing in (:pairwise, :none, :listwise)
                for f in (corkendall, corspearman)
                    if !((arg1 isa Vector) && (arg2 isa Vector) && skipmissing == :listwise)
                        copy1 = copy(arg1)
                        copy2 = copy(arg2)
                        res = f(arg1, arg2; skipmissing)
                        @test isequal(arg1, copy1)
                        @test isequal(arg2, copy2)
                    end
                end
            end
        end
    end
end

@testset "corkendall and corspearman allocations" begin

    Random.seed!(1)
    x = rand(1000, 10)
    xm = ifelse.(x .< 0.1, missing, x)
    #compile
    corkendall(x)
    corkendall(xm, skipmissing=:listwise)
    corkendall(xm, skipmissing=:pairwise)
    corspearman(x)
    corspearman(xm, skipmissing=:listwise)
    corspearman(xm, skipmissing=:pairwise)
    x = rand(1000, 100)
    xm = ifelse.(x .< 0.01, missing, x)

    #=When executing code such as corkendall(x) allocations are approximately affine in the
    number of threads, thanks to use of task-local storage. The tests below have a "safety
    factor" of 1.2 against the expected size of allocations.
    =#
    @test (@allocated corkendall(x)) < (897_808 + Threads.nthreads() * 58_044) * 1.2
    @test (@allocated corkendall(xm, skipmissing=:listwise)) < (1_119_392 + Threads.nthreads() * 22_172) * 1.2
    @test (@allocated corkendall(xm, skipmissing=:pairwise)) < (892_112 + Threads.nthreads() * 61_116) * 1.2
    @test (@allocated corspearman(x)) < (2_678_448 + Threads.nthreads() * 9_128) * 1.2
    @test (@allocated corspearman(xm, skipmissing=:listwise)) < (1_803_712 + Threads.nthreads() * 3_992) * 1.2
    @test (@allocated corspearman(xm, skipmissing=:pairwise)) < (1_692_208 + Threads.nthreads() * 67_172) * 1.2

end
# COV_EXCL_STOP