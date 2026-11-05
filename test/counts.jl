using StatsBase
using Test
using OffsetArrays

n = 5000

@testset "abstract unit ranges as count levels" begin
    x = [1, 1, 2, 3]
    w = weights([1.0, 2.0, 3.0, 4.0])
    @test counts(x, Base.OneTo(4)) == [2, 1, 1, 0]
    @test counts(x, Base.OneTo(5)) == [2, 1, 1, 0, 0]
    @test isempty(counts(x, Base.OneTo(0)))

    lo = typemin(Int)
    limits = Base.IdentityUnitRange(lo:lo+1)
    wide = BigInt[lo, big(lo) + 1]
    @test counts(wide, limits) == [1, 1]
    @test counts(wide, limits, weights([2.0, 3.0])) == [2.0, 3.0]

    # Bin positions depend on the values, even for a range with offset axes.
    levels = Base.IdentityUnitRange(2:5)
    @test counts(x, levels) == [1, 1, 0, 0]
    @test addcounts!(ones(Int, 4), x, levels) == [2, 2, 1, 1]
    @test counts(x, levels, w) == [3.0, 4.0, 0.0, 0.0]
    @test proportions(x, levels) == [0.25, 0.25, 0.0, 0.0]
    @test proportions(x, levels, w) ≈ [0.3, 0.4, 0.0, 0.0]
end

@testset "integer count bin arithmetic" begin
    # The 201st bin must not wrap during Int8 subtraction.
    levels = Int8(-100):Int8(100)
    x = Int8[-101, -100, 100, 101]
    w = weights([1.0, 2.0, 3.0, 4.0])
    expected = zeros(Int, 201)
    expected[[1, 201]] .= 1
    expected_w = zeros(201)
    expected_w[[1, 201]] = [2.0, 3.0]
    @test counts(x, levels) == expected
    @test counts(x, levels, w) == expected_w
    @test addcounts!(ones(Int, 201), x, levels) == expected .+ 1
    @test addcounts!(ones(201), x, levels, w) == expected_w .+ 1
    @test proportions(x, levels) == expected ./ 4
    @test proportions(x, levels, w) ≈ expected_w ./ 10
    @test x == Int8[-101, -100, 100, 101]
    @test w == weights([1.0, 2.0, 3.0, 4.0])

    # Bool is an Integer but does not have a widen method.
    bools = Bool[false, true, true]
    @test counts(bools, 0:1) == [1, 2]
    @test counts(bools, 0:1, weights([2.0, 3.0, 5.0])) == [2.0, 8.0]
    @test addcounts!(ones(Int, 2), bools, 0:1) == [2, 3]
    @test addcounts!(ones(2), bools, 0:1, weights([2.0, 3.0, 5.0])) == [3.0, 9.0]

    # Use ordinal bins with mixed signedness and abstract eltypes.
    mixed = Integer[Int8(-100), UInt64(100), big(0)]
    expected[101] = 1
    @test counts(mixed, levels) == expected
    @test counts(mixed, levels, weights([2.0, 3.0, 4.0])) ==
        [i == 1 ? 2.0 : i == 101 ? 4.0 : i == 201 ? 3.0 : 0.0 for i in 1:201]
    unsigned = UInt64[0, 100, 101]
    expected_u = zeros(Int, 201)
    expected_u[[101, 201]] .= 1
    @test counts(unsigned, levels) == expected_u
    @test counts(unsigned, levels, weights([2.0, 3.0, 4.0])) ==
        [i == 101 ? 2.0 : i == 201 ? 3.0 : 0.0 for i in 1:201]

    # Convert both the range bound and input before the index arithmetic.
    wider = typemin(Int16):typemax(Int16)
    expected_wider = zeros(Int, 65536)
    expected_wider[32869] = 1
    @test counts(Int8[100], wider) == expected_wider
    @test counts(Int8[100], wider, weights([2.0])) == expected_wider .* 2

    for T in (UInt64, UInt128)
        hi = typemax(T)
        limits = (hi - one(T)):hi
        values = T[hi - one(T), hi]
        @test counts(values, limits) == [1, 1]
        @test counts(values, limits, weights([2.0, 3.0])) == [2.0, 3.0]
        @test counts(BigInt.(values), limits) == [1, 1]
        @test counts(BigInt.(values), limits, weights([2.0, 3.0])) == [2.0, 3.0]
        zero_limits = zero(T):one(T)
        @test counts(BigInt[0, 1], zero_limits) == [1, 1]
        @test counts(BigInt[0, 1], zero_limits, weights([2.0, 3.0])) == [2.0, 3.0]
    end

    huge = big(10)^50
    @test counts([huge, huge + 1], huge:huge+1) == [1, 1]
    @test counts([huge, huge + 1], huge:huge+1, weights([2.0, 3.0])) == [2.0, 3.0]
end

@testset "1D integer counts" begin
    x = rand(1:5, n)
    w = weights(rand(n))
    x0 = deepcopy(x)
    w0 = deepcopy(w)

    c0 = Int[count(v->v == i, x) for i in 1:5]
    @test counts(x, 5)                      == c0
    @test counts(x .+ 1, 2:6)               == c0
    @test proportions(x, 1:5)               ≈ (c0 ./ n)
    @test counts(reshape(x, 10, 50, 10), 5) == c0

    @test counts(x)                      == c0
    @test proportions(x)                 ≈ (c0 ./ n)
    @test counts(reshape(x, 10, 50, 10)) == c0

    c0 = reshape(c0, 1, 5)
    @test addcounts!(fill(0, 1, 5),                     x,  1:5) == c0
    @test addcounts!(fill(0, 1, 5), reshape(x, 10, 50, 10), 1:5) == c0

    c0 = Float64[sum(w.values[x .== i]) for i in 1:5]
    @test counts(x, 5, w)                      ≈ c0
    @test counts(x .+ 1, 2:6, w)               ≈ c0
    @test proportions(x, 1:5, w)               ≈ (c0 ./ sum(w))
    @test counts(reshape(x, 10, 50, 10), 5, w) ≈ c0 # Perhaps this should not be allowed

    @test counts(x, w)                      ≈ c0
    @test counts(x .+ 1, 2:6, w)            ≈ c0
    @test proportions(x, w)                 ≈ (c0 ./ sum(w))
    @test counts(reshape(x, 10, 50, 10), w) ≈ c0 # Perhaps this should not be allowed

    #addcounts! to row matrix
    c0 = reshape(c0, 1, 5)
    @test addcounts!(fill(0.0, 1, 5),                     x,  1:5, w) ≈ c0
    @test addcounts!(fill(0.0, 1, 5), reshape(x, 10, 50, 10), 1:5, w) ≈ c0 # Perhaps this should not be allowed

    @test x == x0
    @test w == w0
end


@testset "2D integer counts" begin
    x = rand(1:4, n)
    y = rand(1:5, n)
    w = weights(rand(n))
    x0 = deepcopy(x)
    y0 = deepcopy(y)
    w0 = deepcopy(w)

    c0 = Int[count(t->t != 0,  (x .== i) .& (y .== j)) for i in 1:4, j in 1:5]
    @test counts(x, y, (4, 5))                                           == c0
    @test counts(x .+ 2, y .+ 3, (3:6, 4:8))                             == c0
    @test proportions(x, y, (1:4, 1:5))                                  ≈ (c0 ./ n)
    @test counts(reshape(x, 10, 50, 10), reshape(y, 10, 50, 10), (4, 5)) == c0

    @test counts(x, y)                                           == c0
    @test counts(x .+ 2, y .+ 3, (3:6, 4:8))                     == c0
    @test proportions(x, y,)                                     ≈ (c0 ./ n)
    @test counts(reshape(x, 10, 50, 10), reshape(y, 10, 50, 10)) == c0

    c0 = Float64[sum(w.values[(x .== i) .& (y .== j)]) for i in 1:4, j in 1:5]
    @test counts(x, y, (4, 5), w)                                           ≈ c0
    @test counts(x .+ 2, y .+ 3, (3:6, 4:8), w)                             ≈ c0
    @test proportions(x, y, (1:4, 1:5), w)                                  ≈ (c0 ./ sum(w))
    @test counts(reshape(x, 10, 50, 10), reshape(y, 10, 50, 10), (4, 5), w) ≈ c0 # Perhaps this should not be allowed

    @test counts(x, y, w)                                           ≈ c0
    @test counts(x .+ 2, y .+ 3, (3:6, 4:8), w)                     ≈ c0
    @test proportions(x, y, w)                                      ≈ (c0 ./ sum(w))
    @test counts(reshape(x, 10, 50, 10), reshape(y, 10, 50, 10), w) ≈ c0 # Perhaps this should not be allowed

    @test x == x0
    @test y == y0
    @test w == w0
end

@testset "count map" begin
    x = ["a", "b", "a", "a", "b", "c"]
    w = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

    cm = countmap(x)
    @test cm["a"] == 3
    @test cm["b"] == 2
    @test cm["c"] == 1

    # iterator, non-radixsort
    cm_missing = countmap(skipmissing(x))
    cm_any_itr = countmap((i for i in x))
    @test cm_missing == cm_any_itr == cm
    @test cm_missing isa Dict{String, Int}
    @test cm_any_itr isa Dict{Any, Int}

    pm = proportionmap(x)
    @test pm["a"] ≈ (1/2)
    @test pm["b"] ≈ (1/3)
    @test pm["c"] ≈ (1/6)


    # testing the radixsort branch of countmap
    xx = repeat([6, 1, 3, 1], outer=100_000)
    cm = countmap(xx)
    @test cm == Dict(1 => 200_000, 3 => 100_000, 6 => 100_000)

    # with iterator
    cm_missing = countmap(skipmissing(xx))
    @test cm_missing isa Dict{Int, Int}
    @test cm_missing == cm

    cm_any_itr = countmap((i for i in xx))
    @test cm_any_itr isa Dict{Any,Int} # no knowledge about type
    @test cm_any_itr == cm

    # with multidimensional array
    @test countmap(reshape(xx, 20, 100, 20, 10); alg=:radixsort) == cm
    @test countmap(reshape(xx, 20, 100, 20, 10); alg=:dict)      == cm

    # with empty array
    @test countmap(Int[]) == Dict{Int, Int}()

    # testing the radixsort-based addcounts
    xx = repeat([6, 1, 3, 1], outer=100_000)
    cm = Dict{Int, Int}()
    StatsBase.addcounts_radixsort!(cm,xx)
    @test cm == Dict(1 => 200_000, 3 => 100_000, 6 => 100_000)
    xx2 = repeat([7, 1, 3, 1], outer=100_000)
    StatsBase.addcounts_radixsort!(cm,xx2)
    @test cm == Dict(1 => 400_000, 3 => 200_000, 6 => 100_000, 7 => 100_000)
    # with iterator
    cm_missing = Dict{Int, Int}()
    StatsBase.addcounts_radixsort!(cm_missing,skipmissing(xx))
    @test cm_missing == Dict(1 => 200_000, 3 => 100_000, 6 => 100_000)
    StatsBase.addcounts_radixsort!(cm_missing,skipmissing(xx2))
    @test cm_missing == Dict(1 => 400_000, 3 => 200_000, 6 => 100_000, 7 => 100_000)

    # testing the Dict-based addcounts
    cm = Dict{Int, Int}()
    cm_itr = Dict{Int, Int}()
    StatsBase.addcounts_dict!(cm,xx)
    StatsBase.addcounts_dict!(cm_itr,skipmissing(xx))
    @test cm_itr == cm == Dict(1 => 200_000, 3 => 100_000, 6 => 100_000)
    @test cm_itr isa Dict{Int, Int}

    cm = countmap(x, weights(w))
    @test cm["a"] == 5.5
    @test cm["b"] == 4.5
    @test cm["c"] == 3.5

    @test cm == countmap(x, w)
    @test cm == countmap(reshape(x, 2, 3), w)

    pm = proportionmap(x, weights(w))
    @test pm["a"] ≈ (5.5 / 13.5)
    @test pm["b"] ≈ (4.5 / 13.5)
    @test pm["c"] ≈ (3.5 / 13.5)

    # testing small bits type
    bx = [true, false, true, true, false]
    cm_bx_missing = countmap(skipmissing(bx))
    @test cm_bx_missing == countmap(bx) == Dict(true => 3, false => 2)
    @test cm_bx_missing isa Dict{Bool, Int}

    for T in [UInt8, UInt16, Int8, Int16]
        tx = T[typemin(T), 8, typemax(T), 19, 8]
        tx_missing = skipmissing(T[typemin(T), 8, typemax(T), 19, 8])
        cm_tx_missing = countmap(tx_missing)
        @test cm_tx_missing == countmap(tx) == Dict(typemin(T) => 1, typemax(T) => 1, 8 => 2, 19 => 1)
        @test cm_tx_missing isa Dict{T, Int}
    end

    # -0.0 and NaN
    @test countmap([0.0, -0.0, 0.0, -0.0, -0.0], alg=:dict) ==
        countmap([0.0, -0.0, 0.0, -0.0, -0.0], alg=:radixsort) ==
        Dict(0.0 => 2, -0.0 => 3)
    @test countmap([NaN, NaN], alg=:dict) ==
        countmap([NaN, NaN], alg=:radixsort) ==
        Dict(NaN => 2)
end

@testset "views" begin
    X = view([1,1,1,2,2], 1:5)
    @test countmap(X) == countmap(copy(X))
end

@testset "offset arrays" begin
    x = rand(1:5, n)
    w = rand(n)
    xw = weights(w)
    y = OffsetArray(x, n÷2)
    yw = weights(OffsetArray(w, n÷2))
    z = OffsetArray(x, -2n)
    zw = weights(OffsetArray(w, -2n))

    # proportions calls counts which calls addcounts!
    @test proportions(x)       == proportions(y)       == proportions(z)
    @test proportions(x, xw)   == proportions(y, yw)   == proportions(z, zw)
    @test proportionmap(x)     == proportionmap(y)     == proportionmap(z)
    @test proportionmap(x, xw) == proportionmap(y, yw) == proportionmap(z, zw)
    @test countmap(x) == countmap(x; alg = :dict) == countmap(x; alg = :radixsort) ==
          countmap(y) == countmap(y; alg = :dict) == countmap(y; alg = :radixsort) ==
          countmap(z) == countmap(z; alg = :dict) == countmap(z; alg = :radixsort)
    @test proportionmap(x, xw) == proportionmap(y, yw) == proportionmap(z, zw)
    # countmap and proportionmap only support the :dict algorithm for weighted sums.
end
