using StatsBase
using LinearAlgebra, Random, Test

@testset "StatsBase.Histogram" begin


@testset "Histogram binindex and binvolume" begin
    edg1 = -2:0.5:9
    edg1f0 = -2:0.5f0:9
    edg2 = [-2, -1, 2, 7, 19]
    h1 = Histogram(edg1)
    h2 = Histogram((edg1, edg2))
    h3 = Histogram((edg1f0, edg2))

    @test h1 == Histogram(edg1, :left, false)

    @test @inferred StatsBase.binindex(h1, -0.5) == 4
    @test @inferred StatsBase.binindex(h2, (1.5, 2)) == (8, 3)

    @test [StatsBase.binvolume(h1, i) for i in axes(h1.weights, 1)] ≈ diff(edg1)
    @test [StatsBase.binvolume(h2, (i,j)) for i in axes(h2.weights, 1), j in axes(h2.weights, 2)] ≈ diff(edg1) * diff(edg2)'

    @test typeof(@inferred(StatsBase.binvolume(h2, (1,1)))) == Float64
    @test typeof(@inferred(StatsBase.binvolume(h3, (1,1)))) == Float32
    @test typeof(@inferred(StatsBase.binvolume(Float64, h3, (1,1)))) == Float64
end


@testset "Histogram append" begin
    h = Histogram(0:20:100, Float64, :left, false)
    @test @inferred(append!(h, 0:0.5:99.99)) == h
    @test append!(Histogram(0:20:100, Float64, :left, false), 0:0.5:99.99).weights ≈ [40,40,40,40,40]
    @test append!(Histogram(0:20:100, Float64, :left, true), 0:0.5:99.99).weights ≈ [2,2,2,2,2]
    @test append!(Histogram(0:20:100, Float64, :left, false), 0:0.5:99.99, fill(2, 200)).weights ≈ [80,80,80,80,80]
    @test append!(Histogram(0:20:100, Float64, :left, true), 0:0.5:99.99, fill(2, 200)).weights ≈ [4,4,4,4,4]
end


@testset "Histogram fit" begin
    @test sum(fit(Histogram,[1,2,3]).weights) == 3
    @test fit(Histogram,Int[]).weights == Int[]
    @test fit(Histogram,[1]).weights == [1]
    @test fit(Histogram,[1,2,3],[0,2,4]) == Histogram([0,2,4],[1,2], :left)
    @test fit(Histogram,[1,2,3],[0,2,4]) != Histogram([0,2,4],[1,1], :left)
    @test fit(Histogram,[1,2,3],0:2:4) == Histogram(0:2:4,[1,2], :left)
    @test all(fit(Histogram,[0:99;]/100,0.0:0.01:1.0).weights .==1)
    @test fit(Histogram,[1,1,1,1,1]).weights[1] == 5
    @test sum(fit(Histogram,(rand(100),rand(100))).weights) == 100
    @test fit(Histogram,1:100,nbins=5,closed=:right).weights == [20,20,20,20,20]
    @test fit(Histogram,1:100,nbins=5,closed=:left).weights == [19,20,20,20,20,1]
    @test fit(Histogram,0:99,nbins=5,closed=:right).weights == [1,20,20,20,20,19]
    @test fit(Histogram,0:99,nbins=5,closed=:left).weights == [20,20,20,20,20]

    @test fit(Histogram,(0:99,0:99),nbins=5).weights == Matrix(Diagonal([20,20,20,20,20]))
    @test fit(Histogram,(0:99,0:99),nbins=(5,5)).weights == Matrix(Diagonal([20,20,20,20,20]))

    @test fit(Histogram,0:99,weights(ones(100)),nbins=5).weights == [20,20,20,20,20]
    @test fit(Histogram,0:99,weights(2*ones(100)),nbins=5).weights == [40,40,40,40,40]
    @test fit(Histogram{Int32},0:99,weights(2*ones(100)),nbins=5).weights == [40,40,40,40,40]
    @test fit(Histogram{Float32},0:99,weights(2*ones(100)),nbins=5).weights == [40,40,40,40,40]

    d = collect(0:99)
    v = view(d, fill(true, 100))
    @test fit(Histogram{Float32},v,weights(2*ones(100)),nbins=5).weights == [40,40,40,40,40]

    # Issue 972: OutOfMemoryError for identical values of large magnitude
    @test sum(fit(Histogram, [-5.603325961434038e25, -5.603325961434038e25]).weights) == 2
end


@testset "Histogram element type" begin
    @test eltype(@inferred(fit(Histogram,1:100,weights(ones(Int,100)),nbins=5)).weights) == Int
    @test eltype(@inferred(fit(Histogram{Float32},1:100,weights(ones(Int,100)),nbins=5)).weights) == Float32
    @test eltype(@inferred(fit(Histogram,1:100,weights(ones(Float64,100)),nbins=5)).weights) == Float64
    @test eltype(@inferred(fit(Histogram{Float32},1:100,weights(ones(Float64,100)),nbins=5)).weights) == Float32
end


@testset "histrange" begin
    # Note: atm histrange must be qualified
    @test @inferred(StatsBase.histrange(Float64[], 0, :left)) == [0.0]
    @test StatsBase.histrange(Float64[1:5;], 1, :left) == 0.0:5.0:10.0
    @test StatsBase.histrange(Float64[1:10;], 1, :left) == 0.0:10.0:20.0
    @test StatsBase.histrange(1.0, 10.0, 1, :left) == 0.0:10.0:20.0

    @test StatsBase.histrange([0.201,0.299], 10, :left) == 0.2:0.01:0.3
    @test StatsBase.histrange([0.2,0.299], 10, :left) == 0.2:0.01:0.3
    @test StatsBase.histrange([0.2,0.3], 10, :left)  == 0.2:0.01:0.31
    @test StatsBase.histrange(0.2, 0.3,  10, :left)  == 0.2:0.01:0.31
    @test StatsBase.histrange([0.2,0.3], 10, :right) == 0.19:0.01:0.3
    @test StatsBase.histrange(0.2, 0.3,  10, :right) == 0.19:0.01:0.3

    @test StatsBase.histrange([200.1,299.9], 10, :left) == 200.0:10.0:300.0
    @test StatsBase.histrange([200.0,299.9], 10, :left) == 200.0:10.0:300.0
    @test StatsBase.histrange([200.0,300.0], 10, :left) == 200.0:10.0:310.0
    @test StatsBase.histrange([200.0,300.0], 10, :right) == 190.0:10.0:300.0

    @test @inferred(StatsBase.histrange(Int64[1:5;], 1, :left)) == 0:5:10
    @test StatsBase.histrange([1:5;], 1, :left) isa StatsBase.UniformEdges{Float64}
    @test step(StatsBase.histrange([1:5;], 1, :left)) == 5.0
    @test step(StatsBase.histrange([0.2, 0.3], 10, :left)) == 0.01
    @test step(StatsBase.histrange([0.0, 0.2], 9, :right)) == 0.05
    @test StatsBase.histrange(Int64[1:10;], 1, :left) == 0:10:20

    @test StatsBase.histrange([0, 1, 2, 3], 4, :left) == 0.0:1.0:4.0
    @test StatsBase.histrange([0, 1, 1, 3], 4, :left) == 0.0:1.0:4.0
    @test StatsBase.histrange([0, 9], 4, :left) == 0.0:5.0:10.0
    @test StatsBase.histrange([0, 19], 4, :left) == 0.0:5.0:20.0
    @test StatsBase.histrange([0, 599], 4, :left) == 0.0:200.0:600.0
    @test StatsBase.histrange([-1, -1000], 4, :left) == -1000.0:500.0:0.0

    # Base issue #13326
    l,h = extrema(StatsBase.histrange([typemin(Int),typemax(Int)], 10, :left))
    @test l <= typemin(Int)
    @test h >= typemax(Int)

    # Issue 616/667
    @test StatsBase.histrange([1.0 for i in 1:100], 10, :left) == 1.0:1.0:2.0
    @test StatsBase.histrange([1.05 for i in 1:100], 10, :left) == [1.0, 2.0]
    @test StatsBase.histrange([1.05 for i in 1:100], 10, :right) == [1.0, 2.0]
    @test StatsBase.histrange([1.0 for i in 1:100], 10, :right) == [0.0, 1.0]
    @test StatsBase.histrange([0.0], 3, :left) == [0.0, 1.0]
    @test StatsBase.histrange([0.0], 3, :right) == [-1.0, 0.0]
    @test StatsBase.histrange([-0.7], 3, :left) == [-1.0, 0.0]

    # Edges have the floating point type of the data
    @test StatsBase.histrange(Float32[0.7, 0.8], 12, :left) isa StatsBase.UniformEdges{Float32}
    @test StatsBase.histrange(Float16[0.7, 0.8], 12, :left) isa StatsBase.UniformEdges{Float16}
    @test StatsBase.histrange(BigFloat[0.7, 0.8], 12, :left) isa StatsBase.UniformEdges{BigFloat}
    @test step(StatsBase.histrange(Float32[0.7, 0.8], 12, :left)) === 0.01f0
    @test step(StatsBase.histrange(Float16[0.001, 0.002], 12, :left)) === Float16(0.0001)
    @test StatsBase.histrange(Float32[0.7, 0.8], 12, :left) == Float32.(0.7:0.01:0.81)
    @test StatsBase.histrange(Float32[0.7, 0.8], 12, :right) == Float32.(0.69:0.01:0.8)
    @test StatsBase.histrange(Float16[0.001, 0.002], 12, :left) == Float16.(0.001:0.0001:0.0021)
    # BigFloat[0.7, 0.8] would convert the Float64 literals, which lie below the decimals
    @test StatsBase.histrange(BigFloat.(["0.7", "0.8"]), 12, :left) ==
        BigFloat.(["0.7", "0.71", "0.72", "0.73", "0.74", "0.75", "0.76", "0.77", "0.78", "0.79", "0.8", "0.81"])

    # Beyond the exactly representable powers of ten, edges are within two ulps of the decimal
    for (v, expected) in (([1e300, 3e300], [1e300, 1.5e300, 2e300, 2.5e300, 3e300, 3.5e300]),
                          ([1e-300, 3e-300], [1e-300, 1.5e-300, 2e-300, 2.5e-300, 3e-300, 3.5e-300]),
                          (Float32[1e-33, 3e-33], Float32[1e-33, 1.5e-33, 2e-33, 2.5e-33, 3e-33, 3.5e-33]),
                          (Float16[1e-4, 3e-4], Float16[1e-4, 1.5e-4, 2e-4, 2.5e-4, 3e-4, 3.5e-4]))
        r = StatsBase.histrange(v, 4, :left)
        @test length(r) == length(expected)
        @test all(abs(x - y) <= 2 * eps(y) for (x, y) in zip(r, expected))
    end
    @test StatsBase.histrange(Float16[1e-4, 3e-4], 4, :left) == Float16[1e-4, 1.5e-4, 2e-4, 2.5e-4, 3e-4, 3.5e-4]

    # Requested bins finer than the resolution of the data give fewer, strictly increasing edges
    r = StatsBase.histrange([1e17, 1e17 + 16], 4, :right)
    @test issorted(r, lt = <=) && first(r) < 1e17 && 1e17 + 16 <= last(r)
    @test length(r) == 3  # 2 bins rather than the 4 requested: the data spans one ulp
    for F in (Float16, Float32, Float64), closed in (:left, :right), n in (1, 3, 10, 1000)
        for x in (F(0.7), F(1), prevfloat(F(2)), F(-1000), floatmax(F) / 2), k in (0, 1, 2, 5, 40)
            y = nextfloat(x, k)
            r = StatsBase.histrange(x, y, n, closed)
            @test issorted(r, lt = <=)
            @test length(r) <= max(k, 1) + 2
            if closed == :right
                @test first(r) < x && y <= last(r)
            else
                @test first(r) <= x && y < last(r)
            end
        end
    end

    # binindex on UniformEdges agrees with the generic search for all inputs
    for closed in (:left, :right), r in (StatsBase.histrange([0.2, 0.3], 10, closed),
                                         StatsBase.histrange(Float32[0.0, 1.0], 7, closed),
                                         StatsBase.histrange([-3.0, 5.0], 1, closed))
        edges = collect(r)
        xs = vcat(edges, prevfloat.(edges), nextfloat.(edges), rand(eltype(r), 200) .* (last(r) - first(r) + 1) .+ first(r) .- 0.5,
                  eltype(r)[-0.0, 0.0, Inf, -Inf], [1, 0, -1])
        for x in xs
            @test StatsBase._edge_binindex(r, closed, x) == StatsBase._edge_binindex(edges, closed, x)
        end
        # NaN is outside the bins on both paths (which side differs, as in Base's searchsorted)
        @test StatsBase._edge_binindex(r, closed, NaN) ∉ 1:length(r) - 1
        @test StatsBase._edge_binindex(edges, closed, NaN) ∉ 1:length(r) - 1
        h = fit(Histogram, [first(r)], r, closed=closed)
        @test all(StatsBase.binvolume(h, i) == step(r) for i in 1:length(r) - 1)
    end

    @test_throws ArgumentError StatsBase.histrange([1.0, Inf], 4, :left)
    @test_throws ArgumentError StatsBase.histrange([NaN, 1.0], 4, :left)

    # Issue #1009: observations at the extremes were dropped because the endpoint checks
    # were done in the element type while the edges were Float64, or because the last
    # element of the edge range evaluated an ulp below the value that was checked
    for T in (Float16, Float32, Float64, BigFloat), closed in (:left, :right), nbins in 1:12
        for v in ([0.7, 0.8], [0.81, 0.87], [0.0, 0.2], [8.2, 8.93], [6.34, 6.44],
                  [2.5, 3.6, 4.4], [0.001, 0.002], [0.1, 0.1, 1.0, 0.3, 0.6], [0.3, 0.4])
            vT = T.(v)
            r = StatsBase.histrange(vT, nbins, closed)
            @test eltype(r) == T
            lo, hi = extrema(vT)
            if closed == :right
                @test first(r) < lo && hi <= last(r)
            else
                @test first(r) <= lo && hi < last(r)
            end
            h = fit(Histogram, vT, nbins=nbins, closed=closed)
            @test sum(h.weights) == length(vT)
        end
    end

    # Issue 972: values whose magnitude is so large that a bin width of
    # one (or the computed width) is below the floating-point spacing
    let x = -5.603325961434038e25
        for closed in (:left, :right)
            r = StatsBase.histrange([x, x], 10, closed)
            @test length(r) == 2
            @test r[1] < r[2]
            @test closed == :left ? (r[1] <= x < r[2]) : (r[1] < x <= r[2])
        end
        for closed in (:left, :right), n in (10, 100, 1000)
            r = StatsBase.histrange(x, nextfloat(x, 3), n, closed)
            @test length(r) < 100
            @test minimum(diff(r)) >= eps(x)
            if closed == :left
                @test first(r) <= x
                @test last(r) > nextfloat(x, 3)
            else
                @test first(r) < x
                @test last(r) >= nextfloat(x, 3)
            end
        end
    end
    # tiny spans relative to the magnitude used to stall in the fractional
    # bin width branch
    let r = StatsBase.histrange(1.0e15, 1.0e15 + 0.001, 10, :left)
        @test length(r) < 100
        @test first(r) <= 1.0e15
        @test last(r) > 1.0e15 + 0.001
    end

    @test_throws ArgumentError StatsBase.histrange([1, 10], 0, :left)
    @test_throws ArgumentError StatsBase.histrange([1, 10], -1, :left)
    @test_throws ArgumentError StatsBase.histrange([1.0, 10.0], 0, :left)
    @test_throws ArgumentError StatsBase.histrange([1.0, 10.0], -1, :left)
    @test_throws ArgumentError StatsBase.histrange(Float64[],-1, :left)
    @test_throws ArgumentError StatsBase.histrange([0.], 0, :left)
end


@testset "Histogram show" begin
    # hist show
    show_h = sprint(show, fit(Histogram,[0,1,2]))
    @test occursin("edges:\n  0.0:1.0:3.0", show_h)
    # uniform edges print as first edge, width, last edge however many there are
    show_big = sprint(show, fit(Histogram, rand(1000), 0.0:0.001:1.0))
    @test occursin("edges:\n  0.0:0.001:1.0", show_big)
    @test sprint(show, StatsBase.histrange([0.0, 1.0], 1000, :left)) == "0.0:0.001:1.001"
    # and user-supplied vectors are abbreviated
    show_vec = sprint(show, fit(Histogram, rand(1000), collect(0.0:0.001:1.0)))
    @test occursin("…", split(show_vec, '\n')[3])
    @test occursin("weights: $([1,1,1])", show_h)
    @test occursin("closed: left", show_h)
    @test occursin("isdensity: false", show_h)
end


@testset "Histogram norm and normalize" begin
    rng = MersenneTwister(345678)
    edges = (
        cumsum(rand(rng) * rand(rng, 9)),
        cumsum(rand(rng, 1:10) * rand(rng, 1:100, 11)),
        cumsum(5 * rand(rng) * rand(rng, 14))
    )

    n = 100000

    data = (
        maximum(edges[1]) .* (randn(rng, n) ./ 6 .+ 0.5),
        rand(rng, 1:maximum(edges[2]), n),
        maximum(edges[3]) .* rand(rng, n)
    )

    h = fit(Histogram, data, edges, closed = :left)

    weight_sum = sum(h.weights)
    bin_vols = [ x * y * z for x in diff(edges[1]), y in diff(edges[2]), z in diff(edges[3])]

    @test norm(h) ≈ sum(h.weights .* bin_vols)

    @test @inferred(normalize(h, mode = :none)) == h


    h_pdf = normalize(h, mode = :pdf)
    @test h_pdf.weights ≈ h.weights ./ bin_vols ./ weight_sum
    @test h_pdf.isdensity == true
    @test @inferred(norm(h_pdf)) ≈ 1
#    @test @inferred(normalize(h_pdf, mode = :pdf)) == h_pdf
    @test @inferred(normalize(h_pdf, mode = :density)) == h_pdf
#    @test @inferred(normalize(h_pdf, mode = :probability)) == h_pdf

    h_density = normalize(h, mode = :density)
    @test h_density.weights ≈ h.weights ./ bin_vols
    @test h_density.isdensity == true
    @test @inferred(norm(h_density)) ≈ weight_sum
    @test @inferred(normalize(h_density, mode = :pdf)) ==
        Histogram(h_density.edges, h_density.weights .* (1/norm(h_density)), h_density.closed, true)
    @test normalize(h_density, mode = :pdf).weights ≈ h_pdf.weights
    @test normalize(h_density, mode = :density) == h_density
    @test normalize(h_density, mode = :probability).weights ≈ h_pdf.weights

    h_fraction = normalize(h, mode = :probability)
    @test sum(h_fraction.weights) ≈ 1
    @test h_fraction.isdensity == false
    @test normalize(h_fraction, mode = :pdf).weights ≈ h_pdf.weights
    @test normalize(h_fraction, mode = :density).weights ≈ h_pdf.weights
    @test normalize(h_fraction, mode = :probability).weights ≈ h_fraction.weights

    h_copy = deepcopy(float(h))
    @test @inferred(normalize!(h_copy, mode = :density)) == h_copy

    h2 = deepcopy(float(h))
    mod_h2 = normalize!(h2, mode = :density)
    @test mod_h2 === h2 && mod_h2.weights === h2.weights
    @test h2.weights == h_density.weights

    aux_weights = sqrt.(h.weights)
    divor0 = (a,b) -> (a == 0 && b == 0) ? 0 : a/b
    divor0_cmp = (a_n, a_d, b_n, b_d) -> maximum(abs.(map(divor0, a_n, a_d) - map(divor0, b_n, b_d))) < 1e-10

    h_pdf2, h_pdf2_aux = normalize(float(h), aux_weights, mode = :pdf)
    @test divor0_cmp(h_pdf2_aux, aux_weights, h_pdf2.weights, h.weights)

    h_density2, h_density2_aux = normalize(float(h), aux_weights, mode = :density)
    @test divor0_cmp(h_density2_aux, aux_weights, h_density2.weights, h.weights)

    h_density3, h_density3_aux = normalize(h_density2, h_density2_aux, mode = :pdf)
    @test divor0_cmp(h_density3_aux, h_density2_aux, h_density3.weights, h_density2.weights)
end


@testset "Histogram zero" begin
    h = fit(Histogram, (rand(100), rand(100)))
    h2 = @inferred zero(h)
    @test all(x -> x≈0, h2.weights)
    @test !(h.weights === h2.weights)
    @test h.edges == h2.edges
    @test h.closed == h2.closed
    @test h.isdensity == h2.isdensity
end


@testset "Histogram merge" begin
    histograms = [fit(Histogram, (rand(100), 10 * rand(100)), (0:0.1:1, 0:1:10)) for _ in 1:10]
    h = zero(histograms[1])
    merge!(h, histograms ...)
    @test h.weights == (+).((x->x.weights).(histograms)...)
    @test (@inferred merge(histograms...)) == h
end

@testset "midpoints" begin
    @test StatsBase.midpoints([1, 2, 4]) == [1.5, 3.0]
    @test StatsBase.midpoints(range(0, stop = 1, length = 5)) == 0.125:0.25:0.875
end

@testset "histogram with -0.0" begin
    @test fit(Histogram, [-0.0, 1.0]) == fit(Histogram, [0.0, 1.0])
    @test fit(Histogram, [-0.0, 1.0], closed=:right) ==
        fit(Histogram, [0.0, 1.0], closed=:right)
    @test fit(Histogram, [-0.0, -1.0]) == fit(Histogram, [0.0, -1.0])
    @test fit(Histogram, [-0.0, -1.0], closed=:right) ==
        fit(Histogram, [0.0, -1.0], closed=:right)

    @test fit(Histogram, [-0.0, 1.0], [-0.0, 0.5]) ==
        fit(Histogram, [0.0, 1.0], [0.0, 0.5]) ==
        fit(Histogram, [-0.0, 1.0], [0.0, 0.5]) ==
        fit(Histogram, [0.0, 1.0], [-0.0, 0.5]) ==
        fit(Histogram, [0.0, 1.0], 0.0:0.5:0.5) ==
        fit(Histogram, [-0.0, 1.0], 0.0:0.5:0.5)
    @test fit(Histogram, [-0.0, 1.0], [-0.5, -0.0]) ==
        fit(Histogram, [0.0, 1.0], [-0.5, -0.0]) ==
        fit(Histogram, [-0.0, 1.0], [-0.5, 0.0]) ==
        fit(Histogram, [0.0, 1.0], [-0.5, 0.0]) ==
        fit(Histogram, [-0.0, 1.0], -0.5:0.5:0.0) ==
        fit(Histogram, [0.0, 1.0], -0.5:0.5:0.0)
    @test fit(Histogram, [-0.0, 1.0], [-0.5, -0.0], closed=:right) ==
        fit(Histogram, [0.0, 1.0], [-0.5, 0.0], closed=:right) ==
        fit(Histogram, [0.0, 1.0], -0.5:0.5:0.0, closed=:right)
    @test fit(Histogram, [-0.0, 1.0], [-0.0, 0.5], closed=:right) ==
        fit(Histogram, [0.0, 1.0], [0.0, 0.5], closed=:right) ==
        fit(Histogram, [0.0, 1.0], [-0.0, 0.5], closed=:right) ==
        fit(Histogram, [-0.0, 1.0], [0.0, 0.5], closed=:right) ==
        fit(Histogram, [0.0, 1.0], 0.0:0.5:0.5, closed=:right) ==
        fit(Histogram, [-0.0, 1.0], 0.0:0.5:0.5, closed=:right)
    @test fit(Histogram, [-0.0, 1.0], [-0.5, -0.0], closed=:right) ==
        fit(Histogram, [0.0, 1.0], [-0.5, 0.0], closed=:right) ==
        fit(Histogram, [0.0, 1.0], [-0.5, -0.0], closed=:right) ==
        fit(Histogram, [-0.0, 1.0], [-0.5, 0.0], closed=:right) ==
        fit(Histogram, [0.0, 1.0], -0.5:0.5:0.0, closed=:right) ==
        fit(Histogram, [-0.0, 1.0], -0.5:0.5:0.0, closed=:right)

    # edges containing both -0.0 and 0.0, in either order, behave like two 0.0 edges:
    # the bin between them is empty and the neighbouring bins get the same observations
    for closed in (:left, :right)
        obs = [-0.5, -0.0, 0.0, 0.5]
        w = fit(Histogram, obs, [-1.0, 0.0, 0.0, 1.0], closed=closed).weights
        @test fit(Histogram, obs, [-1.0, -0.0, 0.0, 1.0], closed=closed).weights == w
        @test fit(Histogram, obs, [-1.0, 0.0, -0.0, 1.0], closed=closed).weights == w
        @test w == (closed == :left ? [1, 0, 3] : [3, 0, 1])
    end
    # ranges containing -0.0 behave like the corresponding ranges with 0.0
    @test fit(Histogram, [-0.5, 0.0], LinRange(-1.0, -0.0, 3)) ==
        fit(Histogram, [-0.5, 0.0], LinRange(-1.0, 0.0, 3))
    @test fit(Histogram, [-0.5, 0.0], LinRange(-1.0, -0.0, 3), closed=:right) ==
        fit(Histogram, [-0.5, -0.0], LinRange(-1.0, 0.0, 3), closed=:right)
    @test fit(Histogram, [0.0, 0.5], UnitRange(-0.0, 1.0)) ==
        fit(Histogram, [-0.0, 0.5], UnitRange(0.0, 1.0))
    @test fit(Histogram, [0.0, 0.5], UnitRange(-0.0, 1.0), closed=:right) ==
        fit(Histogram, [-0.0, 0.5], UnitRange(0.0, 1.0), closed=:right)
    # NaN is dropped
    @test fit(Histogram, [NaN, 0.5], 0.0:0.5:1.0).weights == [0, 1]
    @test fit(Histogram, [NaN, 0.5], 0.0:0.5:1.0, closed=:right).weights == [1, 0]
    @test fit(Histogram, [NaN, 0.5], [0.0, 0.5, 1.0]).weights == [0, 1]
    @test fit(Histogram, [NaN, 0.5], [0.0, 0.5, 1.0], closed=:right).weights == [1, 0]
end

end # @testset "StatsBase.Histogram"
