# See src/hist.jl for meaning of "FIXME: closed" comments.

using StatsBase
using Base.Test

@testset "StatsBase.Histogram" begin


@testset "Histogram binindex and binvolume" begin
    edg1 = -2:0.5:9
    edg1f0 = -2:0.5f0:9
    edg2 = [-2, -1, 2, 7, 19]
    h1 = Histogram(edg1, :left)
    h2 = Histogram((edg1, edg2), :left)
    h3 = Histogram((edg1f0, edg2), :left)
    @test @inferred StatsBase.binindex(h1, -0.5) == 4
    @test @inferred StatsBase.binindex(h2, (1.5, 2)) == (8, 3)

    @test [StatsBase.binvolume(h1, i) for i in indices(h1.weights, 1)] ≈ diff(edg1)
    @test [StatsBase.binvolume(h2, (i,j)) for i in indices(h2.weights, 1), j in indices(h2.weights, 2)] ≈ diff(edg1) * diff(edg2)'

    @test typeof(@inferred(StatsBase.binvolume(h2, (1,1)))) == Float64
    @test typeof(@inferred(StatsBase.binvolume(h3, (1,1)))) == Float32
    @test typeof(@inferred(StatsBase.binvolume(Float64, h3, (1,1)))) == Float64
end


@testset "Histogram append" begin
    # FIXME: closed (all lines in this block):
    h = Histogram(0:20:100, Float64, :left, false)
    @test @inferred(append!(h, 0:0.5:99.99)) == h
    @test append!(Histogram(0:20:100, Float64, :left, false), 0:0.5:99.99).weights ≈ [40,40,40,40,40]
    @test append!(Histogram(0:20:100, Float64, :left, true), 0:0.5:99.99).weights ≈ [2,2,2,2,2]
    @test append!(Histogram(0:20:100, Float64, :left, false), 0:0.5:99.99, fill(2, 200)).weights ≈ [80,80,80,80,80]
    @test append!(Histogram(0:20:100, Float64, :left, true), 0:0.5:99.99, fill(2, 200)).weights ≈ [4,4,4,4,4]
end


@testset "Histogram fit" begin
    # FIXME: closed (all lines in this block):
    @test sum(fit(Histogram,[1,2,3], closed=:left).weights) == 3  # FIXME: closed
    @test fit(Histogram,Int[], closed=:left).weights == Int[]  # FIXME: closed
    @test fit(Histogram,[1], closed=:left).weights == [1]  # FIXME: closed
    @test fit(Histogram,[1,2,3],[0,2,4], closed=:left) == Histogram([0,2,4],[1,2], :left)  # FIXME: closed
    @test fit(Histogram,[1,2,3],[0,2,4], closed=:left) != Histogram([0,2,4],[1,1], :left)  # FIXME: closed
    @test fit(Histogram,[1,2,3],0:2:4, closed=:left) == Histogram(0:2:4,[1,2], :left)  # FIXME: closed
    @test all(fit(Histogram,[0:99;]/100,0.0:0.01:1.0, closed=:left).weights .==1)  # FIXME: closed
    @test fit(Histogram,[1,1,1,1,1], closed=:left).weights[1] == 5  # FIXME: closed
    @test sum(fit(Histogram,(rand(100),rand(100)), closed=:left).weights) == 100  # FIXME: closed
    @test fit(Histogram,1:100,nbins=5,closed=:right).weights == [20,20,20,20,20]
    @test fit(Histogram,1:100,nbins=5,closed=:left).weights == [19,20,20,20,20,1]
    @test fit(Histogram,0:99,nbins=5,closed=:right).weights == [1,20,20,20,20,19]
    @test fit(Histogram,0:99,nbins=5,closed=:left).weights == [20,20,20,20,20]

    # FIXME: closed (all lines in this block):
    @test fit(Histogram,(0:99,0:99),nbins=5, closed=:left).weights == diagm([20,20,20,20,20])
    @test fit(Histogram,(0:99,0:99),nbins=(5,5), closed=:left).weights == diagm([20,20,20,20,20])

    # FIXME: closed (all lines in this block):
    @test fit(Histogram,0:99,weights(ones(100)),nbins=5, closed=:left).weights == [20,20,20,20,20]
    @test fit(Histogram,0:99,weights(2*ones(100)),nbins=5, closed=:left).weights == [40,40,40,40,40]
    @test fit(Histogram{Int32},0:99,weights(2*ones(100)),nbins=5, closed=:left).weights == [40,40,40,40,40]
    @test fit(Histogram{Float32},0:99,weights(2*ones(100)),nbins=5, closed=:left).weights == [40,40,40,40,40]
end


@testset "Histogram element type" begin
    # FIXME: closed (all lines in this block):
    @test eltype(@inferred(fit(Histogram,1:100,weights(ones(Int,100)),nbins=5, closed=:left)).weights) == Int
    @test eltype(@inferred(fit(Histogram{Float32},1:100,weights(ones(Int,100)),nbins=5, closed=:left)).weights) == Float32
    @test eltype(@inferred(fit(Histogram,1:100,weights(ones(Float64,100)),nbins=5, closed=:left)).weights) == Float64
    @test eltype(@inferred(fit(Histogram{Float32},1:100,weights(ones(Float64,100)),nbins=5, closed=:left)).weights) == Float32
end


@testset "histrange" begin
    # Note: atm histrange must be qualified
    @test @inferred(StatsBase.histrange(Float64[], 0, :left)) == 0.0:1.0:0.0
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
    @test StatsBase.histrange(Int64[1:10;], 1, :left) == 0:10:20

    # FIXME: closed (all lines in this block):
    @test StatsBase.histrange([0, 1, 2, 3], 4, :left) == 0.0:1.0:4.0
    @test StatsBase.histrange([0, 1, 1, 3], 4, :left) == 0.0:1.0:4.0
    @test StatsBase.histrange([0, 9], 4, :left) == 0.0:5.0:10.0
    @test StatsBase.histrange([0, 19], 4, :left) == 0.0:5.0:20.0
    @test StatsBase.histrange([0, 599], 4, :left) == 0.0:200.0:600.0
    @test StatsBase.histrange([-1, -1000], 4, :left) == -1000.0:500.0:0.0

    # Base issue #13326
    l,h = extrema(StatsBase.histrange([typemin(Int),typemax(Int)], 10, :left))  # FIXME: closed
    @test l <= typemin(Int)
    @test h >= typemax(Int)

    # FIXME: closed (all lines in this block):
    @test_throws ArgumentError StatsBase.histrange([1, 10], 0, :left)
    @test_throws ArgumentError StatsBase.histrange([1, 10], -1, :left)
    @test_throws ArgumentError StatsBase.histrange([1.0, 10.0], 0, :left)
    @test_throws ArgumentError StatsBase.histrange([1.0, 10.0], -1, :left)
    @test_throws ArgumentError StatsBase.histrange(Float64[],-1, :left)
    @test_throws ArgumentError StatsBase.histrange([0.], 0, :left)
end


@testset "Histogram show" begin
    # hist show
    show_h = sprint(show, fit(Histogram,[0,1,2], closed=:left))  # FIXME: closed
    @test contains(show_h, "edges:\n  0.0:1.0:3.0")
    @test contains(show_h, "weights: $([1,1,1])")
    @test contains(show_h, "closed: left")
    @test contains(show_h, "isdensity: false")
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
        maximum(edges[1]) * (randn(rng, n) / 6 + 0.5),
        rand(rng, 1:maximum(edges[2]), n),
        maximum(edges[3]) * rand(rng, n)
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
    @test @inferred(normalize(h_pdf, mode = :pdf)) == h_pdf
    @test @inferred(normalize(h_pdf, mode = :density)) == h_pdf

    h_density = normalize(h, mode = :density)
    @test h_density.weights ≈ h.weights ./ bin_vols
    @test h_density.isdensity == true
    @test @inferred(norm(h_density)) ≈ weight_sum
    @test @inferred(normalize(h_density, mode = :pdf)) ==
        Histogram(h_density.edges, h_density.weights .* (1/norm(h_density)), h_density.closed, true)
    @test normalize(h_density, mode = :pdf).weights ≈ h_pdf.weights
    @test normalize(h_density, mode = :density) == h_density

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


end # @testset "StatsBase.Histogram"
