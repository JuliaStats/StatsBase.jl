using StatsBase
using Base.Test
import Base: maxabs
import StatsBase: norepeat, randi

srand(1234)

#### randi

n = 10^5

x = [randi(10) for i = 1:n]
@test isa(x, Vector{Int})
@test extrema(x) == (1, 10)
@test isapprox(proportions(x, 1:10), fill(0.1, 10), atol=5.0e-3)
@test randi(MersenneTwister(1), 1000) == randi(MersenneTwister(1), 1000)


x = [randi(3, 12) for i = 1:n]
@test isa(x, Vector{Int})
@test extrema(x) == (3, 12)
@test isapprox(proportions(x, 3:12), fill(0.1, 10), atol=5.0e-3)

#### sample with replacement

function check_sample_wrep(a::AbstractArray, vrgn, ptol::Real; ordered::Bool=false)
    vmin, vmax = vrgn
    (amin, amax) = extrema(a)
    @test vmin <= amin <= amax <= vmax
    n = vmax - vmin + 1
    p0 = fill(1/n, n)
    if ordered
        @test issorted(a)
        if ptol > 0
            @test isapprox(proportions(a, vmin:vmax), p0, atol=ptol)
        end
    else
        @test !issorted(a)
        ncols = size(a,2)
        if ncols == 1
            @test isapprox(proportions(a, vmin:vmax), p0, atol=ptol)
        else
            for j = 1:ncols
                aj = view(a, :, j)
                @test isapprox(proportions(aj, vmin:vmax), p0, atol=ptol)
            end
        end
    end
end

import StatsBase: direct_sample!

a = direct_sample!(1:10, zeros(Int, n, 3))
check_sample_wrep(a, (1, 10), 5.0e-3; ordered=false)

a = direct_sample!(3:12, zeros(Int, n, 3))
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=false)

a = direct_sample!([11:20;], zeros(Int, n, 3))
check_sample_wrep(a, (11, 20), 5.0e-3; ordered=false)

@test direct_sample!(MersenneTwister(1), 1:10, zeros(Int, 6)) ==
      direct_sample!(MersenneTwister(1), 1:10, zeros(Int, 6))

a = sample(3:12, n)
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=false)

a = sample(3:12, n; ordered=true)
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=true)

a = sample(3:12, 10; ordered=true)
check_sample_wrep(a, (3, 12), 0; ordered=true)

@test sample(MersenneTwister(1), 1:10, 10) == sample(MersenneTwister(1), 1:10, 10)

#### sampling pairs

srand(1);

@test samplepair(2)  === (1, 2)
@test samplepair(10) === (7, 3)

@test samplepair([3, 4, 2, 6, 8]) === (4, 3)
@test samplepair([1, 2])          === (1, 2)

@test samplepair(MersenneTwister(1), 1000) == samplepair(MersenneTwister(1), 1000)

#### sample without replacement

function check_sample_norep(a::AbstractArray, vrgn, ptol::Real; ordered::Bool=false)
    # each column of a for one run

    vmin, vmax = vrgn
    (amin, amax) = extrema(a)
    @test vmin <= amin <= amax <= vmax
    n = vmax - vmin + 1

    for j = 1:size(a,2)
        aj = view(a,:,j)
        @assert norepeat(aj)
        if ordered
            @assert issorted(aj)
        end
    end

    if ptol > 0
        p0 = fill(1/n, n)
        if ordered
            @test isapprox(proportions(a, vmin:vmax), p0, atol=ptol)
        else
            b = transpose(a)
            for j = 1:size(b,2)
                bj = view(b,:,j)
                @test isapprox(proportions(bj, vmin:vmax), p0, atol=ptol)
            end
        end
    end
end

import StatsBase: knuths_sample!, fisher_yates_sample!, self_avoid_sample!
import StatsBase: seqsample_a!, seqsample_c!

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    knuths_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)
@test knuths_sample!(MersenneTwister(1), 1:10, zeros(Int, 6)) == 
      knuths_sample!(MersenneTwister(1), 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    fisher_yates_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)
@test fisher_yates_sample!(MersenneTwister(1), 1:10, zeros(Int, 6)) == 
      fisher_yates_sample!(MersenneTwister(1), 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    self_avoid_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)
@test self_avoid_sample!(MersenneTwister(1), 1:10, zeros(Int, 6)) == 
      self_avoid_sample!(MersenneTwister(1), 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    seqsample_a!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=true)
@test seqsample_a!(MersenneTwister(1), 1:10, zeros(Int, 6)) == 
      seqsample_a!(MersenneTwister(1), 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    seqsample_c!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=true)
@test seqsample_c!(MersenneTwister(1), 1:10, zeros(Int, 6)) == 
      seqsample_c!(MersenneTwister(1), 1:10, zeros(Int, 6))

a = sample(3:12, 5; replace=false)
check_sample_norep(a, (3, 12), 0; ordered=false)

a = sample(3:12, 5; replace=false, ordered=true)
check_sample_norep(a, (3, 12), 0; ordered=true)

# test of weighted sampling without replacement
a = [1:10;]
wv = WeightVec([zeros(6); 1:4])
x = vcat([sample(a, wv, 1, replace=false) for j in 1:100000]...)
@test minimum(x) == 7
@test maximum(x) == 10
@test maximum(abs(proportions(x) - (1:4)/10)) < 0.01

x = vcat([sample(a, wv, 2, replace=false) for j in 1:50000]...)
exact2 = [0.117261905, 0.220634921, 0.304166667, 0.357936508]
@test minimum(x) == 7
@test maximum(x) == 10
@test maximum(abs(proportions(x) - exact2)) < 0.01

x = vcat([sample(a, wv, 4, replace=false) for j in 1:10000]...)
@test minimum(x) == 7
@test maximum(x) == 10
@test maximum(abs(proportions(x) - 0.25)) == 0

@test_throws DimensionMismatch sample(a, wv, 5, replace=false)

wv = WeightVec([zeros(5); 1:4; -1])
@test_throws ErrorException sample(a, wv, 1, replace=false)

