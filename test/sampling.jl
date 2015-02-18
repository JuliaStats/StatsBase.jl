# Test sample functions

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
@test_approx_eq_eps proportions(x, 1:10) fill(0.1, 10) 5.0e-3


x = [randi(3, 12) for i = 1:n]
@test isa(x, Vector{Int})
@test extrema(x) == (3, 12)
@test_approx_eq_eps proportions(x, 3:12) fill(0.1, 10) 5.0e-3

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
            @test_approx_eq_eps proportions(a, vmin:vmax) p0 ptol
        end
    else
        @test !issorted(a)
        ncols = size(a,2)
        if ncols == 1
            @test_approx_eq_eps proportions(a, vmin:vmax) p0 ptol
        else
            for j = 1:ncols
                aj = view(a, :, j)
                @test_approx_eq_eps proportions(aj, vmin:vmax) p0 ptol
            end
        end
    end
end

import StatsBase: direct_sample!, xmultinom_sample!

a = direct_sample!(1:10, zeros(Int, n, 3))
check_sample_wrep(a, (1, 10), 5.0e-3; ordered=false)

a = direct_sample!(3:12, zeros(Int, n, 3))
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=false)

a = direct_sample!([11:20;], zeros(Int, n, 3))
check_sample_wrep(a, (11, 20), 5.0e-3; ordered=false)

a = xmultinom_sample!(3:12, zeros(Int, n))
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=true)

a = xmultinom_sample!(101:200, zeros(Int, 10))
check_sample_wrep(a, (101, 200), 0; ordered=true)

a = sample(3:12, n)
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=false)

a = sample(3:12, n; ordered=true)
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=true)

a = sample(3:12, 10; ordered=true)
check_sample_wrep(a, (3, 12), 0; ordered=true)


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
            @test_approx_eq_eps proportions(a, vmin:vmax) p0 ptol
        else
            b = transpose(a)
            for j = 1:size(b,2)
                bj = view(b,:,j)
                @test_approx_eq_eps proportions(bj, vmin:vmax) p0 ptol
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

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    fisher_yates_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    self_avoid_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    seqsample_a!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=true)

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    seqsample_c!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=true)

a = sample(3:12, 5; replace=false)
check_sample_norep(a, (3, 12), 0; ordered=false)

a = sample(3:12, 5; replace=false, ordered=true)
check_sample_norep(a, (3, 12), 0; ordered=true)



