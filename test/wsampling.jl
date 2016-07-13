using StatsBase
using Base.Test
import Base: maxabs
import StatsBase: norepeat

srand(1234)

#### weighted sample with replacement

function check_wsample_wrep(a::AbstractArray, vrgn, wv::WeightVec, ptol::Real; ordered::Bool=false)
    K = length(wv)
    (vmin, vmax) = vrgn
    (amin, amax) = extrema(a)
    @test vmin <= amin <= amax <= vmax
    p0 = values(wv) ./ sum(wv)
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

import StatsBase: direct_sample!, alias_sample!

n = 10^5
wv = weights([0.2, 0.8, 0.4, 0.6])

a = direct_sample!(4:7, wv, zeros(Int, n, 3))
check_wsample_wrep(a, (4, 7), wv, 5.0e-3; ordered=false)

a = alias_sample!(4:7, wv, zeros(Int, n, 3))
check_wsample_wrep(a, (4, 7), wv, 5.0e-3; ordered=false)

a = sample(4:7, wv, n; ordered=false)
check_wsample_wrep(a, (4, 7), wv, 5.0e-3; ordered=false)

a = sample(4:7, wv, n; ordered=true)
check_wsample_wrep(a, (4, 7), wv, 5.0e-3; ordered=true)


#### weighted sampling without replacement

function check_wsample_norep(a::AbstractArray, vrgn, wv::WeightVec, ptol::Real; ordered::Bool=false)
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
        p0 = values(wv) ./ sum(wv)
        @test_approx_eq_eps proportions(a[1,:], vmin:vmax) p0 ptol
    end
end

import StatsBase: naive_wsample_norep!

n = 10^5
wv = weights([0.2, 0.8, 0.4, 0.6])

a = zeros(Int, 3, n)
for j = 1:n
    naive_wsample_norep!(4:7, wv, view(a,:,j))
end
check_wsample_norep(a, (4, 7), wv, 5.0e-3; ordered=false)
