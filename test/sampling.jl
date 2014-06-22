# Test sample functions

using StatsBase
using Base.Test
import Base: maxabs

srand(1234)

#### auxiliary 

function norepeat(a::AbstractArray)
    sa = sort(a)
    for i = 2:length(a)
        if a[i] == a[i-1]
            return false
        end
    end
    return true
end


#### randi

const randi = StatsBase.randi
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

a = direct_sample!([11:20], zeros(Int, n, 3))
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



#### sample with replacement

# ## individual

# x = [sample(1:5) for _ = 1:10^5]
# p = fill(0.2, 5)
# @test isa(x, Vector{Int})
# @test_approx_eq_eps est_p(x, 5) p 0.02

# ## unordered

# x = sample(1:5, 10^5)
# @test isa(x, Vector{Int})
# @test length(x) == 10^5
# @test !issorted(x)  # extremely unlikely to be sorted for n=10^5
# @test_approx_eq_eps est_p(x, 5) p 0.02
# @test_approx_eq_eps est_p(x[1:50000], 5) p 0.03

# ## ordered

# x = sample(1:5, 10^5; ordered=true)
# @test isa(x, Vector{Int})
# @test length(x) == 10^5
# @test issorted(x)
# @test_approx_eq_eps est_p(x, 5) p 0.02


# #### sample without replacement

# ## unordered (using samplepair)

# x = Array(Int, 2, n)
# for i = 1:n
#   x[:,i] = sample(1:5, 2; replace=false)
# end

# @test all(x[1,:] .!= x[2,:])
# @test_approx_eq_eps est_p(x[1,:], 5) p 0.02
# @test_approx_eq_eps est_p(x[2,:], 5) p 0.02

# ## unordered (using Fisher-Yates)

# x = Array(Int, 3, n)
# for i = 1:n
#   x[:,i] = sample(1:5, 3; replace=false)
# end

# @test all(x[1,:] .!= x[2,:])
# @test all(x[1,:] .!= x[3,:])
# @test all(x[2,:] .!= x[3,:])
# @test_approx_eq_eps est_p(x[1,:], 5) p 0.02
# @test_approx_eq_eps est_p(x[2,:], 5) p 0.02
# @test_approx_eq_eps est_p(x[3,:], 5) p 0.02

# ## unordered (using Self-avoid method)

# x = Array(Int, 3, 10000)
# for i = 1:10000
#   x[:,i] = sample(1:1000, 3; replace=false)
# end

# @test all(x[1,:] .!= x[2,:])
# @test all(x[1,:] .!= x[3,:])
# @test all(x[2,:] .!= x[3,:])

# ## ordered 

# x = Array(Int, 3, n)
# for i = 1:n
#   x[:,i] = sample(1:5, 3; replace=false, ordered=true)
# end

# @test all(x[1,:] .< x[2,:])
# @test all(x[2,:] .< x[3,:])
# @test_approx_eq_eps est_p(vec(x), 5) p 0.02


# #### weighted sample with replacement

# wv = weights([1.0, 2.0, 4.5, 2.5])
# p = [0.10, 0.20, 0.45, 0.25]

# ## individual

# x = Int[sample(wv) for _ = 1:10^5]
# @test isa(x, Vector{Int})
# @test_approx_eq_eps est_p(x, 4) p 0.02

# ## unordered

# x = sample(1:5, wv, 10^5)
# @test isa(x, Vector{Int})
# @test length(x) == 10^5
# @test !issorted(x)
# @test_approx_eq_eps est_p(x, 4) p 0.02

# ## ordered

# x = sample(1:5, wv, 10^5; ordered=true)
# @test isa(x, Vector{Int})
# @test length(x) == 10^5
# @test issorted(x)
# @test_approx_eq_eps est_p(x, 4) p 0.02

