# Test sample functions

using StatsBase
using Base.Test
import Base: maxabs

srand(1234)

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

import StatsBase: direct_sample!, xmultinom_sample!

a = direct_sample!(1:10, zeros(Int, n, 3))
for j = 1:size(a,2)
    aj = a[:,j]
    @test extrema(aj) == (1, 10)
    @test_approx_eq_eps proportions(aj, 1:10) fill(0.1, 10) 5.0e-3
end

a = direct_sample!(3:12, zeros(Int, n, 3))
for j = 1:size(a,2)
    aj = a[:,j]
    @test extrema(aj) == (3, 12)
    @test_approx_eq_eps proportions(aj, 3:12) fill(0.1, 10) 5.0e-3
end

a = direct_sample!([11:20], zeros(Int, n, 3))
for j = 1:size(a,2)
    aj = a[:,j]
    @test extrema(aj) == (11, 20)
    @test_approx_eq_eps proportions(aj, 11:20) fill(0.1, 10) 5.0e-3
end

a = xmultinom_sample!(3:12, zeros(Int, n))
@test issorted(a)
@test a[1] == 3 && a[end] == 12
@test_approx_eq_eps proportions(a, 3:12) fill(0.1, 10) 5.0e-3

a = xmultinom_sample!(101:200, zeros(Int, 10))
@test issorted(a)
@test 101 <= minimum(a) <= maximum(a) <= 200

a = sample(3:12, n)
@test !issorted(a)
@test extrema(a) == (3, 12)
@test_approx_eq_eps proportions(a, 3:12) fill(0.1, 10) 5.0e-3

a = sample(3:12, n; ordered=true)
@test issorted(a)
@test extrema(a) == (3, 12)
@test_approx_eq_eps proportions(a, 3:12) fill(0.1, 10) 5.0e-3

a = sample(3:12, 10; ordered=true)
@test issorted(a)
@test 3 <= minimum(a) <= maximum(a) <= 12






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

