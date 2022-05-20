using StatsBase
using Test, Random, StableRNGs

Random.seed!(1234)

n = 100000

# test that rng specification is working correctly
# a) if the same rng is passed to a sample function twice,
#    the results should be the same (repeatability)
# b) not specifying a rng should be the same as specifying Random.GLOBAL_RNG
function test_rng_use(func, non_rng_args...)
    # some sampling methods mutate a passed array and return it
    # so that the tests don't pass trivially, we need to copy those
    # pre-allocated storage arrays

    # repeatability
    @test func(MersenneTwister(1), deepcopy(non_rng_args)...) ==
          func(MersenneTwister(1), deepcopy(non_rng_args)...)
    # default RNG is Random.GLOBAL_RNG
    Random.seed!(47)
    x = func(deepcopy(non_rng_args)...)
    Random.seed!(47)
    y = func(Random.GLOBAL_RNG, deepcopy(non_rng_args)...)
    @test x == y
end

#### sample with replacement

function check_sample_wrep(a::AbstractArray, vrgn, ptol::Real; ordered::Bool=false, rev::Bool=false)
    vmin, vmax = vrgn
    (amin, amax) = extrema(a)
    @test vmin <= amin <= amax <= vmax
    n = vmax - vmin + 1
    p0 = fill(1/n, n)
    if ordered
        @test issorted(a; rev=rev)
        if ptol > 0
            @test isapprox(proportions(a, vmin:vmax), p0, atol=ptol)
        end
    else
        @test !issorted(a; rev=rev)
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

test_rng_use(direct_sample!, 1:10, zeros(Int, 6))

a = sample(3:12, n)
check_sample_wrep(a, (3, 12), 5.0e-3; ordered=false)

for rev in (true, false), T in (Int, Int16, Float64, Float16, BigInt, ComplexF64, Rational{Int})
    r = rev ? reverse(3:12) : (3:12)
    r = T===Int ? r : T.(r)
    aa = Int.(sample(r, n; ordered=true))
    check_sample_wrep(aa, (3, 12), 5.0e-3; ordered=true, rev=rev)

    aa = Int.(sample(r, 10; ordered=true))
    check_sample_wrep(aa, (3, 12), 0; ordered=true, rev=rev)
end

@test StatsBase._storeindices(1, 1, BigFloat) == StatsBase._storeindices(1, 1, BigFloat) == false

test_rng_use(sample, 1:10, 10)

@testset "sampling pairs" begin

    rng = StableRNG(1)

    @test samplepair(rng, 2)  ===  (2, 1)
    @test samplepair(rng, 10) === (5, 6)

    @test samplepair(rng, [3, 4, 2, 6, 8]) === (3, 8)
    @test samplepair(rng, [1, 2])          === (1, 2)
end

test_rng_use(samplepair, 1000)

#### sample without replacement

function check_sample_norep(a::AbstractArray, vrgn, ptol::Real; ordered::Bool=false, rev::Bool=false)
    # each column of a for one run

    vmin, vmax = vrgn
    (amin, amax) = extrema(a)
    @test vmin <= amin <= amax <= vmax
    n = vmax - vmin + 1

    for j = 1:size(a,2)
        aj = view(a,:,j)
        @assert allunique(aj)
        if ordered
            @assert issorted(aj, rev=rev)
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
import StatsBase: seqsample_a!, seqsample_c!, seqsample_d!

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    knuths_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)

test_rng_use(knuths_sample!, 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    fisher_yates_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)

test_rng_use(fisher_yates_sample!, 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    self_avoid_sample!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=false)

test_rng_use(self_avoid_sample!, 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    seqsample_a!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=true)

test_rng_use(seqsample_a!, 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    seqsample_c!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=true)

test_rng_use(seqsample_c!, 1:10, zeros(Int, 6))

a = zeros(Int, 5, n)
for j = 1:size(a,2)
    seqsample_d!(3:12, view(a,:,j))
end
check_sample_norep(a, (3, 12), 5.0e-3; ordered=true)

test_rng_use(seqsample_d!, 1:10, zeros(Int, 6))

a = sample(3:12, 5; replace=false)
check_sample_norep(a, (3, 12), 0; ordered=false)

a = sample(3:12, 5; replace=false, ordered=true)
check_sample_norep(a, (3, 12), 0; ordered=true)

a = sample(reverse(3:12), 5; replace=false, ordered=true)
check_sample_norep(a, (3, 12), 0; ordered=true, rev=true)

# tests of multidimensional sampling

a = sample(3:12, (2, 2); replace=false)
check_sample_norep(a, (3, 12), 0; ordered=false)

@test sample(1:1, (2, 2); replace=true) == ones(Int, 2, 2)

# test of weighted sampling without replacement
a = [1:10;]
wv = Weights([zeros(6); 1:4])
x = vcat([sample(a, wv, 1, replace=false) for j in 1:100000]...)
@test minimum(x) == 7
@test maximum(x) == 10
@test maximum(abs, proportions(x) .- (1:4)/10) < 0.01

x = vcat([sample(a, wv, 2, replace=false) for j in 1:50000]...)
exact2 = [0.117261905, 0.220634921, 0.304166667, 0.357936508]
@test minimum(x) == 7
@test maximum(x) == 10
@test maximum(abs, proportions(x) .- exact2) < 0.01

x = vcat([sample(a, wv, 4, replace=false) for j in 1:10000]...)
@test minimum(x) == 7
@test maximum(x) == 10
@test maximum(abs, proportions(x) .- 0.25) == 0

@test_throws DimensionMismatch sample(a, wv, 5, replace=false)

wv = Weights([zeros(5); 1:4; -1])
@test_throws ErrorException sample(a, wv, 1, replace=false)

#### weighted sampling with dimension

# weights respected; this works because of the 0-weight
@test sample([1, 2], Weights([0, 1]), (2,2)) == [2 2 ; 2 2]
wm =  sample(collect(1:4), Weights(1:4), (2,2), replace=false)
@test size(wm) == (2, 2) # correct shape
@test length(Set(wm)) == 4 # no duplicates in elements


#### check that sample and sample! do the same thing
function test_same(;kws...)
    wv = Weights(rand(20))
    Random.seed!(1)
    x1 = sample(1:20, wv, 10; kws...)
    Random.seed!(1)
    x2 = zeros(Int, 10)
    sample!(1:20, wv, x2; kws...)
    @test x1 == x2
end

test_same()
test_same(replace=true)
test_same(replace=false)
test_same(replace=true, ordered=true)
test_same(replace=false, ordered=true)
test_same(replace=true, ordered=false)
test_same(replace=false, ordered=false)

# Test that sample! throws when output shares memory with inputs
x = rand(10)
y = rand(10)
@test_throws ArgumentError sample!(x, x)
@test_throws ArgumentError sample!(x, weights(y), x)
@test_throws ArgumentError sample!(x, weights(x), y)
@test_throws ArgumentError sample!(x, weights(x), x)
@test_throws ArgumentError sample!(view(x, 2:4), view(x, 3:5))
@test_throws ArgumentError sample!(view(x, 2:4), weights(view(x, 3:5)), y)
@test_throws ArgumentError sample!(view(x, 2:4), weights(view(x, 3:5)), view(x, 1:2))
# These corner cases should theoretically succeed
# but the second currently fails as Base.mightalias is not smart enough
sample!(view(x, 2:4), view(x, 5:6))
@test_broken sample!(view(x, 2:4), weights(view(x, 5:6)), y)
