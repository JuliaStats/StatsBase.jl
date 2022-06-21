using StatsBase
using Random, Test

Random.seed!(1234)

#### weighted sample with replacement

function check_wsample_wrep(a::AbstractArray, vrgn, wv::AbstractWeights, ptol::Real;
                            ordered::Bool=false, rev::Bool=false)
    K = length(wv)
    (vmin, vmax) = vrgn
    (amin, amax) = extrema(a)
    @test vmin <= amin <= amax <= vmax
    p0 = wv ./ sum(wv)
    rev && reverse!(p0)
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

import StatsBase: direct_sample!, alias_sample!

n = 10^6
wv = weights([0.2, 0.8, 0.4, 0.6])

for wv in (
    weights([0.2, 0.8, 0.4, 0.6]),
    weights([2, 8, 4, 6]),
    weights(Float32[0.2, 0.8, 0.4, 0.6]),
    Weights(Float32[0.2, 0.8, 0.4, 0.6], 2),
    Weights([2, 8, 4, 6], 20.0),
)
    a = direct_sample!(4:7, wv, zeros(Int, n, 3))
    check_wsample_wrep(a, (4, 7), wv, 5.0e-3; ordered=false)
    test_rng_use(direct_sample!, 4:7, wv, zeros(Int, 100))

    a = alias_sample!(4:7, wv, zeros(Int, n, 3))
    check_wsample_wrep(a, (4, 7), wv, 5.0e-3; ordered=false)

    a = sample(4:7, wv, n; ordered=false)
    check_wsample_wrep(a, (4, 7), wv, 5.0e-3; ordered=false)
end

for rev in (true, false), T in (Int, Int16, Float64, Float16, BigInt, ComplexF64, Rational{Int})
    r = rev ? reverse(4:7) : (4:7)
    r = T===Int ? r : T.(r)
    aa = Int.(sample(r, wv, n; ordered=true))
    check_wsample_wrep(aa, (4, 7), wv, 5.0e-3; ordered=true, rev=rev)
    aa = Int.(sample(r, wv, 10; ordered=true))
    check_wsample_wrep(aa, (4, 7), wv, -1; ordered=true, rev=rev)
end

#### weighted sampling without replacement

function check_wsample_norep(a::AbstractArray, vrgn, wv::AbstractWeights, ptol::Real;
                             ordered::Bool=false, rev::Bool=false)
    # each column of a for one run

    vmin, vmax = vrgn
    (amin, amax) = extrema(a)
    @test vmin <= amin <= amax <= vmax
    n = vmax - vmin + 1

    for j = 1:size(a,2)
        aj = view(a,:,j)
        @assert allunique(aj)
        if ordered
            @assert issorted(aj; rev=rev)
        end
    end

    if ptol > 0
        p0 = wv ./ sum(wv)
        rev && reverse!(p0)
        @test isapprox(proportions(a[1,:], vmin:vmax), p0, atol=ptol)
    end
end

import StatsBase: naive_wsample_norep!, efraimidis_a_wsample_norep!,
                  efraimidis_ares_wsample_norep!, efraimidis_aexpj_wsample_norep!

n = 10^5
wv = weights([0.2, 0.8, 0.4, 0.6])

a = zeros(Int, 3, n)
for j = 1:n
    naive_wsample_norep!(4:7, wv, view(a,:,j))
end
check_wsample_norep(a, (4, 7), wv, 5.0e-3; ordered=false)
test_rng_use(naive_wsample_norep!, 4:7, wv, zeros(Int, 2))

a = zeros(Int, 3, n)
for j = 1:n
    efraimidis_a_wsample_norep!(4:7, wv, view(a,:,j))
end
check_wsample_norep(a, (4, 7), wv, 5.0e-3; ordered=false)
test_rng_use(efraimidis_a_wsample_norep!, 4:7, wv, zeros(Int, 2))

a = zeros(Int, 3, n)
for j = 1:n
    efraimidis_ares_wsample_norep!(4:7, wv, view(a,:,j))
end
check_wsample_norep(a, (4, 7), wv, 5.0e-3; ordered=false)
test_rng_use(efraimidis_ares_wsample_norep!, 4:7, wv, zeros(Int, 2))

a = zeros(Int, 3, n)
for j = 1:n
    efraimidis_aexpj_wsample_norep!(4:7, wv, view(a,:,j))
end
check_wsample_norep(a, (4, 7), wv, 5.0e-3; ordered=false)
test_rng_use(efraimidis_aexpj_wsample_norep!, 4:7, wv, zeros(Int, 2))

a = sample(4:7, wv, 3; replace=false, ordered=false)
check_wsample_norep(a, (4, 7), wv, -1; ordered=false)

for rev in (true, false), T in (Int, Int16, Float64, Float16, BigInt, ComplexF64, Rational{Int})
    r = rev ? reverse(4:7) : (4:7)
    r = T===Int ? r : T.(r)
    aa = Int.(sample(r, wv, 3; replace=false, ordered=true))
    check_wsample_norep(aa, (4, 7), wv, -1; ordered=true, rev=rev)
end
