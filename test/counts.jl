using StatsBase
using Test

n = 5000

# 1D integer counts

x = rand(1:5, n)
w = weights(rand(n))

c = counts(x, 5)
@test size(c) == (5,)
c0 = Int[count(v->v == i, x) for i in 1:5]
@test c == c0
@test counts(x .+ 1, 2:6) == c0
@test proportions(x, 1:5) ≈ (c0 ./ n)

c = counts(x)
@test size(c) == (5,)
c0 = Int[count(v->v == i, x) for i in 1:5]
@test c == c0
@test counts(x .+ 1, 2:6) == c0
@test proportions(x) ≈ (c0 ./ n)

c = counts(x, 5, w)
@test size(c) == (5,)
c0 = Float64[sum(w.values[x .== i]) for i in 1:5]
@test c                      ≈ c0
@test counts(x .+ 1, 2:6, w) ≈ c0
@test proportions(x, 1:5, w) ≈ (c0 ./ sum(w))

c = counts(x, w)
@test size(c) == (5,)
c0 = Float64[sum(w.values[x .== i]) for i in 1:5]
@test c                      ≈ c0
@test counts(x .+ 1, 2:6, w) ≈ c0
@test proportions(x, w)      ≈ (c0 ./ sum(w))

# 2D integer counts

x = rand(1:4, n)
y = rand(1:5, n)
w = weights(rand(n))

c = counts(x, y, (4, 5))
@test size(c) == (4, 5)
c0 = Int[count(t->t != 0,  (x .== i) .& (y .== j)) for i in 1:4, j in 1:5]
@test c == c0
@test counts(x .+ 2, y .+ 3, (3:6, 4:8)) == c0
@test proportions(x, y, (1:4, 1:5)) ≈ (c0 ./ n)

c = counts(x, y)
@test size(c) == (4, 5)
c0 = Int[count(t->t != 0, (x .== i) .& (y .== j)) for i in 1:4, j in 1:5]
@test c == c0
@test counts(x .+ 2, y .+ 3, (3:6, 4:8)) == c0
@test proportions(x, y,) ≈ (c0 ./ n)

c = counts(x, y, (4, 5), w)
@test size(c) == (4, 5)
c0 = Float64[sum(w.values[(x .== i) .& (y .== j)]) for i in 1:4, j in 1:5]
@test c                                     ≈ c0
@test counts(x .+ 2, y .+ 3, (3:6, 4:8), w) ≈ c0
@test proportions(x, y, (1:4, 1:5), w)      ≈ (c0 ./ sum(w))

c = counts(x, y, w)
@test size(c) == (4, 5)
c0 = Float64[sum(w.values[(x .== i) .& (y .== j)]) for i in 1:4, j in 1:5]
@test c                                     ≈ c0
@test counts(x .+ 2, y .+ 3, (3:6, 4:8), w) ≈ c0
@test proportions(x, y, w)                  ≈ (c0 ./ sum(w))


# count map

x = ["a", "b", "a", "a", "b", "c"]
w = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

cm = countmap(x)
@test cm["a"] == 3
@test cm["b"] == 2
@test cm["c"] == 1
pm = proportionmap(x)
@test pm["a"] ≈ (1/2)
@test pm["b"] ≈ (1/3)
@test pm["c"] ≈ (1/6)


# testing the radixsort branch of countmap
xx = repeat([6, 1, 3, 1], outer=100_000)
cm = countmap(xx)
@test cm == Dict(1 => 200_000, 3 => 100_000, 6 => 100_000)

# testing the radixsort-based addcounts
xx = repeat([6, 1, 3, 1], outer=100_000)
cm = Dict{Int, Int}()
StatsBase.addcounts_radixsort!(cm,xx)
@test cm == Dict(1 => 200_000, 3 => 100_000, 6 => 100_000)
xx2 = repeat([7, 1, 3, 1], outer=100_000)
StatsBase.addcounts_radixsort!(cm,xx2)
@test cm == Dict(1 => 400_000, 3 => 200_000, 6 => 100_000, 7 => 100_000)

# testing the Dict-based addcounts
cm = Dict{Int, Int}()
StatsBase.addcounts_dict!(cm,xx)
@test cm == Dict(1 => 200_000, 3 => 100_000, 6 => 100_000)

cm = countmap(x, weights(w))
@test cm["a"] == 5.5
@test cm["b"] == 4.5
@test cm["c"] == 3.5

@test cm == countmap(x, w)

pm = proportionmap(x, weights(w))
@test pm["a"] ≈ (5.5 / 13.5)
@test pm["b"] ≈ (4.5 / 13.5)
@test pm["c"] ≈ (3.5 / 13.5)

# testing small bits type
bx = [true, false, true, true, false]
@test countmap(bx) == Dict(true => 3, false => 2)

for T in [UInt8, UInt16, Int8, Int16]
    tx = T[typemin(T), 8, typemax(T), 19, 8]
    @test countmap(tx) == Dict(typemin(T) => 1, typemax(T) => 1, 8 => 2, 19 => 1)
end

@testset "views" begin
    X = view([1,1,1,2,2], 1:5)
    @test countmap(X) == countmap(copy(X))
end
