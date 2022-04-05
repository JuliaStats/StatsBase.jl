using StatsBase
using Test
using DelimitedFiles
using Statistics

##### Location

## geomean

@test geomean([1, 2, 3])    ≈ cbrt(6.0)
@test geomean(1:3)          ≈ cbrt(6.0)
@test geomean([2, 8])       ≈ 4.0
@test geomean([4, 1, 1/32]) ≈ 0.5
@test geomean([1, 0, 2]) == 0.0

## harmmean

@test harmmean([1, 2, 3]) ≈ 3 / (1 + 1/2 + 1/3)
@test harmmean(1:3)       ≈ 3 / (1 + 1/2 + 1/3)
@test harmmean([1, 2, 4]) ≈ 12 / 7

## genmean
@test genmean([1,1,2,3], 1)         ≈ 7/4
@test genmean([1,4,2], -1)          ≈ 12/7
@test genmean([1,1,2,3], 0)         ≈ (6.0)^(1/4)
@test genmean([1.2,-0.5,0], 2)      ≈ sqrt(169/300)
@test genmean([16/9,0.25,1.0], 1.5) ≈ (755/648)^(2/3)
# Test numerical stability for `p` close to 0 (genmean should be close to geometric mean).
@test isapprox(genmean([1,1,2,3], -1e-8), (6.0)^(1/4), atol=1e-8)
# Test numerical stability for large `p` (genmean should be close to max).
@test isapprox(genmean([0.98,1.02], 1e4), 1.02, atol=1e-4)

## mode & modes

@test mode([1, 2, 3, 3, 2, 2, 1], 1:3) == 2
@test modes([1, 2, 3, 3, 2, 2, 1], 1:3) == [2]
@test modes([1, 3, 2, 3, 3, 2, 2, 1], 1:3) == [2, 3]

@test mode([1, 2, 3, 3, 2, 2, 1]) == 2
@test modes([1, 2, 3, 3, 2, 2, 1]) == [2]
@test sort(modes([1, 3, 2, 3, 3, 2, 2, 1])) == [2, 3]

@test mode(skipmissing([1, missing, missing, 3, 2, 2, missing])) == 2
@test modes(skipmissing([1, missing, missing, 3, 2, 2, missing])) == [2]
@test sort(modes(skipmissing([1, missing, 3, 3, 2, 2, missing]))) == [2, 3]

d1 = [1, 2, 3, 3, 4, 5, 5, 3]
d2 = ['a', 'b', 'c', 'c', 'd', 'e', 'e', 'c']
wv = weights([0.1:0.1:0.7; 0.1])
@test mode(d1) == 3
@test mode(d2) == 'c'
@test mode(d1, wv) == 5
@test mode(d2, wv) == 'e'
@test sort(modes(d1[1:end-1], weights(ones(7)))) == [3, 5]
@test sort(modes(d1, weights([.9, .1, .1, .1, .9, .1, .1, .1]))) == [1, 4]

@test_throws ArgumentError mode(Int[])
@test_throws ArgumentError modes(Int[])
@test_throws ArgumentError mode(Any[])
@test_throws ArgumentError modes(Any[])
@test_throws ArgumentError mode([], weights(Float64[]))
@test_throws ArgumentError modes([], weights(Float64[]))
@test_throws ArgumentError mode([1, 2, 3], weights([0.1, 0.3]))
@test_throws ArgumentError modes([1, 2, 3], weights([0.1, 0.3]))

## zscores

@test zscore([-3:3;], 1.5, 0.5) == [-9.0:2.0:3.0;]

a = [3 4 5 6; 7 8 1 2; 6 9 3 0]
z1 = [4. 6. 8. 10.; 5. 6. -1. 0.; 1.5 3.0 0.0 -1.5]
z2 = [8. 2. 3. 1.; 24. 10. -1. -1.; 20. 12. 1. -2.]

@test zscore(a, [1, 2, 3], [0.5, 1.0, 2.0])    ≈ z1
@test zscore(a, [1 3 2 4], [0.25 0.5 1.0 2.0]) ≈ z2

@test zscore!(collect(-3.0:3.0), 1.5, 0.5) == [-9.0:2.0:3.0;]
@test zscore!(float(a), [1, 2, 3], [0.5, 1.0, 2.0])    ≈ z1
@test zscore!(float(a), [1 3 2 4], [0.25 0.5 1.0 2.0]) ≈ z2

@test zscore!(zeros(7), [-3:3;], 1.5, 0.5) == [-9.0:2.0:3.0;]
@test zscore!(zeros(size(a)), a, [1, 2, 3], [0.5, 1.0, 2.0])    ≈ z1
@test zscore!(zeros(size(a)), a, [1 3 2 4], [0.25 0.5 1.0 2.0]) ≈ z2

@test zscore(a)    ≈ zscore(a, mean(a), std(a))
@test zscore(a, 1) ≈ zscore(a, mean(a, dims=1), std(a, dims=1))
@test zscore(a, 2) ≈ zscore(a, mean(a, dims=2), std(a, dims=2))


###### quantile & friends

@test nquantile(1:5, 2) ≈ [1, 3, 5]
@test nquantile(1:5, 4) ≈ [1:5;]
@test nquantile(skipmissing([missing, 2, 5, missing]), 2) ≈ [2.0, 3.5, 5.0]

@test percentile([1:5;], 25)           ≈  2.0
@test percentile([1:5;], [25, 50, 75]) ≈ [2.0, 3.0, 4.0]
@test percentile(skipmissing([missing, 2, 5, missing]), 25) ≈ 2.75
@test percentile(skipmissing([missing, 2, 5, missing]), [25, 50, 75]) ≈ [2.75, 3.5, 4.25]

@testset "quantilerank and percentilerank" begin
     @testset "value as number and array" begin
         @testset ":inc and :exc" begin
             v1 = [1, 1, 1, 2, 3, 4, 8, 11, 12, 13]
             v2 = [1, 2, 3, 6, 6, 6, 7, 8, 9]
             v3 = [1, 2, 4, 3, 4]
             v4 = [1, 2, 1, 3, 4]
             @test quantilerank(v1, 2, method=:inc)    == 1/3
             @test quantilerank(v1, 4, method=:inc)    == 5/9
             @test quantilerank(v1, 8, method=:inc)    == 2/3
             @test quantilerank(v1, 5, method=:inc)    == 7/12        
             @test quantilerank(v2, 7, method=:exc)    == 0.7
             @test quantilerank(v2, 5.43, method=:exc) == 0.381
             @test quantilerank(v3, 4, method=:exc)    == 6/9
             @test quantilerank(v3, 4, method=:inc)    == 3/4
             @test quantilerank(v4, 1, method=:exc)    == 1/6
             @test quantilerank(v4, -100, method=:inc) == 0.0
             @test quantilerank(v4,  100, method=:inc) == 1.0
             @test quantilerank(v4, -100, method=:exc) == 0.0
             @test quantilerank(v4,  100, method=:exc) == 1.0
             @test percentilerank(v1, 2)               == 100 * quantilerank(v1, 2)
             @test percentilerank(v2, 7, method=:exc)  == 100 * quantilerank(v2, 7, method=:exc)
         end
         @testset ":compete" begin
             v = [0, 0, 1, 1, 2, 2, 2, 2, 4, 4]
             @test quantilerank(v, 1, method=:compete)    == 2/9
             @test quantilerank(v, 2, method=:compete)    == 4/9
             @test quantilerank(v, 4, method=:compete)    == 8/9
             @test quantilerank(v, -100, method=:compete) == 0.0
             @test quantilerank(v,  100, method=:compete) == 1.0
         end
         @testset ":strict, :weak and :tied" begin
             v = [7, 8, 2, 1, 3, 4, 5, 4, 6, 9]
             for (method, res1, res2) in [(:tied, .4, [.4, .85]),
                                          (:strict, .3, [.3, .8]),
                                          (:weak, .5, [.5, .9])]
                 @test quantilerank(v, 4, method=method) == res1
             end
         end
     end
     @testset "errors" begin
         v1 = [1, 2, 3, 5, 6, missing, 8]
         v2 = [missing, missing]
         v3 = [1, 2, 3, 5, 6, NaN, 8]
         v4 = [1, 2, 3, 3, 4]
         for method in (:tied, :strict, :weak)
             @test_throws ArgumentError quantilerank(v1, 4, method=method)
             @test_throws ArgumentError quantilerank(v2, 4, method=method)
             @test_throws ArgumentError quantilerank(v3, 4, method=method)
         end
         @test_throws ArgumentError quantilerank(v4, 3, method=:wrongargument)
         @test_throws ArgumentError quantilerank(v4, NaN)
         @test_throws ArgumentError quantilerank(v4, missing)
         @test_throws ArgumentError quantilerank([], 3)
         @test_throws ArgumentError quantilerank([1], 3)
     end
 end
 
##### Dispersion

@test span([3, 4, 5, 6, 2]) == (2:6)
@test span(skipmissing([1, missing, 5, missing])) == 1:5

@test variation([1:5;]) ≈ 0.527046276694730
@test variation(skipmissing([missing; 1:5; missing])) ≈ 0.527046276694730

@test @inferred(sem([1:5;])) ≈ 0.707106781186548
@test @inferred(sem(skipmissing([missing; 1:5; missing]))) ≈ 0.707106781186548
@test @inferred(sem(skipmissing([missing; 1:5; missing]), mean=3.0)) ≈ 0.707106781186548
@test @inferred(sem([1:5;], UnitWeights{Int}(5))) ≈ 0.707106781186548
@test @inferred(sem([1:5;], UnitWeights{Int}(5); mean=mean(1:5))) ≈ 0.707106781186548
@test_throws DimensionMismatch sem(1:5, UnitWeights{Int}(4))
@test @inferred(sem([1:5;], ProbabilityWeights([1:5;]))) ≈ 0.6166 rtol=.001
μ = mean(1:5, ProbabilityWeights([1:5;]))
@test @inferred(sem([1:5;], ProbabilityWeights([1:5;]); mean=μ)) ≈ 0.6166 rtol=.001
@test @inferred(sem([10; 1:5;], ProbabilityWeights([0; 1:5;]); mean=μ)) ≈ 0.6166 rtol=.001
x = sort!(vcat([5:-1:i for i in 1:5]...))
μ = mean(x)
@test @inferred(sem([1:5;], FrequencyWeights([1:5;]))) ≈ sem(x)
@test @inferred(sem([1:5;], FrequencyWeights([1:5;]); mean=μ)) ≈ sem(x)

@inferred sem([1:5f0;]; mean=μ) ≈ sem(x)
@inferred sem([1:5f0;], ProbabilityWeights([1:5;]); mean=μ)
@inferred sem([1:5f0;], FrequencyWeights([1:5;]); mean=μ)
# Broken: Bug to do with Statistics.jl's implementation of `var`
# @inferred sem([1:5f0;], UnitWeights{Int}(5); mean=μ)

@test @inferred(isnan(sem(Int[])))
@test @inferred(isnan(sem(Int[], FrequencyWeights(Int[]))))
@test @inferred(isnan(sem(Int[], ProbabilityWeights(Int[]))))

@test @inferred(isnan(sem(Int[]; mean=0f0)))
@test @inferred(isnan(sem(Int[], FrequencyWeights(Int[]); mean=0f0)))
@test @inferred(isnan(sem(Int[], ProbabilityWeights(Int[]); mean=0f0)))

@test @inferred(isnan(sem(skipmissing(Union{Int,Missing}[missing, missing]))))
@test_throws MethodError sem(Any[])
@test_throws MethodError sem(skipmissing([missing]))

@test mad(1:5; center=3, normalize=true) ≈ 1.4826022185056018
@test mad(skipmissing([missing; 1:5; missing]); center=3, normalize=true) ≈ 1.4826022185056018
@test StatsBase.mad!([1:5;]; center=3, normalize=true) ≈ 1.4826022185056018
@test mad(1:5, normalize=true) ≈ 1.4826022185056018
@test mad(1:5, normalize=false) ≈ 1.0
@test mad(skipmissing([missing; 1:5; missing]), normalize=true) ≈ 1.4826022185056018
@test mad(skipmissing([missing; 1:5; missing]), normalize=false) ≈ 1.0
@test StatsBase.mad!([1:5;], normalize=false) ≈ 1.0
@test mad(1:5, center=3, normalize=false) ≈ 1.0
@test mad(skipmissing([missing; 1:5; missing]), center=3, normalize=false) ≈ 1.0
@test StatsBase.mad!([1:5;], center=3, normalize=false) ≈ 1.0
@test mad((x for x in (1, 2.1)), normalize=false) ≈ 0.55
@test mad(Any[1, 2.1], normalize=false) ≈ 0.55
@test mad(Union{Int,Missing}[1, 2], normalize=false) ≈ 0.5
@test_throws ArgumentError mad(Int[], normalize = true)

# Issue 197
@test mad(1:2, normalize=true) ≈ 0.7413011092528009

@test iqr(1:5) ≈ 2.0

nutrient = readdlm(joinpath(@__DIR__, "data", "nutrient.txt"))[:,2:end]
@test @inferred(genvar(nutrient)) ≈ 2.8310418e19 rtol=1e-6
@test @inferred(totalvar(nutrient)) ≈ 2.83266877e6 rtol=1e-6

X = [1 2 5
     4 1 6
     4 0 4]
@test @inferred(genvar(X)) ≈ 0.0
@test @inferred(totalvar(X)) ≈ 5.0

x = rand(Float32, 10)
@test genvar(x) == totalvar(x) == var(x)

it = (xᵢ for xᵢ in x)
@test genvar(it) == totalvar(it) == var(it)


##### entropy

@test @inferred(entropy([0.5, 0.5]))      ≈ 0.6931471805599453
@test @inferred(entropy([1//2, 1//2]))    ≈ 0.6931471805599453
@test @inferred(entropy([0.5f0, 0.5f0])) isa Float32
@test @inferred(entropy([0.2, 0.3, 0.5])) ≈ 1.0296530140645737
@test iszero(@inferred(entropy([0, 1])))
@test iszero(@inferred(entropy([0.0, 1.0])))

@test @inferred(entropy([0.5, 0.5], 2))      ≈ 1.0
@test @inferred(entropy([1//2, 1//2], 2))    ≈ 1.0
@test @inferred(entropy([0.2, 0.3, 0.5], 2)) ≈ 1.4854752972273344

@test_throws ArgumentError @inferred(entropy(Float64[]))
@test_throws ArgumentError @inferred(entropy(Int[]))

##### Renyi entropies
# Generate a random probability distribution
nindiv = 50
dist = rand(nindiv)
dist /= sum(dist)

# Check Shannon entropy against Renyi entropy of order 1
@test entropy(dist)         ≈ renyientropy(dist, 1)
@test renyientropy(dist, 1) ≈ renyientropy(dist, 1.0)

# Check Renyi entropy of order 0 is the natural log of the count of non-zeros
@test renyientropy(dist, 0) ≈ log(mapreduce(x -> x > 0 ? 1 : 0, +, dist))

# And is therefore not affected by the addition of non-zeros
zdist = dist
zdist = append!(dist, zeros(50))
@test renyientropy(dist, 0) ≈ renyientropy(zdist, 0)

# Indeed no Renyi entropy should be
loworder = rand() # Between 0 and 1
@test renyientropy(dist, loworder) ≈ renyientropy(zdist, loworder)
order = rand() * 49 + 1 # And over 1
@test renyientropy(dist, order) ≈ renyientropy(zdist, order)

# Renyi entropy of order infinity is -log(maximum(dist))
@test renyientropy(dist, Inf) ≈ -log(maximum(dist))

# And all Renyi entropies of uniform probabilities should be log n, where n is the length
udist = ones(nindiv) / nindiv
@test renyientropy(udist, order) ≈ log(nindiv)

# And test generalised probability distributions (sum(p) != 1)
scale = rand()
@test renyientropy(udist * scale, 0)     ≈ renyientropy(udist, 0) - log(scale)
@test renyientropy(udist * scale, 1)     ≈ renyientropy(udist, 1) - log(scale)
@test renyientropy(udist * scale, Inf)   ≈ renyientropy(udist, Inf) - log(scale)
@test renyientropy(udist * scale, order) ≈ renyientropy(udist, order) - log(scale)

##### Cross entropy
@test @inferred(crossentropy([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]))     ≈ 1.1176681825904018
@test @inferred(crossentropy([1//5, 3//10, 1//2], [0.3, 0.4, 0.3])) ≈ 1.1176681825904018
@test @inferred(crossentropy([1//5, 3//10, 1//2], [0.3f0, 0.4f0, 0.3f0])) isa Float32
@test @inferred(crossentropy([0.2, 0.3, 0.5], [0.3, 0.4, 0.3], 2))     ≈ 1.6124543443825532
@test @inferred(crossentropy([1//5, 3//10, 1//2], [0.3, 0.4, 0.3], 2)) ≈ 1.6124543443825532
@test @inferred(crossentropy([1//5, 3//10, 1//2], [0.3f0, 0.4f0, 0.3f0], 2f0)) isa Float32

# deprecated, should throw an `ArgumentError` at some point
logpattern = (:warn, "support for empty collections will be removed since they do not represent proper probability distributions")
@test iszero(@test_logs logpattern @inferred(crossentropy(Float64[], Float64[])))
@test iszero(@test_logs logpattern @inferred(crossentropy(Int[], Int[])))

##### KL divergence
@test @inferred(kldivergence([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]))     ≈ 0.08801516852582819
@test @inferred(kldivergence([1//5, 3//10, 1//2], [0.3, 0.4, 0.3])) ≈ 0.08801516852582819
@test @inferred(kldivergence([1//5, 3//10, 1//2], [0.3f0, 0.4f0, 0.3f0])) isa Float32
@test @inferred(kldivergence([0.2, 0.3, 0.5], [0.3, 0.4, 0.3], 2))     ≈ 0.12697904715521868
@test @inferred(kldivergence([1//5, 3//10, 1//2], [0.3, 0.4, 0.3], 2)) ≈ 0.12697904715521868
@test @inferred(kldivergence([1//5, 3//10, 1//2], [0.3f0, 0.4f0, 0.3f0], 2f0)) isa Float32
@test iszero(@inferred(kldivergence([0, 1], [0f0, 1f0])))

# deprecated, should throw an `ArgumentError` at some point
logpattern = (:warn, "support for empty collections will be removed since they do not represent proper probability distributions")
@test iszero(@test_logs logpattern @inferred(kldivergence(Float64[], Float64[])))
@test iszero(@test_logs logpattern @inferred(kldivergence(Int[], Int[])))

##### summarystats

s = summarystats(1:5)
@test isa(s, StatsBase.SummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test s.mean   ≈ 3.0
@test s.median ≈ 3.0
@test s.q25    ≈ 2.0
@test s.q75    ≈ 4.0

# Issue #631
s = summarystats([-2, -1, 0, 1, 2, missing])
@test isa(s, StatsBase.SummaryStats)
@test s.min == -2.0
@test s.max == 2.0
@test s.mean   ≈ 0.0
@test s.median ≈ 0.0
@test s.q25    ≈ -1.0
@test s.q75    ≈ +1.0

# Issue #631
s = summarystats(zeros(10))
@test isa(s, StatsBase.SummaryStats)
@test s.min == 0.0
@test s.max == 0.0
@test s.mean   ≈ 0.0
@test s.median ≈ 0.0
@test s.q25    ≈ 0.0
@test s.q75    ≈ 0.0

# Issue #631
s = summarystats(Union{Float64,Missing}[missing, missing])
@test isa(s, StatsBase.SummaryStats)
@test s.nobs == 2
@test s.nmiss == 2
@test isnan(s.mean)
@test isnan(s.median)
