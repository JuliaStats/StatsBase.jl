using StatsBase
using Base.Test
using Distributions

##### Location

## geomean

@test_approx_eq geomean([1, 2, 3]) (6.0)^(1/3)
@test_approx_eq geomean(1:3) (6.0)^(1/3)
@test_approx_eq geomean([2, 8]) 4.0
@test_approx_eq geomean([4, 1, 1/32]) 0.5
@test geomean([1, 0, 2]) == 0.0

## harmmean

@test_approx_eq harmmean([1, 2, 3]) 3 / (1 + 1/2 + 1/3)
@test_approx_eq harmmean(1:3) 3 / (1 + 1/2 + 1/3)
@test_approx_eq harmmean([1, 2, 4]) 12 / 7

## trimmean

@test_approx_eq trimmean([-100, 2, 3, 7, 200], 0.0) 22.4
@test_approx_eq trimmean([-100, 2, 3, 7, 200], 0.4) 4.0
@test_approx_eq trimmean([-100, 2, 3, 7, 200], 0.8) 3.0
@test_approx_eq trimmean([2, 3, -100, 200, 7], 0.4) 4.0
@test_approx_eq trimmean([2, 3, -100, 200, 7], 0.8) 3.0

## mode & modes

@test mode([1, 2, 3, 3, 2, 2, 1], 1:3) == 2
@test modes([1, 2, 3, 3, 2, 2, 1], 1:3) == [2]
@test modes([1, 3, 2, 3, 3, 2, 2, 1], 1:3) == [2, 3]

@test mode([1, 2, 3, 3, 2, 2, 1]) == 2
@test modes([1, 2, 3, 3, 2, 2, 1]) == [2]
@test sort(modes([1, 3, 2, 3, 3, 2, 2, 1])) == [2, 3]

## zscores

@test zscore([-3:3;], 1.5, 0.5) == [-9.0:2.0:3.0;]

a = [3 4 5 6; 7 8 1 2; 6 9 3 0]
z1 = [4. 6. 8. 10.; 5. 6. -1. 0.; 1.5 3.0 0.0 -1.5]
z2 = [8. 2. 3. 1.; 24. 10. -1. -1.; 20. 12. 1. -2.]

@test_approx_eq zscore(a, [1, 2, 3], [0.5, 1.0, 2.0]) z1
@test_approx_eq zscore(a, [1 3 2 4], [0.25 0.5 1.0 2.0]) z2

@test zscore!(collect(-3.0:3.0), 1.5, 0.5) == [-9.0:2.0:3.0;]
@test_approx_eq zscore!(float(a), [1, 2, 3], [0.5, 1.0, 2.0]) z1
@test_approx_eq zscore!(float(a), [1 3 2 4], [0.25 0.5 1.0 2.0]) z2

@test zscore!(zeros(7), [-3:3;], 1.5, 0.5) == [-9.0:2.0:3.0;]
@test_approx_eq zscore!(zeros(size(a)), a, [1, 2, 3], [0.5, 1.0, 2.0]) z1
@test_approx_eq zscore!(zeros(size(a)), a, [1 3 2 4], [0.25 0.5 1.0 2.0]) z2

@test_approx_eq zscore(a) zscore(a, mean(a), std(a))
@test_approx_eq zscore(a, 1) zscore(a, mean(a,1), std(a,1))
@test_approx_eq zscore(a, 2) zscore(a, mean(a,2), std(a,2))


###### quantile & friends

@test_approx_eq quantile(1:5) [1:5;]
@test_approx_eq nquantile(1:5, 2) [1, 3, 5]
@test_approx_eq nquantile(1:5, 4) [1:5;]

@test_approx_eq percentile([1:5;], 25) 2.0
@test_approx_eq percentile([1:5;], [25, 50, 75]) [2.0, 3.0, 4.0]


##### Dispersion

@test StatsBase.span([3, 4, 5, 6, 2]) == (2:6)

@test_approx_eq variation([1:5;]) 0.527046276694730

@test_approx_eq sem([1:5;]) 0.707106781186548

@test_approx_eq mad([1:5;], 3) 1.4826
@test_approx_eq mad([1:5;]) 1.4826
@test_approx_eq mad(1:5) 1.4826
@test_approx_eq mad([1:5;], constant=1.0) 1.0
@test_approx_eq mad([1:5;], 3, constant=1.0) 1.0

@test_approx_eq iqr(1:5) 2.0


##### entropy

@test_approx_eq entropy([0.5, 0.5]) 0.6931471805599453
@test_approx_eq entropy([0.2, 0.3, 0.5]) 1.0296530140645737

@test_approx_eq entropy([0.5, 0.5],2) 1.0
@test_approx_eq entropy([0.2, 0.3, 0.5], 2) 1.4854752972273344

##### Renyi entropies
# Generate a random probability distribution of variable length
nindiv = rand(DiscreteUniform(2, 100))
dist = rand(Dirichlet(ones(nindiv)))

# Check Shannon entropy against Renyi entropy of order 1
@test_approx_eq entropy(dist) renyientropy(dist, 1)
@test_approx_eq renyientropy(dist, 1) renyientropy(dist, 1.0)

# Check Renyi entropy of order 0 is the natural log of the count of non-zeros
@test_approx_eq renyientropy(dist, 0) log(mapreduce(x -> x > 0 ? 1 : 0, +, dist))

# And is therefore not affected by the addition of non-zeros
zdist = dist
zdist = append!(dist, zeros(rand(DiscreteUniform(1, 100))))
@test_approx_eq renyientropy(dist, 0) renyientropy(zdist, 0)

# Indeed no Renyi entropy should be
order = rand(Uniform(0, 100))
@test_approx_eq renyientropy(dist, order) renyientropy(zdist, order)

# Renyi entropy of order infinity is -log(maximum(dist))
@test_approx_eq renyientropy(dist, Inf) -log(maximum(dist))

# And all Renyi entropies of uniform probabilities should be log n, where n is the length
udist = ones(nindiv) / nindiv
@test_approx_eq renyientropy(udist, order) log(nindiv)

##### Cross entropy
@test_approx_eq crossentropy([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]) 1.1176681825904018
@test_approx_eq crossentropy([0.2, 0.3, 0.5], [0.3, 0.4, 0.3], 2) 1.6124543443825532

##### KL divergence
@test_approx_eq kldivergence([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]) 0.08801516852582819
@test_approx_eq kldivergence([0.2, 0.3, 0.5], [0.3, 0.4, 0.3], 2) 0.12697904715521868

##### summarystats

s = summarystats(1:5)
@test isa(s, StatsBase.SummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test_approx_eq s.mean 3.0
@test_approx_eq s.median 3.0
@test_approx_eq s.q25 2.0
@test_approx_eq s.q75 4.0
