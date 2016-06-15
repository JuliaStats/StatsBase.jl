using StatsBase
using Base.Test

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

@test span([3, 4, 5, 6, 2]) == (2:6)

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

@test_approx_eq crossentropy([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]) 1.1176681825904018
@test_approx_eq crossentropy([0.2, 0.3, 0.5], [0.3, 0.4, 0.3], 2) 1.6124543443825532

@test_approx_eq kldivergence([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]) 0.08801516852582819
@test_approx_eq kldivergence([0.2, 0.3, 0.5], [0.3, 0.4, 0.3], 2) 0.12697904715521868

##### summarystats

s = describe(1:5)
@show s
@test isa(s, StatsBase.SummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test_approx_eq s.mean 3.0
@test_approx_eq s.median 3.0
@test_approx_eq s.p25 2.0
@test_approx_eq s.p75 4.0

s = describe(1:5, detail = true)
@show s
@test isa(s, StatsBase.DetailedSummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test_approx_eq s.mean 3.0
@test_approx_eq s.median 3.0
@test_approx_eq s.p25 2.0
@test_approx_eq s.p75 4.0
@test_approx_eq s.p1 1.04
@test_approx_eq s.p5 1.2
@test_approx_eq s.p95 4.8
@test_approx_eq s.p99 4.96
@test_approx_eq s.skewness 0.0
@test_approx_eq s.kurtosis -1.3

s = describe(1:5, weights(1:5))
@show s
@test isa(s, StatsBase.WeightedSummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test_approx_eq s.mean  3.6666666666666
@test_approx_eq s.median 3.6666666666666
@test_approx_eq s.p25 2.750000
@test_approx_eq s.p75 4.375000

s = describe(1:5, weights(1:5), detail = true)
@show s
@test isa(s, StatsBase.WeightedDetailedSummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test_approx_eq s.mean  3.6666666666666
@test_approx_eq s.median 3.6666666666666
@test_approx_eq s.p25 2.750000
@test_approx_eq s.p75 4.375000
@test_approx_eq s.p1  1.100000
@test_approx_eq s.p5 1.500000
@test_approx_eq s.p95 4.875000
@test_approx_eq s.p99 4.975000
@test_approx_eq s.skewness -0.5879747322073332
@test_approx_eq s.kurtosis  -0.7285714285714291


