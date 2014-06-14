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

## middle

@test middle(1, 3) == 2
@test middle(1:6) == 3.5

## mode & modes

@test mode([1, 2, 3, 3, 2, 2, 1], 1:3) == 2
@test modes([1, 2, 3, 3, 2, 2, 1], 1:3) == [2]
@test modes([1, 3, 2, 3, 3, 2, 2, 1], 1:3) == [2, 3]

@test mode([1, 2, 3, 3, 2, 2, 1]) == 2
@test modes([1, 2, 3, 3, 2, 2, 1]) == [2]
@test sort(modes([1, 3, 2, 3, 3, 2, 2, 1])) == [2, 3]


###### quantile & friends

@test_approx_eq quantile(1:5) [1:5]
@test_approx_eq nquantile(1:5, 2) [1, 3, 5]
@test_approx_eq nquantile(1:5, 4) [1:5] 

@test_approx_eq percentile([1:5], 25) 2.0
@test_approx_eq percentile([1:5], [25, 50, 75]) [2.0, 3.0, 4.0]


##### Dispersion

@test_approx_eq variation([1:5]) 0.527046276694730

@test_approx_eq sem([1:5]) 0.707106781186548

@test_approx_eq mad([1:5], 3) 1.4826
@test_approx_eq mad([1:5]) 1.4826
@test_approx_eq mad(1:5) 1.4826

@test_approx_eq iqr(1:5) 2.0


##### moments

wv = weights(ones(5) * 2.0)

@test_approx_eq skewness(1:5) 0.0
@test_approx_eq skewness([1, 2, 3, 4, 5]) 0.0
@test_approx_eq skewness([1, 2, 2, 2, 5]) 1.1731251294063556
@test_approx_eq skewness([1, 4, 4, 4, 5]) -1.1731251294063556

@test_approx_eq skewness([1, 2, 2, 2, 5], wv) 1.1731251294063556

@test_approx_eq kurtosis(1:5) -1.3
@test_approx_eq kurtosis([1, 2, 3, 4, 5]) -1.3
@test_approx_eq kurtosis([1, 2, 3, 3, 2]) -1.1530612244897953

@test_approx_eq kurtosis([1, 2, 3, 4, 5], wv) -1.3


##### entropy

@test_approx_eq entropy([0.5, 0.5]) 0.6931471805599453
@test_approx_eq entropy([0.2, 0.3, 0.5]) 1.0296530140645737

@test_approx_eq crossentropy([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]) 1.1176681825904018

@test_approx_eq kldivergence([0.2, 0.3, 0.5], [0.3, 0.4, 0.3]) 0.08801516852582819


##### summarystats

s = summarystats(1:5)
@test isa(s, StatsBase.SummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test_approx_eq s.mean 3.0
@test_approx_eq s.median 3.0
@test_approx_eq s.q25 2.0
@test_approx_eq s.q75 4.0


