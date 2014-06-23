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

x = float64([2:8])
@test_approx_eq moment(x, 2) sum((x .- 5).^2) / length(x)
@test_approx_eq moment(x, 3) sum((x .- 5).^3) / length(x)
@test_approx_eq moment(x, 4) sum((x .- 5).^4) / length(x)
@test_approx_eq moment(x, 5) sum((x .- 5).^5) / length(x)

@test_approx_eq moment(x, 2, 4.0) sum((x .- 4).^2) / length(x)
@test_approx_eq moment(x, 3, 4.0) sum((x .- 4).^3) / length(x)
@test_approx_eq moment(x, 4, 4.0) sum((x .- 4).^4) / length(x)
@test_approx_eq moment(x, 5, 4.0) sum((x .- 4).^5) / length(x)

w = weights([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
x2 = float64([2:6])
@test_approx_eq moment(x, 2, w) sum((x2 .- 4).^2) / 5
@test_approx_eq moment(x, 3, w) sum((x2 .- 4).^3) / 5
@test_approx_eq moment(x, 4, w) sum((x2 .- 4).^4) / 5
@test_approx_eq moment(x, 5, w) sum((x2 .- 4).^5) / 5


##### weighted var & std

x = rand(10)
wv = weights(rand(10))
m = mean(x, wv)

@test_approx_eq var(x, wv) sum(abs2(x .- m), wv) ./ sum(wv)
@test_approx_eq var(x, wv; mean=0) sum(abs2(x), wv) ./ sum(wv)
@test_approx_eq var(x, wv; mean=1.0) sum(abs2(x .- 1.0), wv) ./ sum(wv)

@test_approx_eq std(x, wv) sqrt(var(x, wv))
@test_approx_eq std(x, wv; mean=0) sqrt(var(x, wv; mean=0))
@test_approx_eq std(x, wv; mean=1.0) sqrt(var(x, wv; mean=1.0))


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


