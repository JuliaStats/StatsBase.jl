using Stats
using Base.Test

## moments

@test_approx_eq skewness(1:5) 0.0
@test_approx_eq skewness([1, 2, 3, 4, 5]) 0.0
@test_approx_eq skewness([1, 2, 2, 2, 5]) 1.1731251294063556
@test_approx_eq skewness([1, 4, 4, 4, 5]) -1.1731251294063556

@test_approx_eq kurtosis(1:5) -1.3
@test_approx_eq kurtosis([1, 2, 3, 4, 5]) -1.3
@test_approx_eq kurtosis([1, 2, 3, 3, 2]) -1.1530612244897953

## variability

@test_approx_eq variation([1:5]) 0.527046276694730

@test_approx_eq sem([1:5]) 0.707106781186548

@test_approx_eq mad([1:5], 3) 1.4826
@test_approx_eq mad([1:5]) 1.4826
@test_approx_eq mad(1:5) 1.4826

# min/max related

@test minmax(1:5) == (1,5)
@test minmax([3, 2, 7, 1, 4]) == (1,7)
@test minmax([3.0, NaN, 1.0, NaN, 2.0]) == (1.0, 3.0)

@test_approx_eq midrange([1, 2, 3, 4, 5]) 3.0
@test_approx_eq midrange([1, 2, 2, 2, 5]) 3.0
@test_approx_eq midrange([1, 4, 4, 4, 5]) 3.0
@test_approx_eq midrange([1:5]) 3.0

@test_approx_eq range([1, 2, 3, 4, 5]) 4.0
@test_approx_eq range([1, 2, 2, 2, 5]) 4.0
@test_approx_eq range([1, 4, 4, 4, 5]) 4.0
@test_approx_eq range([1:15]) 14.0

# quantile & friends

@test_approx_eq quantile(1:5) [1:5]
@test_approx_eq nquantile(1:5, 2) [1, 3, 5]
@test_approx_eq nquantile(1:5, 4) [1:5] 

@test_approx_eq prctile([1:5], 25) 2.0
@test_approx_eq prctile([1:5], [25, 50, 75]) [2.0, 3.0, 4.0]
@test_approx_eq iqr(1:5) 2.0

# mode & modes

@test mode([1, 2, 3, 3, 2, 2, 1], 1:3) == 2
@test modes([1, 2, 3, 3, 2, 2, 1], 1:3) == [2]
@test modes([1, 3, 2, 3, 3, 2, 2, 1], 1:3) == [2, 3]

@test mode([1, 2, 3, 3, 2, 2, 1]) == 2
@test modes([1, 2, 3, 3, 2, 2, 1]) == [2]
@test sort(modes([1, 3, 2, 3, 3, 2, 2, 1])) == [2, 3]

# summarystats

s = summarystats(1:5)
@test isa(s, Stats.SummaryStats)
@test s.min == 1.0
@test s.max == 5.0
@test_approx_eq s.mean 3.0
@test_approx_eq s.median 3.0
@test_approx_eq s.q25 2.0
@test_approx_eq s.q75 4.0


