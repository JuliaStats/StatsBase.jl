using StatsBase
using Base.Test

# Should this be 4 - 2?
# @test_approx_eq iqr([1, 2, 3, 4, 5])

@test_approx_eq kurtosis([1, 2, 3, 4, 5]) -1.3
@test_approx_eq kurtosis([1, 2, 3, 3, 2]) -1.1530612244897953

@test_approx_eq mad([1, 2, 3, 4, 5]) 1.4826
# TODO: This shouldn't fail
# @test_approx_eq mad(1:5) 1.4826

@test_approx_eq midrange([1, 2, 3, 4, 5]) 3.0
@test_approx_eq midrange([1, 2, 2, 2, 5]) 3.0
@test_approx_eq midrange([1, 4, 4, 4, 5]) 3.0
@test_approx_eq midrange([1:5]) 3.0

@test_approx_eq range([1, 2, 3, 4, 5]) 4.0
@test_approx_eq range([1, 2, 2, 2, 5]) 4.0
@test_approx_eq range([1, 4, 4, 4, 5]) 4.0
@test_approx_eq range([1:15]) 14.0

@test_approx_eq skewness([1, 2, 3, 4, 5]) 0.0
@test_approx_eq skewness([1, 2, 2, 2, 5]) 1.1731251294063556
@test_approx_eq skewness([1, 4, 4, 4, 5]) -1.1731251294063556

# sem()

# variation()
