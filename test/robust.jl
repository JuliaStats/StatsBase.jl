using StatsBase
using Test

### Trimming outliers

@test trim([1,2,3,4,5,6,7,8],   prop=0.1) == [1,2,3,4,5,6,7,8]
@test trim([1,2,3,4,5,6,7,8],   prop=0.2) == [2,3,4,5,6,7]
@test trim([1,2,3,4,5,6,7,8,9], prop=0.4) == [4,5,6]
@test trim([1,2,3,4,5,6,7,8],   count=1)  == [2,3,4,5,6,7]
@test trim([1,2,3,4,5,6,7,8,9], count=3)  == [4,5,6]


@test_throws ArgumentError trim([])
@test_throws ArgumentError trim([1,2,3,4,5], prop=0.5)

@test trim!([1,2,3,4,5,6,7,8],   prop=0.1) == [1,2,3,4,5,6,7,8]
@test trim!([1,2,3,4,5,6,7,8],   prop=0.2) == [2,3,4,5,6,7]
@test trim!([1,2,3,4,5,6,7,8,9], prop=0.4) == [4,5,6]
@test trim!([1,2,3,4,5,6,7,8],   count=1)  == [2,3,4,5,6,7]
@test trim!([1,2,3,4,5,6,7,8,9], count=3)  == [4,5,6]

@test_throws ArgumentError trim!([])
@test_throws ArgumentError trim!([1,2,3,4,5], prop=0.5)

@test winsor([1,2,3,4,5,6,7,8],   prop=0.1) == [1,2,3,4,5,6,7,8]
@test winsor([1,2,3,4,5,6,7,8],   prop=0.2) == [2,2,3,4,5,6,7,7]
@test winsor([1,2,3,4,5,6,7,8,9], prop=0.4) == [4,4,4,4,5,6,6,6,6]
@test winsor([1,2,3,4,5,6,7,8],   count=1)  == [2,2,3,4,5,6,7,7]
@test winsor([1,2,3,4,5,6,7,8,9], count=3)  == [4,4,4,4,5,6,6,6,6]

@test_throws ArgumentError winsor([])
@test_throws ArgumentError winsor([1,2,3,4,5], prop=0.5)

@test winsor!([1,2,3,4,5,6,7,8],   prop=0.1) == [1,2,3,4,5,6,7,8]
@test winsor!([1,2,3,4,5,6,7,8],   prop=0.2) == [2,2,3,4,5,6,7,7]
@test winsor!([1,2,3,4,5,6,7,8,9], prop=0.4) == [4,4,4,4,5,6,6,6,6]
@test winsor!([1,2,3,4,5,6,7,8],   count=1)  == [2,2,3,4,5,6,7,7]
@test winsor!([1,2,3,4,5,6,7,8,9], count=3)  == [4,4,4,4,5,6,6,6,6]

@test_throws ArgumentError winsor!([])
@test_throws ArgumentError winsor!([1,2,3,4,5], prop=0.5)

### Variance

@test trimvar([1,1,1,1,1]) ≈ 0.0
@test trimvar([2,3,4,5,6,7,8,9], prop=0.25) ≈ 1.0

@test_throws ArgumentError trimvar([])
@test_throws ArgumentError trimvar([1,2,3,4,5], prop=0.5)

### Other

@test mean(trim([-Inf,1,2,3,4], count=1)) == 2
@test mean(winsor([-Inf,1,2,3,4], count=1)) == 2
