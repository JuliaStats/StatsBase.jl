using StatsBase
using Base.Test

### Trimming outliers

@test trim([1,2,3,4,5,6,7,8])           == [2,3,4,5,6,7]
@test trim([1,2,3,4,5,6,7,8], 0.1)      == [1,2,3,4,5,6,7,8]
@test trim([1,2,3,4,5,6,7,8], 0.2)      == [2,3,4,5,6,7]
@test trim([1,2,3,4,5,6,7,8,9], 0.4)    == [4,5,6]

@test_throws ArgumentError trim([])
@test_throws ArgumentError trim([1,2,3,4,5], 0.5)

@test trim!([1,2,3,4,5,6,7,8])          == [2,3,4,5,6,7]
@test trim!([1,2,3,4,5,6,7,8], 0.1)     == [1,2,3,4,5,6,7,8]
@test trim!([1,2,3,4,5,6,7,8], 0.2)     == [2,3,4,5,6,7]
@test trim!([1,2,3,4,5,6,7,8,9], 0.4)   == [4,5,6]

@test_throws ArgumentError trim!([])
@test_throws ArgumentError trim!([1,2,3,4,5], 0.5)

@test winsor([1,2,3,4,5,6,7,8])         == [2,2,3,4,5,6,7,7]
@test winsor([1,2,3,4,5,6,7,8], 0.1)    == [1,2,3,4,5,6,7,8]
@test winsor([1,2,3,4,5,6,7,8], 0.2)    == [2,2,3,4,5,6,7,7]
@test winsor([1,2,3,4,5,6,7,8,9], 0.4)  == [4,4,4,4,5,6,6,6,6]

@test_throws ArgumentError winsor([])
@test_throws ArgumentError winsor([1,2,3,4,5], 0.5)

@test winsor!([1,2,3,4,5,6,7,8])        == [2,2,3,4,5,6,7,7]
@test winsor!([1,2,3,4,5,6,7,8], 0.1)   == [1,2,3,4,5,6,7,8]
@test winsor!([1,2,3,4,5,6,7,8], 0.2)   == [2,2,3,4,5,6,7,7]
@test winsor!([1,2,3,4,5,6,7,8,9], 0.4) == [4,4,4,4,5,6,6,6,6]

@test_throws ArgumentError winsor!([])
@test_throws ArgumentError winsor!([1,2,3,4,5], 0.5)

### Variance

@test_approx_eq trimvar([1,1,1,1,1]) 0.0
@test_approx_eq trimvar([2,3,4,5,6,7,8,9], 0.25) 1.0

@test_throws ArgumentError trimvar([])
@test_throws ArgumentError trimvar([1,2,3,4,5], 0.5)