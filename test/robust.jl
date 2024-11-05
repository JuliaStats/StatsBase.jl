using StatsBase
using Test, Random

### Trimming outliers

@test collect(trim([8,2,3,4,5,6,7,1], prop=0.1)) == [8,2,3,4,5,6,7,1]
@test collect(trim([8,2,3,4,5,6,7,1], prop=0.2)) == [2,3,4,5,6,7]
@test collect(trim([1,2,3,4,5,6,7,8,9], prop=0.4)) == [4,5,6]
@test collect(trim([8,7,6,5,4,3,2,1], count=1)) == [7,6,5,4,3,2]
@test collect(trim([1,2,3,4,5,6,7,8,9], count=3)) == [4,5,6]


@test_throws ArgumentError trim([])
@test_throws ArgumentError trim([1,2,3,4,5], prop=0.5)

@test collect(trim!([8,2,3,4,5,6,7,1], prop=0.1)) == [8,2,3,4,5,6,7,1]
@test collect(trim!([8,2,3,4,5,6,7,1], prop=0.2)) == [2,3,4,5,6,7]
@test collect(trim!([1,2,3,4,5,6,7,8,9], prop=0.4)) == [4,5,6]
@test collect(trim!([8,7,6,5,4,3,2,1], count=1)) == [7,6,5,4,3,2]
@test collect(trim!([1,2,3,4,5,6,7,8,9], count=3)) == [4,5,6]

@test_throws ArgumentError trim!([])
@test_throws ArgumentError trim!([1,2,3,4,5], prop=0.5)

@test collect(winsor([8,2,3,4,5,6,7,1], prop=0.1)) == [8,2,3,4,5,6,7,1]
@test collect(winsor([8,2,3,4,5,6,7,1], prop=0.2)) == [7,2,3,4,5,6,7,2]
@test collect(winsor([1,2,3,4,5,6,7,8,9], prop=0.4)) == [4,4,4,4,5,6,6,6,6]
@test collect(winsor([1,2,3,4,5,6,7,8], count=1)) == [2,2,3,4,5,6,7,7]
@test collect(winsor([8,7,6,5,4,3,2,1], count=1)) == [7,7,6,5,4,3,2,2]
@test collect(winsor([1,2,3,4,5,6,7,8,9], count=3)) == [4,4,4,4,5,6,6,6,6]

@test_throws ArgumentError winsor([])
@test_throws ArgumentError winsor([1,2,3,4,5], prop=0.5)

@test collect(winsor!([8,2,3,4,5,6,7,1], prop=0.1)) == [8,2,3,4,5,6,7,1]
@test collect(winsor!([8,2,3,4,5,6,7,1], prop=0.2)) == [7,2,3,4,5,6,7,2]
@test collect(winsor!([1,2,3,4,5,6,7,8,9], prop=0.4)) == [4,4,4,4,5,6,6,6,6]
@test collect(winsor!([8,7,6,5,4,3,2,1], count=1)) == [7,7,6,5,4,3,2,2]
@test collect(winsor!([1,2,3,4,5,6,7,8,9], count=3)) == [4,4,4,4,5,6,6,6,6]

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

Random.seed!(1234)
for n in 2100:2120, c in 0:1000
    x = randperm(n)
    @test sort!(collect(winsor(x, count=c))) ==
          reverse!(collect(winsor(n:-1:1, count=c))) ==
          collect(winsor(1:n, count=c))
end
