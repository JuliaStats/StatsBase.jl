# test hist
@test sum(hist([1,2,3]).weights) == 3
@test hist([]).weights == []
@test hist([1]).weights == [1]
@test hist([1,2,3],[0,2,4]) == Histogram([0,2,4],[2,1])
@test hist([1,2,3],[0,2,4]) != Histogram([0,2,4],[1,1])
@test hist([1,2,3],0:2:4) == Histogram(0:2:4,[2,1])
@test all(hist([1:100]/100,0.0:0.01:1.0).weights .==1)
@test hist([1,1,1,1,1]).weights[1] == 5
@test sum(hist((rand(100),rand(100))).weights) == 100
@test hist(1:100,5;closed=:right).weights == [20,20,20,20,20]
@test hist(1:100,5;closed=:left).weights == [19,20,20,20,20,1]
@test hist(0:99,5;closed=:right).weights == [1,20,20,20,20,19]
@test hist(0:99,5;closed=:left).weights == [20,20,20,20,20]

@test hist((1:100,1:100),5).weights == diagm([20,20,20,20,20])

@test hist(1:100,weights(ones(100)),5).weights == [20,20,20,20,20]
@test hist(1:100,weights(2*ones(100)),5).weights == [40,40,40,40,40]

@test eltype(hist(1:100,weights(ones(Int,100)),5).weights) == Int
@test eltype(hist(1:100,weights(ones(Float64,100)),5).weights) == Float64

import StatsBase.midpoints

@test midpoints(1.0:1.0:10.0) == 1.5:1.0:9.5
@test midpoints(1:10) == 1.5:9.5
@test midpoints(Float64[1.0:1.0:10.0]) == Float64[1.5:1.0:9.5]

