using StatsBase
using Base.Test

# test hist
@test sum(fit(Histogram,[1,2,3]).weights) == 3
@test fit(Histogram,Int[]).weights == Int[]
@test fit(Histogram,[1]).weights == [1]
@test fit(Histogram,[1,2,3],[0,2,4]) == Histogram([0,2,4],[2,1])
@test fit(Histogram,[1,2,3],[0,2,4]) != Histogram([0,2,4],[1,1])
@test fit(Histogram,[1,2,3],0:2:4) == Histogram(0:2:4,[2,1])
@test all(fit(Histogram,[1:100;]/100,0.0:0.01:1.0).weights .==1)
@test fit(Histogram,[1,1,1,1,1]).weights[1] == 5
@test sum(fit(Histogram,(rand(100),rand(100))).weights) == 100
@test fit(Histogram,1:100,nbins=5,closed=:right).weights == [20,20,20,20,20]
@test fit(Histogram,1:100,nbins=5,closed=:left).weights == [19,20,20,20,20,1]
@test fit(Histogram,0:99,nbins=5,closed=:right).weights == [1,20,20,20,20,19]
@test fit(Histogram,0:99,nbins=5,closed=:left).weights == [20,20,20,20,20]

@test fit(Histogram,(1:100,1:100),nbins=5).weights == diagm([20,20,20,20,20])
@test fit(Histogram,(1:100,1:100),nbins=(5,5)).weights == diagm([20,20,20,20,20])

@test fit(Histogram,1:100,weights(ones(100)),nbins=5).weights == [20,20,20,20,20]
@test fit(Histogram,1:100,weights(2*ones(100)),nbins=5).weights == [40,40,40,40,40]

@test eltype(fit(Histogram,1:100,weights(ones(Int,100)),nbins=5).weights) == Int
@test eltype(fit(Histogram,1:100,weights(ones(Float64,100)),nbins=5).weights) == Float64

import StatsBase.midpoints

@test midpoints(1.0:1.0:10.0) == 1.5:1.0:9.5
@test midpoints(1:10) == 1.5:9.5
@test midpoints(Float64[1.0:1.0:10.0;]) == Float64[1.5:1.0:9.5;]
