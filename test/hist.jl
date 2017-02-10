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

# histrange
# Note: atm histrange must be qualified
@test StatsBase.histrange(Float64[], 0, :left) == 0.0:1.0:0.0
@test StatsBase.histrange(Float64[1:5;], 1, :left) == 0.0:5.0:10.0
@test StatsBase.histrange(Float64[1:10;], 1, :left) == 0.0:10.0:20.0
@test StatsBase.histrange(1.0, 10.0, 1, :left) == 0.0:10.0:20.0

@test StatsBase.histrange([0.201,0.299], 10, :left) == 0.2:0.01:0.3
@test StatsBase.histrange([0.2,0.299], 10, :left) == 0.2:0.01:0.3
@test StatsBase.histrange([0.2,0.3], 10, :left)  == 0.2:0.01:0.31
@test StatsBase.histrange(0.2, 0.3,  10, :left)  == 0.2:0.01:0.31
@test StatsBase.histrange([0.2,0.3], 10, :right) == 0.19:0.01:0.3
@test StatsBase.histrange(0.2, 0.3,  10, :right) == 0.19:0.01:0.3

@test StatsBase.histrange([200.1,299.9], 10, :left) == 200.0:10.0:300.0
@test StatsBase.histrange([200.0,299.9], 10, :left) == 200.0:10.0:300.0
@test StatsBase.histrange([200.0,300.0], 10, :left) == 200.0:10.0:310.0
@test StatsBase.histrange([200.0,300.0], 10, :right) == 190.0:10.0:300.0

@test StatsBase.histrange(Int64[1:5;], 1, :left) == 0:5:10
@test StatsBase.histrange(Int64[1:10;], 1, :left) == 0:10:20

@test StatsBase.histrange([1, 2, 3, 4], 4) == 0.0:1.0:4.0
@test StatsBase.histrange([1, 2, 2, 4], 4) == 0.0:1.0:4.0
@test StatsBase.histrange([1, 10], 4) == 0.0:5.0:10.0
@test StatsBase.histrange([1, 20], 4) == 0.0:5.0:20.0
@test StatsBase.histrange([1, 600], 4) == 0.0:200.0:600.0
@test StatsBase.histrange([1, -1000], 4) == -1500.0:500.0:500.0

# Base issue #13326
l,h = extrema(StatsBase.histrange([typemin(Int),typemax(Int)], 10))
@test l <= typemin(Int)
@test h >= typemax(Int)

@test_throws ArgumentError StatsBase.histrange([1, 10], 0)
@test_throws ArgumentError StatsBase.histrange([1, 10], -1)
@test_throws ArgumentError StatsBase.histrange([1.0, 10.0], 0)
@test_throws ArgumentError StatsBase.histrange([1.0, 10.0], -1)
@test_throws ArgumentError StatsBase.histrange(Float64[],-1)
@test_throws ArgumentError StatsBase.histrange([0.], 0)


# hist show
show_h = sprint(show, fit(Histogram,[1,2,3]))
@test contains(show_h, "edges:\n  0.0:1.0:3.0")
@test contains(show_h, "weights: $([1,1,1])")
@test contains(show_h, "closed: right")
