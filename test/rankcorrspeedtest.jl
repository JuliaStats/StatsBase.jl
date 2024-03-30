using Random, StatsBase, Dates, BenchmarkTools

using Random;
x = rand(MersenneTwister(0), 1000, 10);
xm = ifelse.(x .< 0.05, missing, x);

#compile...
res_1 = corkendall(x)
res_2 = corkendall(xm; skipmissing=:pairwise)
res_3 = pairwise(corkendall, eachcol(xm); skipmissing=:pairwise)
res_4 = corkendall(xm; skipmissing=:listwise)
res_5 = pairwise(corkendall, eachcol(xm); skipmissing=:listwise)
res_6 = corkendall(xm; skipmissing=:none)
res_7 = pairwise(corkendall, eachcol(xm), skipmissing=:none)
res_8 = corspearman(x)
res_9 = corspearman(xm; skipmissing=:pairwise)
res_10 = pairwise(corspearman, eachcol(xm); skipmissing=:pairwise)
res_11 = corspearman(xm; skipmissing=:listwise)
res_12 = pairwise(corspearman, eachcol(xm); skipmissing=:listwise)
res_13 = corspearman(xm; skipmissing=:none)
res_14 = pairwise(corspearman, eachcol(xm), skipmissing=:none)

@assert res_2 == res_3
@assert res_4 == res_5
@assert res_9 == res_10
@assert res_11 ≈ res_12

x = rand(MersenneTwister(0), 1000, 1000);
xm = ifelse.(x .< 0.05, missing, x);

println("="^130)
@show(Dates.now())
@show ENV["COMPUTERNAME"]
println("Julia Version $VERSION")
@show Threads.nthreads()
@show size(x)
@show typeof(x)
@show size(xm)
@show typeof(xm)

print("corkendall(\$x)                                           ")
@btime res_1 = corkendall($x);
print("corkendall(\$xm; skipmissing=:pairwise)                   ")
@btime res_2 = corkendall($xm; skipmissing=:pairwise);
print("pairwise(corkendall,eachcol(\$xm); skipmissing=:pairwise) ")
@btime res_3 = pairwise(corkendall, eachcol($xm); skipmissing=:pairwise);
print("corkendall(\$xm; skipmissing=:listwise)                   ")
@btime res_4 = corkendall($xm; skipmissing=:listwise);
print("pairwise(corkendall,eachcol(\$xm); skipmissing=:listwise) ")
@btime res_5 = pairwise(corkendall, eachcol($xm); skipmissing=:listwise);
print("corkendall(\$xm; skipmissing=:none)                       ")
@btime res_6 = corkendall($xm; skipmissing=:none);
print("pairwise(corkendall,eachcol(\$xm),skipmissing=:none)      ")
@btime res_7 = pairwise(corkendall, eachcol($xm), skipmissing=:none);
print("corspearman(\$x)                                          ")
@btime res_8 = corspearman($x);
print("corspearman(\$xm; skipmissing=:pairwise)                  ")
@btime res_9 = corspearman($xm; skipmissing=:pairwise);
print("pairwise(corspearman,eachcol(\$xm); skipmissing=:pairwise)")
@btime res_10 = pairwise(corspearman, eachcol($xm); skipmissing=:pairwise);
print("corspearman(\$xm; skipmissing=:listwise)                  ")
@btime res_11 = corspearman($xm; skipmissing=:listwise);
print("pairwise(corspearman,eachcol(\$xm); skipmissing=:listwise)")
@btime res_12 = pairwise(corspearman, eachcol($xm); skipmissing=:listwise);
print("corspearman(\$xm; skipmissing=:none)                      ")
@btime res_13 = corspearman($xm; skipmissing=:none);
print("pairwise(corspearman,eachcol(\$xm),skipmissing=:none)     ")
@btime res_14 = pairwise(corspearman, eachcol($xm), skipmissing=:none);

println("="^130)

#=
==================================================================================================================================
Dates.now() = DateTime("2024-03-30T10:16:29.303")
ENV["COMPUTERNAME"] = "PHILIP-LAPTOP"
Julia Version 1.10.2
Threads.nthreads() = 8
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
corkendall($x)                                             5.758 s (1193 allocations: 15.84 MiB)
corkendall($xm; skipmissing=:pairwise)                     5.250 s (1193 allocations: 15.51 MiB)
pairwise(corkendall,eachcol($xm); skipmissing=:pairwise)   5.541 s (1161 allocations: 15.51 MiB)
corkendall($xm; skipmissing=:listwise)                     9.694 ms (2194 allocations: 11.82 MiB)
pairwise(corkendall,eachcol($xm); skipmissing=:listwise)   9.863 ms (2162 allocations: 11.82 MiB)
corkendall($xm; skipmissing=:none)                         258.379 ms (1204 allocations: 16.47 MiB)
pairwise(corkendall,eachcol($xm),skipmissing=:none)        238.871 ms (1202 allocations: 16.47 MiB)
corspearman($x)                                            35.905 ms (1127 allocations: 39.31 MiB)
corspearman($xm; skipmissing=:pairwise)                    1.971 s (1260 allocations: 23.19 MiB)
pairwise(corspearman,eachcol($xm); skipmissing=:pairwise)  1.706 s (1228 allocations: 23.18 MiB)
corspearman($xm; skipmissing=:listwise)                    4.466 ms (2113 allocations: 19.43 MiB)
pairwise(corspearman,eachcol($xm); skipmissing=:listwise)  4.432 ms (2081 allocations: 19.43 MiB)
corspearman($xm; skipmissing=:none)                        20.017 ms (1625 allocations: 17.27 MiB)
pairwise(corspearman,eachcol($xm),skipmissing=:none)       20.085 ms (1623 allocations: 17.27 MiB)
==================================================================================================================================

=#