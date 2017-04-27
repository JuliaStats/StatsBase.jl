using StatsBase
using Base.Test

##### weighted var & std

x = rand(10)
wv = fweights(rand(10))
m = mean(x, wv)

@test var(x, wv; corrected=false)           ≈ sum(abs2.(x .- m), wv) ./ sum(wv)
@test var(x, wv; mean=0, corrected=false)   ≈ sum(abs2.(x), wv) ./ sum(wv)
@test var(x, wv; mean=1.0, corrected=false) ≈ sum(abs2.(x .- 1.0), wv) ./ sum(wv)

@test std(x, wv; corrected=false)           ≈ sqrt(var(x, wv; corrected=false))
@test std(x, wv; mean=0, corrected=false)   ≈ sqrt(var(x, wv; mean=0, corrected=false))
@test std(x, wv; mean=1.0, corrected=false) ≈ sqrt(var(x, wv; mean=1.0, corrected=false))

(m, v) = mean_and_var(x; corrected=false)
@test m == mean(x)
@test v == var(x; corrected=false)

(m, s) = mean_and_std(x; corrected=false)
@test m == mean(x)
@test s == std(x; corrected=false)

(m, v) = mean_and_var(x, wv; corrected=false)
@test m == mean(x, wv)
@test v == var(x, wv; corrected=false)

(m, s) = mean_and_std(x, wv; corrected=false)
@test m == mean(x, wv)
@test s == std(x, wv; corrected=false)

x = rand(5, 6)
w1 = rand(5)
w2 = rand(6)
wv1 = fweights(w1)
wv2 = fweights(w2)
m1 = mean(x, wv1, 1)
m2 = mean(x, wv2, 2)

@test var(x, wv1, 1; mean=0, corrected=false) ≈ sum(abs2.(x) .* w1, 1) ./ sum(wv1)
@test var(x, wv2, 2; mean=0, corrected=false) ≈ sum(abs2.(x) .* w2', 2) ./ sum(wv2)

@test var(x, wv1, 1; mean=m1, corrected=false) ≈ sum(abs2.(x .- m1) .* w1, 1) ./ sum(wv1)
@test var(x, wv2, 2; mean=m2, corrected=false) ≈ sum(abs2.(x .- m2) .* w2', 2) ./ sum(wv2)

@test var(x, wv1, 1; corrected=false) ≈ sum(abs2.(x .- m1) .* w1, 1) ./ sum(wv1)
@test var(x, wv2, 2; corrected=false) ≈ sum(abs2.(x .- m2) .* w2', 2) ./ sum(wv2)

@test std(x, wv1, 1; corrected=false)          ≈ sqrt.(var(x, wv1, 1; corrected=false))
@test std(x, wv2, 2; corrected=false)          ≈ sqrt.(var(x, wv2, 2; corrected=false))
@test std(x, wv1, 1; mean=0, corrected=false)  ≈ sqrt.(var(x, wv1, 1; mean=0, corrected=false))
@test std(x, wv2, 2; mean=0, corrected=false)  ≈ sqrt.(var(x, wv2, 2; mean=0, corrected=false))
@test std(x, wv1, 1; mean=m1, corrected=false) ≈ sqrt.(var(x, wv1, 1; mean=m1, corrected=false))
@test std(x, wv2, 2; mean=m2, corrected=false) ≈ sqrt.(var(x, wv2, 2; mean=m2, corrected=false))

for d in 1:2
    (m, v) = mean_and_var(x, d; corrected=false)
    @test m == mean(x, d)
    @test v == var(x, d; corrected=false)

    (m, s) = mean_and_std(x, d; corrected=false)
    @test m == mean(x, d)
    @test s == std(x, d; corrected=false)
end

(m, v) = mean_and_var(x, wv1, 1; corrected=false)
@test m == mean(x, wv1, 1)
@test v == var(x, wv1, 1; corrected=false)

(m, v) = mean_and_var(x, wv2, 2; corrected=false)
@test m == mean(x, wv2, 2)
@test v == var(x, wv2, 2; corrected=false)

(m, s) = mean_and_std(x, wv1, 1; corrected=false)
@test m == mean(x, wv1, 1)
@test s == std(x, wv1, 1; corrected=false)

(m, s) = mean_and_std(x, wv2, 2; corrected=false)
@test m == mean(x, wv2, 2)
@test s == std(x, wv2, 2; corrected=false)



##### skewness & kurtosis

wv = fweights(ones(5) * 2.0)

@test skewness(1:5; corrected=false)             ≈  0.0
@test skewness([1, 2, 3, 4, 5]; corrected=false) ≈  0.0
@test skewness([1, 2, 2, 2, 5]; corrected=false) ≈  1.1731251294063556
@test skewness([1, 4, 4, 4, 5]; corrected=false) ≈ -1.1731251294063556

@test skewness([1, 2, 2, 2, 5], wv; corrected=false) ≈ 1.1731251294063556

@test kurtosis(1:5; corrected=false)             ≈ -1.3
@test kurtosis([1, 2, 3, 4, 5]; corrected=false) ≈ -1.3
@test kurtosis([1, 2, 3, 3, 2]; corrected=false) ≈ -1.1530612244897953

@test kurtosis([1, 2, 3, 4, 5], wv; corrected=false) ≈ -1.3


##### general moments

x = collect(2.0:8.0)
@test moment(x, 2; corrected=false) ≈ sum((x .- 5).^2) / length(x)
@test moment(x, 3; corrected=false) ≈ sum((x .- 5).^3) / length(x)
@test moment(x, 4; corrected=false) ≈ sum((x .- 5).^4) / length(x)
@test moment(x, 5; corrected=false) ≈ sum((x .- 5).^5) / length(x)

@test moment(x, 2, 4.0; corrected=false) ≈ sum((x .- 4).^2) / length(x)
@test moment(x, 3, 4.0; corrected=false) ≈ sum((x .- 4).^3) / length(x)
@test moment(x, 4, 4.0; corrected=false) ≈ sum((x .- 4).^4) / length(x)
@test moment(x, 5, 4.0; corrected=false) ≈ sum((x .- 4).^5) / length(x)

w = fweights([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
x2 = collect(2.0:6.0)
@test moment(x, 2, w; corrected=false) ≈ sum((x2 .- 4).^2) / 5
@test moment(x, 3, w; corrected=false) ≈ sum((x2 .- 4).^3) / 5
@test moment(x, 4, w; corrected=false) ≈ sum((x2 .- 4).^4) / 5
@test moment(x, 5, w; corrected=false) ≈ sum((x2 .- 4).^5) / 5

# Test corrected cases (this will be cleaner in testsets)
x = rand(10)

# AnalyticWeights
@test var(x, aweights(ones(10))) ≈ var(x)

w = aweights(rand(10))
n = length(w) # Could be count(!iszero, w) instead
w = aweights(w .* (n / sum(w)))
sw = sum(w) # This is now equal to n, but maybe we should support non-normalized weights?
xbar = sum(w .* x) ./ sw
expected = sum(w .* (x .- xbar).^2)/(sw - sum(w.^2)/sw)
@test var(x, w) ≈ expected

# FrequencyWeights
@test var(x, fweights(ones(Int, 10))) ≈ var(x)
w = fweights(rand(UInt, 10))
sw = sum(w)
xbar = sum(w .* x) / sw
expected = sum(w .* (x .- xbar).^2) ./ (sum(w) - 1)
@test var(x, w) ≈ expected

# ProbabilityWeights
@test var(x, pweights(ones(10))) ≈ var(x)
w = pweights(rand(10))
n = count(!iszero, w)
sw = sum(w)
xbar = sum(w .* x)/sw
expected = sum(w .* (x .- xbar).^2)/sw * n/(n - 1)
@test var(x, w) ≈ expected
