using StatsBase
using Base.Test

##### weighted var & std

x = rand(10)
wv = weights(rand(10))
m = mean(x, wv)

@test_approx_eq var(x, wv) sum(abs2(x .- m), wv) ./ sum(wv)
@test_approx_eq var(x, wv; mean=0) sum(abs2(x), wv) ./ sum(wv)
@test_approx_eq var(x, wv; mean=1.0) sum(abs2(x .- 1.0), wv) ./ sum(wv)

@test_approx_eq std(x, wv) sqrt(var(x, wv))
@test_approx_eq std(x, wv; mean=0) sqrt(var(x, wv; mean=0))
@test_approx_eq std(x, wv; mean=1.0) sqrt(var(x, wv; mean=1.0))

(m, v) = mean_and_var(x)
@test m == mean(x)
@test v == var(x)

(m, s) = mean_and_std(x)
@test m == mean(x)
@test s == std(x)

(m, v) = mean_and_var(x, wv)
@test m == mean(x, wv)
@test v == var(x, wv)

(m, s) = mean_and_std(x, wv)
@test m == mean(x, wv)
@test s == std(x, wv)

x = rand(5, 6)
w1 = rand(5)
w2 = rand(6)
wv1 = weights(w1)
wv2 = weights(w2)
m1 = mean(x, wv1, 1)
m2 = mean(x, wv2, 2)

@test_approx_eq var(x, wv1, 1; mean=0) sum(abs2(x) .* w1, 1) ./ sum(wv1)
@test_approx_eq var(x, wv2, 2; mean=0) sum(abs2(x) .* w2', 2) ./ sum(wv2)

@test_approx_eq var(x, wv1, 1; mean=m1) sum(abs2(x .- m1) .* w1, 1) ./ sum(wv1)
@test_approx_eq var(x, wv2, 2; mean=m2) sum(abs2(x .- m2) .* w2', 2) ./ sum(wv2)

@test_approx_eq var(x, wv1, 1) sum(abs2(x .- m1) .* w1, 1) ./ sum(wv1)
@test_approx_eq var(x, wv2, 2) sum(abs2(x .- m2) .* w2', 2) ./ sum(wv2)

@test_approx_eq std(x, wv1, 1) sqrt(var(x, wv1, 1))
@test_approx_eq std(x, wv2, 2) sqrt(var(x, wv2, 2))
@test_approx_eq std(x, wv1, 1; mean=0) sqrt(var(x, wv1, 1; mean=0))
@test_approx_eq std(x, wv2, 2; mean=0) sqrt(var(x, wv2, 2; mean=0))
@test_approx_eq std(x, wv1, 1; mean=m1) sqrt(var(x, wv1, 1; mean=m1))
@test_approx_eq std(x, wv2, 2; mean=m2) sqrt(var(x, wv2, 2; mean=m2))

for d in 1:2
    (m, v) = mean_and_var(x, d)
    @test m == mean(x, d)
    @test v == var(x, d)

    (m, s) = mean_and_std(x, d)
    @test m == mean(x, d)
    @test s == std(x, d)
end

(m, v) = mean_and_var(x, wv1, 1)
@test m == mean(x, wv1, 1)
@test v == var(x, wv1, 1)

(m, v) = mean_and_var(x, wv2, 2)
@test m == mean(x, wv2, 2)
@test v == var(x, wv2, 2)

(m, s) = mean_and_std(x, wv1, 1)
@test m == mean(x, wv1, 1)
@test s == std(x, wv1, 1)

(m, s) = mean_and_std(x, wv2, 2)
@test m == mean(x, wv2, 2)
@test s == std(x, wv2, 2)



##### skewness & kurtosis

wv = weights(ones(5) * 2.0)

@test_approx_eq skewness(1:5) 0.0
@test_approx_eq skewness([1, 2, 3, 4, 5]) 0.0
@test_approx_eq skewness([1, 2, 2, 2, 5]) 1.1731251294063556
@test_approx_eq skewness([1, 4, 4, 4, 5]) -1.1731251294063556

@test_approx_eq skewness([1, 2, 2, 2, 5], wv) 1.1731251294063556

@test_approx_eq kurtosis(1:5) -1.3
@test_approx_eq kurtosis([1, 2, 3, 4, 5]) -1.3
@test_approx_eq kurtosis([1, 2, 3, 3, 2]) -1.1530612244897953

@test_approx_eq kurtosis([1, 2, 3, 4, 5], wv) -1.3


##### general moments

x = float64([2:8;])
@test_approx_eq moment(x, 2) sum((x .- 5).^2) / length(x)
@test_approx_eq moment(x, 3) sum((x .- 5).^3) / length(x)
@test_approx_eq moment(x, 4) sum((x .- 5).^4) / length(x)
@test_approx_eq moment(x, 5) sum((x .- 5).^5) / length(x)

@test_approx_eq moment(x, 2, 4.0) sum((x .- 4).^2) / length(x)
@test_approx_eq moment(x, 3, 4.0) sum((x .- 4).^3) / length(x)
@test_approx_eq moment(x, 4, 4.0) sum((x .- 4).^4) / length(x)
@test_approx_eq moment(x, 5, 4.0) sum((x .- 4).^5) / length(x)

w = weights([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
x2 = float64([2:6;])
@test_approx_eq moment(x, 2, w) sum((x2 .- 4).^2) / 5
@test_approx_eq moment(x, 3, w) sum((x2 .- 4).^3) / 5
@test_approx_eq moment(x, 4, w) sum((x2 .- 4).^4) / 5
@test_approx_eq moment(x, 5, w) sum((x2 .- 4).^5) / 5