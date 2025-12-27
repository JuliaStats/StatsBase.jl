using BenchmarkTools
using StatsBase
using Random

println("=" ^ 60)
println("Half-Sample Mode (hsm_mode) Benchmarks")
println("=" ^ 60)

# Small dataset
println("\n1. Small dataset (n=100)")
Random.seed!(123)
small = randn(100)
@btime hsm_mode($small)

# Medium dataset
println("\n2. Medium dataset (n=1,000)")
Random.seed!(123)
medium = randn(1_000)
@btime hsm_mode($medium)

# Large dataset
println("\n3. Large dataset (n=10,000)")
Random.seed!(123)
large = randn(10_000)
@btime hsm_mode($large)

# Very large dataset
println("\n4. Very large dataset (n=100,000)")
Random.seed!(123)
very_large = randn(100_000)
@btime hsm_mode($very_large)

# With outliers
println("\n5. Dataset with outliers (n=1,000)")
Random.seed!(123)
with_outliers = vcat(randn(990), fill(1000.0, 10))
@btime hsm_mode($with_outliers)

# Skewed distribution
println("\n6. Skewed distribution (n=1,000)")
Random.seed!(123)
skewed = abs.(randn(1_000)) .^ 2
@btime hsm_mode($skewed)

println("\n" * "=" ^ 60)
println("Benchmark complete!")
println("=" ^ 60)
