# Internal facilities for fast random number generation

immutable RandIntSampler  # for generating Int samples in [0, K-1]
    a::Int
    Ku::UInt
    U::UInt

    @compat RandIntSampler(K::Int) = (Ku = UInt(K); new(1, Ku, div(typemax(UInt), Ku) * Ku))
    @compat RandIntSampler(a::Int, b::Int) = (Ku = UInt(b-a+1); new(a, Ku, div(typemax(UInt), Ku) * Ku))
end

function rand(s::RandIntSampler)
    x = rand(UInt)
    while x >= s.U
        x = rand(UInt)
    end
    @compat s.a + Int(rem(x, s.Ku))
end

randi(K::Int) = rand(RandIntSampler(K))
randi(a::Int, b::Int) = rand(RandIntSampler(a, b))

# draw a number from a binomial distribution

rand_binom(n::Real, p::Real) =
    @compat Int(ccall((:rbinom, "libRmath-julia"), Float64, (Float64, Float64), n, p))
