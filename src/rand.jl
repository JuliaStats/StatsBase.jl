# Internal facilities for fast random number generation

immutable RandIntSampler  # for generating Int samples in [0, K-1]
    a::Int
    Ku::UInt
    U::UInt

    @compat RandIntSampler(K::Int) = (Ku = UInt(K); new(1, Ku, div(typemax(UInt), Ku) * Ku))
    @compat RandIntSampler(a::Int, b::Int) = (Ku = UInt(b-a+1); new(a, Ku, div(typemax(UInt), Ku) * Ku))
end

function rand(rng::AbstractRNG, s::RandIntSampler)
    x = rand(rng, UInt)
    while x >= s.U
        x = rand(rng, UInt)
    end
    @compat s.a + Int(rem(x, s.Ku))
end
rand(s::RandIntSampler) = rand(Base.GLOBAL_RNG, s)

randi(rng::AbstractRNG, K::Int) = rand(rng, RandIntSampler(K))
randi(K::Int) = randi(Base.GLOBAL_RNG, K)
randi(rng::AbstractRNG, a::Int, b::Int) = rand(rng, RandIntSampler(a, b))
randi(a::Int, b::Int) = randi(Base.GLOBAL_RNG, a, b)
