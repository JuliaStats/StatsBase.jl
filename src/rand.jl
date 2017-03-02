# Internal facilities for fast random number generation

immutable RandIntSampler  # for generating Int samples in [0, K-1]
    a::Int
    Ku::UInt
    U::UInt

    @compat RandIntSampler(K::Int) = (Ku = UInt(K); new(1, Ku, div(typemax(UInt), Ku) * Ku))
    @compat RandIntSampler(a::Int, b::Int) = (Ku = UInt(b-a+1); new(a, Ku, div(typemax(UInt), Ku) * Ku))
end

function rand(s::RandIntSampler, rng::AbstractRNG = Base.GLOBAL_RNG)
    x = rand(rng, UInt)
    while x >= s.U
        x = rand(rng, UInt)
    end
    @compat s.a + Int(rem(x, s.Ku))
end

randi(K::Int, rng::AbstractRNG = Base.GLOBAL_RNG) = rand(RandIntSampler(K), rng)
randi(a::Int, b::Int, rng::AbstractRNG = Base.GLOBAL_RNG) = rand(RandIntSampler(a, b), rng)
