# Internal facilities for fast random number generation

immutable RandIntSampler  # for generating Int samples in [0, K-1]
	a::Int
	Ku::Uint
	U::Uint

	RandIntSampler(K::Int) = (Ku = uint(K); new(1, Ku, div(typemax(Uint), Ku) * Ku))
	RandIntSampler(a::Int, b::Int) = (Ku = uint(b-a+1); new(a, Ku, div(typemax(Uint), Ku) * Ku))
end

function rand(s::RandIntSampler)
	x = rand(Uint)
	while x >= s.U
		x = rand(Uint)
	end
	s.a + int(rem(x, s.Ku))
end

randi(K::Int) = rand(RandIntSampler(K))
randi(a::Int, b::Int) = rand(RandIntSampler(a, b))

# draw a number from a binomial distribution

rand_binom(n::Real, p::Real) = int(ccall((:rbinom, "libRmath-julia"), Float64, (Float64, Float64), n, p))
