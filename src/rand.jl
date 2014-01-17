# Internal facilities for fast random number generation

function _randu(Ku::Uint, U::Uint)   # ~ U[0:Ku-1]
	x = rand(Uint)
	while x > U
		x = rand(Uint)
	end
	rem(x, Ku)
end

immutable RandIntSampler
	a::Int
	Ku::Uint
	U::Uint

	function RandIntSampler(K::Int) # 1:K
		Ku = uint(K)
		U = div(typemax(Uint), Ku) * Ku
		new(1, Ku, U)
	end

	function RandIntSampler(a::Int, b::Int)  # a:b
		Ku = uint(b - a + 1)
		U = div(typemax(Uint), Ku) * Ku
		new(a, Ku, U)
	end
end

rand(s::RandIntSampler) = int(_randu(s.Ku, s.U)) + s.a

randi(K::Int) = rand(RandIntSampler(K))
randi(a::Int, b::Int) = rand(RandIntSampler(a, b))


# draw a number from a binomial distribution

rand_binom(n::Real, p::Real) = int(ccall((:rbinom, :libRmath), Float64, (Float64, Float64), n, p))
