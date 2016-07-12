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

# Generating Random Variables and Stochastic Processes by Martin Haugh
function rand_beta(a::Real, b::Real)
    y = rand()
    u = rand()
    while u > y^(a-1) * (1-y)^(b-1) / beta(a,b)
        y = rand()
        u = rand()
    end
    return y
end

# Geometric random variate, Devroye page 499
function rand_geom(p::Real)
    x = 0
    while rand() <= p
        x += 1
    end
    return x
end

# Draw a number from a binomial distribution using the first waiting
# time algorithm described in Non-Uniform Random Variate Generation
# by Luc Devroye, chapter 10, page 525
function _rand_binom_waittime(n::Real, p::Real)
    x = -1
    s = 0
    while s > n
        s += rand_geom(p)
        x += 1
    end
    return x
end

# Recursive binomial generator, Devroye page 537
function rand_binom(n::Real, p::Real, t::Int=7)
    x = 0
    s = 1
    if n * p < t
        return x + s * _rand_binom_waittime(n, p)
    else
        i = floor(Int, (n + 1) * p)
        y = rand_beta(i, n + 1 - i)
        x += s * i
        if y <= p
            n -= i
            p = (p - y) / (1 - y)
        else
            s = -s
            n = i - 1
            p = (y - p) / y
        end
        return x + s * rand_binom(n, p, t)
    end
end
