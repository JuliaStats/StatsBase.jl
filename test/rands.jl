# test of random shuffle and random sample

using Stats
using Base.Test

srand(1234)

function verify_randsample(r::Matrix{Int}, m::Int, tol::Float64)
    # r is a matrix of size [q, n]
    # where each column of r is random sample without replacement from 1:m
    # this function verify this by looking at statistics
    
    q = size(r, 1)
    n = size(r, 2)
    P = zeros(q, m)
    for i = 1 : q
        for j = 1 : m
            P[i,j] = nnz(r[i,:] .== j) / n
        end
    end
    
    P0 = fill(1/m, (q, m))
    max(abs(P[:] - P0[:])) < tol
end



# rand shuffle

m = 5
n = 10000

r = zeros(Int, (m, n))
for i = 1 : n
    s = [1:m]
    randshuffle!(s, m)
    @assert sort(s) == [1:m]
    r[:,i] = s
end

@test verify_randsample(r, m, 0.02)

# rand pick

r = zeros(Int, (1, n))
for i = 1 : n
    r[:,i] = randsample(1:m, 1)
end
@test verify_randsample(r, m, 0.02)

r = zeros(Int, (2, n))
for i = 1 : n
    r[:,i] = randsample(1:m, 2)
end
@test verify_randsample(r, m, 0.02)

r = zeros(Int, (4, n))
for i = 1 : n
    r[:,i] = randsample(1:m, 4)
end
@test verify_randsample(r, m, 0.02)

m = 10^2
r = zeros(Int, (4, n))
for i = 1 : n
    r[:,i] = randsample(1:m, 4)
end
@test verify_randsample(r, m, 0.02)

