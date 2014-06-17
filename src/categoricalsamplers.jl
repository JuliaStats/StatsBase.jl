
##### Alias Table #####

immutable AliasTable <: Sampler{Univariate,Discrete}
    accept::Vector{Float64}
    alias::Vector{Int}
    isampler::RandIntSampler
end

function AliasTable(probs)

    n = length(probs)
    accept = float64(probs*n)

    alias = Array(Int,n)
    larges = Array(Int,0)
    smalls = Array(Int,0)

    for i = 1:n
        acci = accept[i] 
        if acci > 1.0 
            push!(larges,i)
        elseif acci < 1.0
            push!(smalls,i)
        end
    end
    while !isempty(larges) && !isempty(smalls)
        s = pop!(smalls)
        l = pop!(larges)
        alias[s] = l
        accept[l] = (accept[l] - 1.0) + accept[s]
        if accept[l] > 1
            push!(larges,l)
        else
            push!(smalls,l)
        end
    end

    # this loop should be redundant, except for rounding
    for s = smalls
        accept[s] = 1.0
    end

    AliasTable(accept, alias, RandIntSampler(n))
end

function rand(a::AliasTable)
    i = rand(a.isampler)
    u = rand()
    u < a.accept[i] ? i : a.alias[i]
end

Base.show(io::IO, a::AliasTable) = @printf io "AliasTable with %d entries" length(a.accept)



##### Draw Table #####

# Store an alias table
immutable DiscreteDistributionTable <: Sampler{Univariate,Discrete}
    table::Vector{Vector{Int64}}
    bounds::Vector{Int64}
end

# TODO: Test if bit operations can speed up Base64 mod's and fld's
function DiscreteDistributionTable{T <: Real}(probs::Vector{T})
    # Cache the cardinality of the outcome set
    n = length(probs)

    # Convert all Float64's into integers
    vals = Array(Int64, n)
    for i in 1:n
        vals[i] = int(probs[i] * 64^9)
    end

    # Allocate digit table and digit sums as table bounds
    table = Array(Vector{Int64}, 9)
    bounds = zeros(Int64, 9)

    # Special case for deterministic distributions
    for i in 1:n
        if vals[i] == 64^9
            table[1] = Array(Int64, 64)
            for j in 1:64
                table[1][j] = i
            end
            bounds[1] = 64^9
            for j in 2:9
                table[j] = Array(Int64, 0)
                bounds[j] = 64^9
            end
            return DiscreteDistributionTable(table, bounds)
        end
    end

    # Fill tables
    multiplier = 1
    for index in 9:-1:1
        counts = Array(Int64, 0)
        for i in 1:n
            digit = mod(vals[i], 64)
            # vals[i] = fld(vals[i], 64)
            vals[i] >>= 6
            bounds[index] += digit
            for itr in 1:digit
                push!(counts, i)
            end
        end
        bounds[index] *= multiplier
        table[index] = counts
        # multiplier *= 64
        multiplier <<= 6
    end

    # Make bounds cumulative
    bounds = cumsum(bounds)

    return DiscreteDistributionTable(table, bounds)
end

function rand(table::DiscreteDistributionTable)
    # 64^9 - 1 == 0x003fffffffffffff
    i = rand(1:(64^9 - 1))
    # if i == 64^9
    #   return table.table[9][rand(1:length(table.table[9]))]
    # end
    bound = 1
    while i > table.bounds[bound] && bound < 9
        bound += 1
    end
    if bound > 1
        index = fld(i - table.bounds[bound - 1] - 1, 64^(9 - bound)) + 1
    else
        index = fld(i - 1, 64^(9 - bound)) + 1
    end
    return table.table[bound][index]
end

Base.show(io::IO, table::DiscreteDistributionTable) = @printf io "DiscreteDistributionTable"


##### Huffman Table ######

abstract HuffmanNode{T} <: Sampler{Univariate,Discrete}

immutable HuffmanLeaf{T} <: HuffmanNode{T}
    value::T
    weight::Uint64
end

immutable HuffmanBranch{T} <: HuffmanNode{T}
    left::HuffmanNode{T}
    right::HuffmanNode{T}
    weight::Uint64
end
HuffmanBranch{T}(ha::HuffmanNode{T},hb::HuffmanNode{T}) = HuffmanBranch(ha, hb, ha.weight + hb.weight)

Base.isless{T}(ha::HuffmanNode{T}, hb::HuffmanNode{T}) = isless(ha.weight,hb.weight)
Base.show{T}(io::IO, t::HuffmanNode{T}) = show(io,typeof(t))

function Base.getindex{T}(h::HuffmanBranch{T},u::Uint64) 
    while isa(h,HuffmanBranch{T})
        if u < h.left.weight
            h = h.left
        else 
            u -= h.left.weight
            h = h.right
        end
    end
    h.value
end

# build the huffman tree
# could be slightly more efficient using a Deque.
function huffman{T}(values::AbstractVector{T},weights::AbstractVector{Uint64})
    leafs = [HuffmanLeaf{T}(values[i],weights[i]) for i = 1:length(weights)]
    sort!(leafs; rev=true)
    
    branches = Array(HuffmanBranch{T},0)
        
    while !isempty(leafs) || length(branches) > 1
        left = isempty(branches) || (!isempty(leafs) && first(leafs) < first(branches)) ? pop!(leafs) : pop!(branches)
        right = isempty(branches) || (!isempty(leafs) && first(leafs) < first(branches)) ? pop!(leafs) : pop!(branches)
        unshift!(branches,HuffmanBranch(left,right))
    end
    
    pop!(branches)    
end

function rand{T}(h::HuffmanNode{T})
    w = h.weight
    # generate uniform Uint64 objects on the range 0:(w-1)
    # unfortunately we can't use Range objects, as they don't have sufficient length
    u = rand(Uint64)
    if (w & (w-1)) == 0
        # power of 2
        u = u & (w-1)
    else
        m = typemax(Uint64)
        lim = m - (rem(m,w)+1)
        while u > lim
            u = rand(Uint64)
        end
        u = rem(u,w)
    end
    h[u]
end


