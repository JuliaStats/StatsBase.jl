
##### Alias Table #####

immutable AliasTable <: Sampler{Univariate,Discrete}
    accept::Vector{Float64}
    alias::Vector{Int}
    isampler::RandIntSampler
end
numcategories(s::AliasTable) = length(s.accept)

function AliasTable{T<:Real}(probs::AbstractVector{T})
    n = length(probs)
    n > 0 || error("The input probability vector is empty.")

    accept = Array(Float64, n)
    for i=1:n
        @inbounds accept[i] = probs[i] * n
    end
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
    for s in smalls
        accept[s] = 1.0
    end

    AliasTable(accept, alias, RandIntSampler(n))
end

function rand(s::AliasTable)
    i = rand(s.isampler)
    u = rand()
    u < s.accept[i] ? i : s.alias[i]
end

Base.show(io::IO, s::AliasTable) = @printf(io, "AliasTable with %d entries", numcategories(s))
