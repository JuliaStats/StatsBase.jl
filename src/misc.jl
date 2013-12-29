# Miscelleneous stuff

# run-length encoding
function rle{T}(v::Vector{T})
    n = length(v)
    vals = T[]
    lens = Int[]

    cv = v[1]
    cl = 1

    i = 2
    @inbounds while i <= n
        vi = v[i]
        if vi == cv
            cl += 1
        else
            push!(vals, cv)
            push!(lens, cl)
            cv = vi
            cl = 1
        end
        i += 1
    end

    # the last section
    push!(vals, cv)
    push!(lens, cl)

    return (vals, lens)
end

# inverse run-length encoding
function inverse_rle{T}(vals::AbstractVector{T}, lens::IntegerVector)
    m = length(vals)
    length(lens) == m || raise_dimerror()

    r = Array(T, sum(lens))
    p = 0
    @inbounds for i = 1 : m
        j = lens[i]
        v = vals[i]
        while j > 0
            r[p+=1] = v
            j -=1 
        end
    end
    return r
end


# findat (get positions (within a) for elements in b)

function indexmap{T}(a::AbstractArray{T})
    d = (T=>Int)[]
    for i = 1 : length(a)
        @inbounds k = a[i]
        if !haskey(d, k)
            d[k] = i
        end
    end
    return d
end

function findat!{T}(r::IntegerArray, a::AbstractArray{T}, b::AbstractArray{T})
    length(r) == length(b) || raise_dimerror()
    d = indexmap(a)
    @inbounds for i = 1 : length(b)
        r[i] = get(d, b[i], 0)
    end
    return r
end

findat(a::AbstractArray, b::AbstractArray) = findat!(Array(Int, size(b)), a, b)



function indicators{T}(input::AbstractMatrix{T},
                       categories::Array{Any,1}={};
                       sparse::Bool=false)
    nfeatures, nsamples = size(input)
    if length(categories) != 0 && length(categories) != nfeatures
        error("You must provide either categories for each feature or no categories")
    end
    internal_categories = copy(categories)
    noutrows = 0
    if length(internal_categories) != nfeatures
        for i in 1:nfeatures
            push!(internal_categories, sort(unique(input[i, :])))
        end
    end
    for i in 1:nfeatures
        noutrows += length(internal_categories[i])
    end
    if sparse
        output = spzeros(noutrows, nsamples)
    else
        output = zeros(noutrows, nsamples)
    end
    offset = 1
    for i in 1:nfeatures
        indicators!(output, offset, slice(input, i, :), internal_categories[i])
        offset += length(internal_categories[i])
    end
    return output
end

function indicators{T}(input::AbstractVector{T},
                       categories::Array{T,1}=sort(unique(input));
                       sparse::Bool=false)
    if sparse
        output = spzeros(length(categories), length(input))
    else
        output = zeros(length(categories), length(input))
    end
    indicators!(output, 1, input, categories)
    return output
end
 
function indicators!{S<:Real,T}(output::AbstractArray{S},
                                offset::Integer,
                                input::AbstractVector{T},
                                categories::Array{T,1}=sort(unique(input)))
    indices = (T=>Integer)[categories[i]=>i for i in 1:length(categories)]
    const lo = offset-1
    for i in 1:length(input)
        output[indices[input[i]]+lo, i] = one(S)
    end
    return
end

