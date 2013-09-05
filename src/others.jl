# Other stuff

# run-length encoding
function rle{T}(v::Vector{T})
    n = length(v)
    current_value = v[1]
    current_length = 1
    values = Array(T, n)
    total_values = 1
    lengths = Array(Int, n)
    total_lengths = 1
    for i in 2:n
        if v[i] == current_value
            current_length += 1
        else
            values[total_values] = current_value
            total_values += 1
            lengths[total_lengths] = current_length
            total_lengths += 1
            current_value = v[i]
            current_length = 1
        end
    end
    values[total_values] = current_value
    lengths[total_lengths] = current_length
    return (values[1:total_values], lengths[1:total_lengths])
end

# inverse run-length encoding
function inverse_rle{T}(values::Vector{T}, lengths::Vector{Int})
    total_n = sum(lengths)
    pos = 0
    res = Array(T, total_n)
    n = length(values)
    for i in 1:n
        v = values[i]
        l = lengths[i]
        for j in 1:l
            pos += 1
            res[pos] = v
        end
    end
    return res
end

# TODO: Support slicing along any dimensions
function findat!{T}(indices::Vector{Int},
                    a::AbstractArray,
                    b::AbstractArray{T})
    inds = Dict{T, Int}()
    for i in 1:length(b)
        tmp = b[i]
        if !haskey(inds, tmp)
            inds[tmp] = i
        end
    end
    for i = 1:length(a)
        indices[i] = get(inds, a[i], 0)
    end
    return
end

# TODO: Support slicing along any dimensions
function findat(a::AbstractArray, b::AbstractArray)
    indices = Array(Int, length(a))
    findat!(indices, a, b)
    return indices
end


## Empirical cummulative density function
function ecdf{T}(X::AbstractVector{T})
    Xs = sort(X)
    isnan(Xs[end]) && error("ecdf undefined in presence of NaNs")
    n = length(X)
    e(x::Real) = searchsortedlast(Xs, x) / n
    function e(v::Vector)
        ord = sortperm(v)
        m = length(v)
        r = Array(T, m)
        r0 = 0
        i = 1
        for x in Xs
            while i <= m && x > v[ord[i]]
                r[ord[i]] = r0
                i += 1
            end
            r0 += 1
            if i > m break end
        end
        while i <= m
            r[ord[i]] = n
            i += 1
        end
        return r / n
    end
    return e
end

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

abstract StatisticalModel

coef(obj::StatisticalModel) = error("No method defined")
coeftable(obj::StatisticalModel) = error("No method defined")
confint(obj::StatisticalModel) = error("No method defined")
loglikelihood(obj::StatisticalModel) = error("No method defined")
nobs(obj::StatisticalModel) = size(model_response(obj), 1)
stderr(obj::StatisticalModel) = sqrt(diag(vcov(obj)))
vcov(obj::StatisticalModel) = error("No method defined")

abstract RegressionModel <: StatisticalModel

residuals(obj::RegressionModel) = error("No method defined")
model_response(obj::RegressionModel) = error("No method defined")
predict(obj::RegressionModel) = error("No method defined")
