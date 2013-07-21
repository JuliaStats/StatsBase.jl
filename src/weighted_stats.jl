# Weighted statistics

function wmean{T<:Number,W<:Real}(v::AbstractArray{T}, w::AbstractArray{W})
    sv = 0.0
    sw = 0.0
    for i in 1:length(v)
    	wi = w[i]
        sv += v[i] * wi
        sw += wi
    end
    return sv / sw
end


