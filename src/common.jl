# common utilities

function depcheck(fname::Symbol, varname::Symbol, b::Union{Bool,Nothing})
    if b === nothing
        msg = "$fname will default to $varname=true in the future. Use $varname=false for previous behaviour."
        Base.depwarn(msg, fname)
        false
    else
        b
    end
end

_add((x1, x2)::Tuple{<:Real,<:Real}, (y1, y2)::Tuple{<:Real,<:Real}) = (x1 + y1, x2 + y2)
