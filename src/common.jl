# common utilities

function depcheck(fname::Symbol, varname::Symbol, b::Union{Bool, Nothing})
    if b === nothing
        msg = "$fname will default to $varname=true in the future. Use $varname=false for previous behaviour."
        Base.depwarn(msg, fname)
        false
    else
        b
    end
end
