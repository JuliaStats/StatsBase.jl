using Base: @deprecate, @deprecate_binding, depwarn

if !isdefined(Base, :stderr)
    @deprecate stderr(obj::StatisticalModel) stderror(obj)
else
    function (io::typeof(stderr))(obj::StatisticalModel)
        Base.depwarn("stderr(obj::StatisticalModel) is deprecated, use stderror(obj) instead", :stderr)
        io === stderr ? stderror(obj) : throw(MethodErrror(io, (obj,)))
    end
end

@deprecate model_response(obj::StatisticalModel) response(obj)

@deprecate norepeats(a::AbstractArray) allunique(a)

@deprecate(mad!(v::AbstractArray{<:Real}, center;
                constant::Real = BigFloat("1.482602218505601860547076529360423431326703202590312896536266275245674447622701")),
           mad!(v, center=center, constant=constant))
