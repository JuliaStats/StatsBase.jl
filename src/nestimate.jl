
"""
    nestimate(x::AbstractArray; method = "mvue")

Estimate ``N`` using an ``n`` sample of integer IDs. If the method is mvue 
(Minimum-variance unbiased estimator) then 

```math
m + \\frac{m}{n} - 1
```

is the estimate of ``N`` where ``n`` is the sample size and  ``m`` is the sample maximum. 
If the method is selected as "mom" then the the estimate of ``N``is 

```math
2 * \\bar{x} - 1
```

where ``\\bar{x}`` is the sample mean.

* `x::AbstractArray`: Vector of samples. 
* `method::String`: Either "mvue" or "mom". Default is "mvue".

# Description 

The problem is well-known as the German Tank problem, which is basically a special 
example of the problem of estimating ``N`` using a sample. 

# References: 
- https://en.wikipedia.org/wiki/German_tank_problem

"""
function nestimate(sampl::AbstractArray; method = "mvue")::Float64
    
    function mvue(sampl::AbstractArray)::Float64
        mx = maximum(sampl)
        estimate = mx + mx / length(sampl) - 1.0
        return estimate
    end 

    function mom(sampl::AbstractArray)::Float64
        return 2.0 * mean(sampl) - 1.0
    end 
   
    if method == "mvue"
        return mvue(sampl)
    elseif method == "mom"
        return mom(sampl)
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end
