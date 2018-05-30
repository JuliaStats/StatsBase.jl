# This file contains code formerly part of Julia. License is MIT: https://julialang.org/license

# deprecations from base/deprecated.jl
@deprecate varm(A::AbstractArray, m::AbstractArray, dims; kwargs...) varm(A, m; kwargs..., dims=dims)
@deprecate var(A::AbstractArray, dims; kwargs...)                    var(A; kwargs..., dims=dims)
@deprecate std(A::AbstractArray, dims; kwargs...)                    std(A; kwargs..., dims=dims)
@deprecate cov(X::AbstractMatrix, dim::Int; kwargs...)               cov(X; kwargs..., dims=dim)
@deprecate cov(x::AbstractVecOrMat, y::AbstractVecOrMat, dim::Int; kwargs...) cov(x, y; kwargs..., dims=dim)
@deprecate cor(X::AbstractMatrix, dim::Int)                          cor(X, dims=dim)
@deprecate cor(x::AbstractVecOrMat, y::AbstractVecOrMat, dim::Int)   cor(x, y, dims=dim)

# PR #21709
@deprecate cov(x::AbstractVector, corrected::Bool) cov(x, corrected=corrected)
@deprecate cov(x::AbstractMatrix, vardim::Int, corrected::Bool) cov(x, dims=vardim, corrected=corrected)
@deprecate cov(X::AbstractVector, Y::AbstractVector, corrected::Bool) cov(X, Y, corrected=corrected)
@deprecate cov(X::AbstractVecOrMat, Y::AbstractVecOrMat, vardim::Int, corrected::Bool) cov(X, Y, dims=vardim, corrected=corrected)

##### variances #####

# faster computation of real(conj(x)*y)
realXcY(x::Real, y::Real) = x*y
realXcY(x::Complex, y::Complex) = real(x)*real(y) + imag(x)*imag(y)

var(iterable; corrected::Bool=true, mean=nothing) = _var(iterable, corrected, mean)

function _var(iterable, corrected::Bool, mean)
    y = iterate(iterable)
    if y === nothing
        throw(ArgumentError("variance of empty collection undefined: $(repr(iterable))"))
    end
    count = 1
    value, state = y
    y = iterate(iterable, state)
    if mean === nothing
        # Use Welford algorithm as seen in (among other places)
        # Knuth's TAOCP, Vol 2, page 232, 3rd edition.
        M = value / 1
        S = real(zero(M))
        while y !== nothing
            value, state = y
            y = iterate(iterable, state)
            count += 1
            new_M = M + (value - M) / count
            S = S + realXcY(value - M, value - new_M)
            M = new_M
        end
        return S / (count - Int(corrected))
    elseif isa(mean, Number) # mean provided
        # Cannot use a compensated version, e.g. the one from
        # "Updating Formulae and a Pairwise Algorithm for Computing Sample Variances."
        # by Chan, Golub, and LeVeque, Technical Report STAN-CS-79-773,
        # Department of Computer Science, Stanford University,
        # because user can provide mean value that is different to mean(iterable)
        sum2 = abs2(value - mean::Number)
        while y !== nothing
            value, state = y
            y = iterate(iterable, state)
            count += 1
            sum2 += abs2(value - mean)
        end
        return sum2 / (count - Int(corrected))
    else
        throw(ArgumentError("invalid value of mean, $(mean)::$(typeof(mean))"))
    end
end

centralizedabs2fun(m) = x -> abs2.(x - m)
centralize_sumabs2(A::AbstractArray, m) =
    mapreduce(centralizedabs2fun(m), +, A)
centralize_sumabs2(A::AbstractArray, m, ifirst::Int, ilast::Int) =
    Base.mapreduce_impl(centralizedabs2fun(m), +, A, ifirst, ilast)

function centralize_sumabs2!(R::AbstractArray{S}, A::AbstractArray, means::AbstractArray) where S
    # following the implementation of _mapreducedim! at base/reducedim.jl
    lsiz = Base.check_reducedims(R,A)
    isempty(R) || fill!(R, zero(S))
    isempty(A) && return R

    if Base.has_fast_linear_indexing(A) && lsiz > 16
        nslices = div(Base._length(A), lsiz)
        ibase = first(LinearIndices(A))-1
        for i = 1:nslices
            @inbounds R[i] = centralize_sumabs2(A, means[i], ibase+1, ibase+lsiz)
            ibase += lsiz
        end
        return R
    end
    indsAt, indsRt = Base.safe_tail(axes(A)), Base.safe_tail(axes(R)) # handle d=1 manually
    keep, Idefault = Broadcast.shapeindexer(indsRt)
    if Base.reducedim1(R, A)
        i1 = first(Base.indices1(R))
        @inbounds for IA in CartesianIndices(indsAt)
            IR = Broadcast.newindex(IA, keep, Idefault)
            r = R[i1,IR]
            m = means[i1,IR]
            @simd for i in axes(A, 1)
                r += abs2(A[i,IA] - m)
            end
            R[i1,IR] = r
        end
    else
        @inbounds for IA in CartesianIndices(indsAt)
            IR = Broadcast.newindex(IA, keep, Idefault)
            @simd for i in axes(A, 1)
                R[i,IR] += abs2(A[i,IA] - means[i,IR])
            end
        end
    end
    return R
end

function varm!(R::AbstractArray{S}, A::AbstractArray, m::AbstractArray; corrected::Bool=true) where S
    if isempty(A)
        fill!(R, convert(S, NaN))
    else
        rn = div(Base._length(A), Base._length(R)) - Int(corrected)
        centralize_sumabs2!(R, A, m)
        R .= R .* (1 // rn)
    end
    return R
end

"""
    varm(v, m; dims, corrected::Bool=true)

Compute the sample variance of a collection `v` with known mean(s) `m`,
optionally over the given dimensions. `m` may contain means for each dimension of
`v`. If `corrected` is `true`, then the sum is scaled with `n-1`,
whereas the sum is scaled with `n` if `corrected` is `false` where `n = length(x)`.

!!! note
    Julia does not ignore `NaN` values in the computation. Use the [`missing`](@ref) type
    to represent missing values, and the [`skipmissing`](@ref) function to omit them.
"""
varm(A::AbstractArray, m::AbstractArray; corrected::Bool=true, dims=:) = _varm(A, m, corrected, dims)

_varm(A::AbstractArray{T}, m, corrected::Bool, region) where {T} =
    varm!(Base.reducedim_init(t -> abs2(t)/2, +, A, region), A, m; corrected=corrected)

varm(A::AbstractArray, m; corrected::Bool=true) = _varm(A, m, corrected, :)

function _varm(A::AbstractArray{T}, m, corrected::Bool, ::Colon) where T
    n = Base._length(A)
    n == 0 && return typeof((abs2(zero(T)) + abs2(zero(T)))/2)(NaN)
    return centralize_sumabs2(A, m) / (n - Int(corrected))
end


"""
    var(v; dims, corrected::Bool=true, mean=nothing)

Compute the sample variance of a vector or array `v`, optionally along the given dimensions.
The algorithm will return an estimator of the generative distribution's variance
under the assumption that each entry of `v` is an IID drawn from that generative
distribution. This computation is equivalent to calculating `sum(abs2, v - mean(v)) /
(length(v) - 1)`. If `corrected` is `true`, then the sum is scaled with `n-1`,
whereas the sum is scaled with `n` if `corrected` is `false` where `n = length(x)`.
The mean `mean` over the region may be provided.

!!! note
    Julia does not ignore `NaN` values in the computation. Use the [`missing`](@ref) type
    to represent missing values, and the [`skipmissing`](@ref) function to omit them.
"""
var(A::AbstractArray; corrected::Bool=true, mean=nothing, dims=:) = _var(A, corrected, mean, dims)

_var(A::AbstractArray, corrected::Bool, mean, dims) =
    varm(A, coalesce(mean, Base.mean(A, dims=dims)); corrected=corrected, dims=dims)

_var(A::AbstractArray, corrected::Bool, mean, ::Colon) =
    real(varm(A, coalesce(mean, Base.mean(A)); corrected=corrected))

varm(iterable, m; corrected::Bool=true) = _var(iterable, corrected, m)

## variances over ranges

varm(v::AbstractRange, m::AbstractArray) = range_varm(v, m)
varm(v::AbstractRange, m) = range_varm(v, m)

function range_varm(v::AbstractRange, m)
    f  = first(v) - m
    s  = step(v)
    l  = length(v)
    vv = f^2 * l / (l - 1) + f * s * l + s^2 * l * (2 * l - 1) / 6
    if l == 0 || l == 1
        return typeof(vv)(NaN)
    end
    return vv
end

function var(v::AbstractRange)
    s  = step(v)
    l  = length(v)
    vv = abs2(s) * (l + 1) * l / 12
    if l == 0 || l == 1
        return typeof(vv)(NaN)
    end
    return vv
end


##### standard deviation #####

function sqrt!(A::AbstractArray)
    for i in eachindex(A)
        @inbounds A[i] = sqrt(A[i])
    end
    A
end

stdm(A::AbstractArray, m; corrected::Bool=true) =
    sqrt.(varm(A, m; corrected=corrected))

"""
    std(v; corrected::Bool=true, mean=nothing, dims)

Compute the sample standard deviation of a vector or array `v`, optionally along the given
dimensions. The algorithm returns an estimator of the generative distribution's standard
deviation under the assumption that each entry of `v` is an IID drawn from that generative
distribution. This computation is equivalent to calculating `sqrt(sum((v - mean(v)).^2) /
(length(v) - 1))`. A pre-computed `mean` may be provided. If `corrected` is `true`,
then the sum is scaled with `n-1`, whereas the sum is scaled with `n` if `corrected` is
`false` where `n = length(x)`.

!!! note
    Julia does not ignore `NaN` values in the computation. Use the [`missing`](@ref) type
    to represent missing values, and the [`skipmissing`](@ref) function to omit them.
"""
std(A::AbstractArray; corrected::Bool=true, mean=nothing, dims=:) = _std(A, corrected, mean, dims)

_std(A::AbstractArray, corrected::Bool, mean, dims) =
    sqrt.(var(A; corrected=corrected, mean=mean, dims=dims))

_std(A::AbstractArray, corrected::Bool, mean, ::Colon) =
    sqrt.(var(A; corrected=corrected, mean=mean))

_std(A::AbstractArray{<:AbstractFloat}, corrected::Bool, mean, dims) =
    sqrt!(var(A; corrected=corrected, mean=mean, dims=dims))

_std(A::AbstractArray{<:AbstractFloat}, corrected::Bool, mean, ::Colon) =
    sqrt.(var(A; corrected=corrected, mean=mean))

std(iterable; corrected::Bool=true, mean=nothing) =
    sqrt(var(iterable, corrected=corrected, mean=mean))

"""
    stdm(v, m; corrected::Bool=true)

Compute the sample standard deviation of a vector `v`
with known mean `m`. If `corrected` is `true`,
then the sum is scaled with `n-1`, whereas the sum is
scaled with `n` if `corrected` is `false` where `n = length(x)`.

!!! note
    Julia does not ignore `NaN` values in the computation. Use the [`missing`](@ref) type
    to represent missing values, and the [`skipmissing`](@ref) function to omit them.
"""
stdm(iterable, m; corrected::Bool=true) =
    std(iterable, corrected=corrected, mean=m)


###### covariance ######

# auxiliary functions

_conj(x::AbstractArray{<:Real}) = x
_conj(x::AbstractArray) = conj(x)

_getnobs(x::AbstractVector, vardim::Int) = Base._length(x)
_getnobs(x::AbstractMatrix, vardim::Int) = size(x, vardim)

function _getnobs(x::AbstractVecOrMat, y::AbstractVecOrMat, vardim::Int)
    n = _getnobs(x, vardim)
    _getnobs(y, vardim) == n || throw(DimensionMismatch("dimensions of x and y mismatch"))
    return n
end

_vmean(x::AbstractVector, vardim::Int) = mean(x)
_vmean(x::AbstractMatrix, vardim::Int) = mean(x, dims=vardim)

# core functions

unscaled_covzm(x::AbstractVector{<:Number})    = sum(abs2, x)
unscaled_covzm(x::AbstractVector)              = sum(t -> t*t', x)
unscaled_covzm(x::AbstractMatrix, vardim::Int) = (vardim == 1 ? _conj(x'x) : x * x')

unscaled_covzm(x::AbstractVector, y::AbstractVector) = dot(y, x)
unscaled_covzm(x::AbstractVector, y::AbstractMatrix, vardim::Int) =
    (vardim == 1 ? *(transpose(x), _conj(y)) : *(transpose(x), transpose(_conj(y))))
unscaled_covzm(x::AbstractMatrix, y::AbstractVector, vardim::Int) =
    (c = vardim == 1 ? *(transpose(x), _conj(y)) :  x * _conj(y); reshape(c, length(c), 1))
unscaled_covzm(x::AbstractMatrix, y::AbstractMatrix, vardim::Int) =
    (vardim == 1 ? *(transpose(x), _conj(y)) : *(x, adjoint(y)))

# covzm (with centered data)

covzm(x::AbstractVector; corrected::Bool=true) = unscaled_covzm(x) / (Base._length(x) - Int(corrected))
function covzm(x::AbstractMatrix, vardim::Int=1; corrected::Bool=true)
    C = unscaled_covzm(x, vardim)
    T = promote_type(typeof(first(C) / 1), eltype(C))
    A = convert(AbstractMatrix{T}, C)
    b = 1//(size(x, vardim) - corrected)
    A .= A .* b
    return A
end
covzm(x::AbstractVector, y::AbstractVector; corrected::Bool=true) =
    unscaled_covzm(x, y) / (Base._length(x) - Int(corrected))
function covzm(x::AbstractVecOrMat, y::AbstractVecOrMat, vardim::Int=1; corrected::Bool=true)
    C = unscaled_covzm(x, y, vardim)
    T = promote_type(typeof(first(C) / 1), eltype(C))
    A = convert(AbstractArray{T}, C)
    b = 1//(_getnobs(x, y, vardim) - corrected)
    A .= A .* b
    return A
end

# covm (with provided mean)
## Use map(t -> t - xmean, x) instead of x .- xmean to allow for Vector{Vector}
## which can't be handled by broadcast
covm(x::AbstractVector, xmean; corrected::Bool=true) =
    covzm(map(t -> t - xmean, x); corrected=corrected)
covm(x::AbstractMatrix, xmean, vardim::Int=1; corrected::Bool=true) =
    covzm(x .- xmean, vardim; corrected=corrected)
covm(x::AbstractVector, xmean, y::AbstractVector, ymean; corrected::Bool=true) =
    covzm(map(t -> t - xmean, x), map(t -> t - ymean, y); corrected=corrected)
covm(x::AbstractVecOrMat, xmean, y::AbstractVecOrMat, ymean, vardim::Int=1; corrected::Bool=true) =
    covzm(x .- xmean, y .- ymean, vardim; corrected=corrected)

# cov (API)
"""
    cov(x::AbstractVector; corrected::Bool=true)

Compute the variance of the vector `x`. If `corrected` is `true` (the default) then the sum
is scaled with `n-1`, whereas the sum is scaled with `n` if `corrected` is `false` where `n = length(x)`.
"""
cov(x::AbstractVector; corrected::Bool=true) = covm(x, Base.mean(x); corrected=corrected)

"""
    cov(X::AbstractMatrix; dims::Int=1, corrected::Bool=true)

Compute the covariance matrix of the matrix `X` along the dimension `dims`. If `corrected`
is `true` (the default) then the sum is scaled with `n-1`, whereas the sum is scaled with `n`
if `corrected` is `false` where `n = size(X, dims)`.
"""
cov(X::AbstractMatrix; dims::Int=1, corrected::Bool=true) =
    covm(X, _vmean(X, dims), dims; corrected=corrected)

"""
    cov(x::AbstractVector, y::AbstractVector; corrected::Bool=true)

Compute the covariance between the vectors `x` and `y`. If `corrected` is `true` (the
default), computes ``\\frac{1}{n-1}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where
``*`` denotes the complex conjugate and `n = length(x) = length(y)`. If `corrected` is
`false`, computes ``\\frac{1}{n}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*``.
"""
cov(x::AbstractVector, y::AbstractVector; corrected::Bool=true) =
    covm(x, Base.mean(x), y, Base.mean(y); corrected=corrected)

"""
    cov(X::AbstractVecOrMat, Y::AbstractVecOrMat; dims::Int=1, corrected::Bool=true)

Compute the covariance between the vectors or matrices `X` and `Y` along the dimension
`dims`. If `corrected` is `true` (the default) then the sum is scaled with `n-1`, whereas
the sum is scaled with `n` if `corrected` is `false` where `n = size(X, dims) = size(Y, dims)`.
"""
cov(X::AbstractVecOrMat, Y::AbstractVecOrMat; dims::Int=1, corrected::Bool=true) =
    covm(X, _vmean(X, dims), Y, _vmean(Y, dims), dims; corrected=corrected)

##### correlation #####

"""
    clampcor(x)

Clamp a real correlation to between -1 and 1, leaving complex correlations unchanged
"""
clampcor(x::Real) = clamp(x, -1, 1)
clampcor(x) = x

# cov2cor!

function cov2cor!(C::AbstractMatrix{T}, xsd::AbstractArray) where T
    nx = length(xsd)
    size(C) == (nx, nx) || throw(DimensionMismatch("inconsistent dimensions"))
    for j = 1:nx
        for i = 1:j-1
            C[i,j] = adjoint(C[j,i])
        end
        C[j,j] = oneunit(T)
        for i = j+1:nx
            C[i,j] = clampcor(C[i,j] / (xsd[i] * xsd[j]))
        end
    end
    return C
end
function cov2cor!(C::AbstractMatrix, xsd, ysd::AbstractArray)
    nx, ny = size(C)
    length(ysd) == ny || throw(DimensionMismatch("inconsistent dimensions"))
    for (j, y) in enumerate(ysd)   # fixme (iter): here and in all `cov2cor!` we assume that `C` is efficiently indexed by integers
        for i in 1:nx
            C[i,j] = clampcor(C[i, j] / (xsd * y))
        end
    end
    return C
end
function cov2cor!(C::AbstractMatrix, xsd::AbstractArray, ysd)
    nx, ny = size(C)
    length(xsd) == nx || throw(DimensionMismatch("inconsistent dimensions"))
    for j in 1:ny
        for (i, x) in enumerate(xsd)
            C[i,j] = clampcor(C[i,j] / (x * ysd))
        end
    end
    return C
end
function cov2cor!(C::AbstractMatrix, xsd::AbstractArray, ysd::AbstractArray)
    nx, ny = size(C)
    (length(xsd) == nx && length(ysd) == ny) ||
        throw(DimensionMismatch("inconsistent dimensions"))
    for (i, x) in enumerate(xsd)
        for (j, y) in enumerate(ysd)
            C[i,j] = clampcor(C[i,j] / (x * y))
        end
    end
    return C
end

# corzm (non-exported, with centered data)

corzm(x::AbstractVector{T}) where {T} = one(real(T))
function corzm(x::AbstractMatrix, vardim::Int=1)
    c = unscaled_covzm(x, vardim)
    return cov2cor!(c, collect(sqrt(c[i,i]) for i in 1:min(size(c)...)))
end
corzm(x::AbstractVector, y::AbstractMatrix, vardim::Int=1) =
    cov2cor!(unscaled_covzm(x, y, vardim), sqrt(sum(abs2, x)), sqrt!(sum(abs2, y, dims=vardim)))
corzm(x::AbstractMatrix, y::AbstractVector, vardim::Int=1) =
    cov2cor!(unscaled_covzm(x, y, vardim), sqrt!(sum(abs2, x, dims=vardim)), sqrt(sum(abs2, y)))
corzm(x::AbstractMatrix, y::AbstractMatrix, vardim::Int=1) =
    cov2cor!(unscaled_covzm(x, y, vardim), sqrt!(sum(abs2, x, dims=vardim)), sqrt!(sum(abs2, y, dims=vardim)))

# corm

corm(x::AbstractVector{T}, xmean) where {T} = one(real(T))
corm(x::AbstractMatrix, xmean, vardim::Int=1) = corzm(x .- xmean, vardim)
function corm(x::AbstractVector, mx, y::AbstractVector, my)
    n = length(x)
    length(y) == n || throw(DimensionMismatch("inconsistent lengths"))
    n > 0 || throw(ArgumentError("correlation only defined for non-empty vectors"))

    @inbounds begin
        # Initialize the accumulators
        xx = zero(sqrt(abs2(x[1])))
        yy = zero(sqrt(abs2(y[1])))
        xy = zero(x[1] * y[1]')

        @simd for i in eachindex(x, y)
            xi = x[i] - mx
            yi = y[i] - my
            xx += abs2(xi)
            yy += abs2(yi)
            xy += xi * yi'
        end
    end
    return clampcor(xy / max(xx, yy) / sqrt(min(xx, yy) / max(xx, yy)))
end

corm(x::AbstractVecOrMat, xmean, y::AbstractVecOrMat, ymean, vardim::Int=1) =
    corzm(x .- xmean, y .- ymean, vardim)

# cor
"""
    cor(x::AbstractVector)

Return the number one.
"""
cor(x::AbstractVector) = one(real(eltype(x)))

"""
    cor(X::AbstractMatrix; dims::Int=1)

Compute the Pearson correlation matrix of the matrix `X` along the dimension `dims`.
"""
cor(X::AbstractMatrix; dims::Int=1) = corm(X, _vmean(X, dims), dims)

"""
    cor(x::AbstractVector, y::AbstractVector)

Compute the Pearson correlation between the vectors `x` and `y`.
"""
cor(x::AbstractVector, y::AbstractVector) = corm(x, Base.mean(x), y, Base.mean(y))

"""
    cor(X::AbstractVecOrMat, Y::AbstractVecOrMat; dims=1)

Compute the Pearson correlation between the vectors or matrices `X` and `Y` along the dimension `dims`.
"""
cor(x::AbstractVecOrMat, y::AbstractVecOrMat; dims::Int=1) =
    corm(x, _vmean(x, dims), y, _vmean(y, dims), dims)

function cov(X::SparseMatrixCSC; dims::Int=1, corrected::Bool=true)
    vardim = dims
    a, b = size(X)
    n, p = vardim == 1 ? (a, b) : (b, a)

    # The covariance can be decomposed into two terms
    # 1/(n - 1) ∑ (x_i - x̄)*(x_i - x̄)' = 1/(n - 1) (∑ x_i*x_i' - n*x̄*x̄')
    # which can be evaluated via a sparse matrix-matrix product

    # Compute ∑ x_i*x_i' = X'X using sparse matrix-matrix product
    out = Matrix(unscaled_covzm(X, vardim))

    # Compute x̄
    x̄ᵀ = mean(X, dims=vardim)

    # Subtract n*x̄*x̄' from X'X
    @inbounds for j in 1:p, i in 1:p
        out[i,j] -= x̄ᵀ[i] * x̄ᵀ[j]' * n
    end

    # scale with the sample size n or the corrected sample size n - 1
    return rmul!(out, inv(n - corrected))
end
# This is the function that does the reduction underlying var/std
function centralize_sumabs2!(R::AbstractArray{S}, A::SparseMatrixCSC{Tv,Ti}, means::AbstractArray) where {S,Tv,Ti}
    lsiz = Base.check_reducedims(R,A)
    size(means) == size(R) || error("size of means must match size of R")
    isempty(R) || fill!(R, zero(S))
    isempty(A) && return R

    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    m = size(A, 1)
    n = size(A, 2)

    if size(R, 1) == size(R, 2) == 1
        # Reduction along both columns and rows
        R[1, 1] = centralize_sumabs2(A, means[1])
    elseif size(R, 1) == 1
        # Reduction along rows
        @inbounds for col = 1:n
            mu = means[col]
            r = convert(S, (m-colptr[col+1]+colptr[col])*abs2(mu))
            @simd for j = colptr[col]:colptr[col+1]-1
                r += abs2(nzval[j] - mu)
            end
            R[1, col] = r
        end
    elseif size(R, 2) == 1
        # Reduction along columns
        rownz = fill(convert(Ti, n), m)
        @inbounds for col = 1:n
            @simd for j = colptr[col]:colptr[col+1]-1
                row = rowval[j]
                R[row, 1] += abs2(nzval[j] - means[row])
                rownz[row] -= 1
            end
        end
        for i = 1:m
            R[i, 1] += rownz[i]*abs2(means[i])
        end
    else
        # Reduction along a dimension > 2
        @inbounds for col = 1:n
            lastrow = 0
            @simd for j = colptr[col]:colptr[col+1]-1
                row = rowval[j]
                for i = lastrow+1:row-1
                    R[i, col] = abs2(means[i, col])
                end
                R[row, col] = abs2(nzval[j] - means[row, col])
                lastrow = row
            end
            for i = lastrow+1:m
                R[i, col] = abs2(means[i, col])
            end
        end
    end
    return R
end


"""
    linreg(x, y)

Perform simple linear regression using Ordinary Least Squares. Returns `a` and `b` such
that `a + b*x` is the closest straight line to the given points `(x, y)`, i.e., such that
the squared error between `y` and `a + b*x` is minimized.

# Examples
```julia-repl
julia> x = 1:10
1:10

julia> y = 2 .* x .+ 5
7:2:25

julia> a, b = linreg(x, y)
(5.000000000000002, 1.9999999999999996)
```
"""
function linreg(x::AbstractVector, y::AbstractVector)
    # Least squares given
    # Y = a + b*X
    # where
    # b = cov(X, Y)/var(X)
    # a = mean(Y) - b*mean(X)
    if size(x) != size(y)
        throw(DimensionMismatch("x has size $(size(x)) and y has size $(size(y)), " *
            "but these must be the same size"))
    end
    mx = mean(x)
    my = mean(y)
    # don't need to worry about the scaling (n vs n - 1) since they cancel in the ratio
    b = covm(x, mx, y, my)/varm(x, mx)
    a = my - b*mx
    return (a, b)
end
