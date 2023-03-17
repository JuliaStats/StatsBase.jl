
## coefficient tables with specialized show method

mutable struct CoefTable
    cols::Vector
    colnms::Vector
    rownms::Vector
    pvalcol::Int
    teststatcol::Int
    function CoefTable(cols::Vector,colnms::Vector,rownms::Vector,
                       pvalcol::Int=0,teststatcol::Int=0)
        nc = length(cols)
        nrs = map(length,cols)
        nr = nrs[1]
        length(colnms) in [0,nc] || throw(ArgumentError("colnms should have length 0 or $nc"))
        length(rownms) in [0,nr] || throw(ArgumentError("rownms should have length 0 or $nr"))
        all(nrs .== nr) || throw(ArgumentError("Elements of cols should have equal lengths, but got $nrs"))
        pvalcol in 0:nc || throw(ArgumentError("pvalcol should be between 0 and $nc"))
        teststatcol in 0:nc || throw(ArgumentError("teststatcol should be between 0 and $nc"))
        new(cols,colnms,rownms,pvalcol,teststatcol)
    end

    function CoefTable(mat::Matrix,colnms::Vector,rownms::Vector,
                       pvalcol::Int=0,teststatcol::Int=0)
        nc = size(mat,2)
        cols = Any[mat[:, i] for i in 1:nc]
        CoefTable(cols,colnms,rownms,pvalcol,teststatcol)
    end
end

Base.length(ct::CoefTable) = length(ct.cols[1])
function Base.eltype(ct::CoefTable)
    names = isempty(ct.rownms) ?
        tuple(Symbol.(ct.colnms)...) :
        tuple(Symbol("Name"), Symbol.(ct.colnms)...)
    types = isempty(ct.rownms) ?
        Tuple{eltype.(ct.cols)...} :
        Tuple{eltype(ct.rownms), eltype.(ct.cols)...}
    NamedTuple{names, types}
end

function Base.iterate(ct::CoefTable, i::Integer=1)
    if i in 1:length(ct)
        cols = getindex.(ct.cols, Ref(i))
        nt = isempty(ct.rownms) ?
            eltype(ct)(tuple(cols...)) :
            eltype(ct)(tuple(ct.rownms[i], cols...))
        (nt, i+1)
    else
        nothing
    end
end

"""
Show a p-value using 6 characters, either using the standard 0.XXXX
representation or as <Xe-YY.
"""
struct PValue <: Real
    v::Real
    function PValue(v::Real)
        0 <= v <= 1 || isnan(v) || error("p-values must be in [0; 1]")
        new(v)
    end
end
PValue(p::PValue) = p

function show(io::IO, pv::PValue)
    v = pv.v
    if isnan(v)
        @printf(io,"%d", v)
    elseif v >= 1e-4
        @printf(io,"%.4f", v)
    else
        @printf(io,"<1e%2.2d", ceil(Integer, max(nextfloat(log10(v)), -99)))
    end
end

"""Show a test statistic using 2 decimal digits"""
struct TestStat <: Real
    v::Real
end

show(io::IO, x::TestStat) = @printf(io, "%.2f", x.v)
TestStat(x::TestStat) = x

float(x::Union{TestStat, PValue}) = float(x.v)

for op in [:(==), :<, :≤, :(isless), :(isequal)] # isless and < to place nice with NaN
    @eval begin
        Base.$op(x::Union{TestStat, PValue}, y::Real) = $op(x.v, y)
        Base.$op(y::Real, x::Union{TestStat, PValue}) = $op(y, x.v)
        Base.$op(x1::Union{TestStat, PValue}, x2::Union{TestStat, PValue}) = $op(x1.v, x2.v)
    end
end

Base.hash(x::Union{TestStat, PValue}, h::UInt) = hash(x.v, h)

# necessary to avoid a method ambiguity with isless(::TestStat, NaN)
Base.isless(x::Union{TestStat, PValue}, y::AbstractFloat) = isless(x.v, y)
Base.isless(y::AbstractFloat, x::Union{TestStat, PValue},) = isless(y, x.v)
Base.isequal(y::AbstractFloat, x::Union{TestStat, PValue}) = isequal(y, x.v)
Base.isequal(x::Union{TestStat, PValue}, y::AbstractFloat) = isequal(x.v, y)

Base.isapprox(x::Union{TestStat, PValue}, y::Real; kwargs...) = isapprox(x.v, y; kwargs...)
Base.isapprox(y::Real, x::Union{TestStat, PValue}; kwargs...) = isapprox(y, x.v; kwargs...)
Base.isapprox(x1::Union{TestStat, PValue}, x2::Union{TestStat, PValue}; kwargs...) = isapprox(x1.v, x2.v; kwargs...)


"""Wrap a string so that show omits quotes"""
struct NoQuote
    s::String
end

show(io::IO, n::NoQuote) = print(io, n.s)

function show(io::IO, ct::CoefTable)
    cols = ct.cols; rownms = ct.rownms; colnms = ct.colnms;
    nc = length(cols)
    nr = length(cols[1])
    if length(rownms) == 0
        rownms = [lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    mat = [j == 1 ? NoQuote(rownms[i]) :
           j-1 == ct.pvalcol ? NoQuote(sprint(show, PValue(cols[j-1][i]))) :
           j-1 in ct.teststatcol ? TestStat(cols[j-1][i]) :
           cols[j-1][i] isa AbstractString ? NoQuote(cols[j-1][i]) : cols[j-1][i]
           for i in 1:nr, j in 1:nc+1]
    # Code inspired by print_matrix in Base
    io = IOContext(io, :compact=>true, :limit=>false)
    A = Base.alignment(io, mat, 1:size(mat, 1), 1:size(mat, 2),
                       typemax(Int), typemax(Int), 3)
    nmswidths = pushfirst!(length.(colnms), 0)
    A = [nmswidths[i] > sum(A[i]) ? (A[i][1]+nmswidths[i]-sum(A[i]), A[i][2]) : A[i]
         for i in 1:length(A)]
    totwidth = sum(sum.(A)) + 2 * (length(A) - 1)
    println(io, repeat('─', totwidth))
    print(io, repeat(' ', sum(A[1])))
    for j in 1:length(colnms)
        print(io, "  ", lpad(colnms[j], sum(A[j+1])))
    end
    println(io, '\n', repeat('─', totwidth))
    for i in 1:size(mat, 1)
        Base.print_matrix_row(io, mat, A, i, 1:size(mat, 2), "  ")
        i != size(mat, 1) && println(io)
    end
    print(io, '\n', repeat('─', totwidth))
    nothing
end

function show(io::IO, ::MIME"text/markdown", ct::CoefTable)
    cols = ct.cols; rownms = ct.rownms; colnms = ct.colnms;
    nc = length(cols)
    nr = length(cols[1])
    if length(rownms) == 0
        rownms = [lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    mat = [j == 1 ? NoQuote(rownms[i]) :
           j-1 == ct.pvalcol ? NoQuote(sprint(show, PValue(cols[j-1][i]))) :
           j-1 in ct.teststatcol ? TestStat(cols[j-1][i]) :
           cols[j-1][i] isa AbstractString ? NoQuote(cols[j-1][i]) : cols[j-1][i]
           for i in 1:nr, j in 1:nc+1]
    # Code inspired by print_matrix in Base
    io = IOContext(io, :compact=>true, :limit=>false)
    A = Base.alignment(io, mat, 1:size(mat, 1), 1:size(mat, 2),
                       typemax(Int), typemax(Int), 3)
    nmswidths = pushfirst!(length.(colnms), 0)
    A = [nmswidths[i] > sum(A[i]) ? (A[i][1]+nmswidths[i]-sum(A[i]), A[i][2]) : A[i]
         for i in 1:length(A)]

    # not using Markdown stdlib here because that won't give us nice decimal
    # alignment (even if that is lost when rendering to HTML, it's still nice
    # when looking at the markdown itself)

    print(io, '|', ' '^(sum(A[1])+1))
    for j in 1:length(colnms)
        print(io, " | ", lpad(colnms[j], sum(A[j+1])))
    end

    println(io, " |")
    print(io, '|', rpad(':', sum(A[1])+2, '-'))
    for j in 1:length(colnms)
        _pad = j-1 in [ct.teststatcol; ct.pvalcol] ? rpad : lpad
        print(io, '|', _pad(':', sum(A[j+1])+2, '-'))
    end
    println(io, '|')

    for i in 1:size(mat, 1)
        print(io, "| ")
        Base.print_matrix_row(io, mat, A, i, 1:size(mat, 2), " | ")
        print(io, " |")
        i != size(mat, 1) && println(io)
    end

    nothing
end

"""
    ConvergenceException(iters::Int, lastchange::Real=NaN, tol::Real=NaN)

The fitting procedure failed to converge in `iters` number of iterations,
i.e. the `lastchange` between the cost of the final and penultimate iteration was greater than
specified tolerance `tol`.
"""
struct ConvergenceException{T<:Real} <: Exception
    iters::Int
    lastchange::T
    tol::T
    msg::String
    function ConvergenceException{T}(iters, lastchange::T, tol::T, msg::String) where T<:Real
        if tol > lastchange
            throw(ArgumentError("Change must be greater than tol."))
        else
            new(iters, lastchange, tol, msg)
        end
    end
end

ConvergenceException(iters, lastchange::T=NaN, tol::T=NaN,
                     msg::AbstractString="") where {T<:Real} =
    ConvergenceException{T}(iters, lastchange, tol, String(msg))

function Base.showerror(io::IO, ce::ConvergenceException)
    print(io, "failure to converge after $(ce.iters) iterations.")
    if !isnan(ce.lastchange)
        print(io, " Last change ($(ce.lastchange)) was greater than tolerance ($(ce.tol)).")
    end
    if !isempty(ce.msg)
        print(io, ' ', ce.msg)
    end
end
