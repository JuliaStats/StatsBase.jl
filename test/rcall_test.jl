using RCall

_r_ccf(x, y, lm) = rcopy(rcall(:ccf, x, y, plot=false, lag=lm, type="covariance"))[:acf]
_r_acf(x, lm) = rcopy(rcall(:acf, x, plot=false, lag=lm, type="covariance"))[:acf]
_r_pacf(x, lm) = rcopy(rcall(:acf, x, plot=false, lag=lm))[:acf]

function _ccf(is_auto::Bool, x::AbstractVector{<:Real}...)
    lx = size(x[1], 1)
    lag_max =  min(lx - 1, round(Int, 10*log10(lx)))
    cor = reverse(dropdims(_r_ccf(x[1], x[end], lag_max), dims = (2, 3)))
    is_auto ? cor[lx:end] : cor
end

r_autocov(x::AbstractVector) = r_crosscov(x, x, true)

function r_crosscov(x::AbstractVector, y::AbstractVector, is_auto=false)
    cv(x,y) = _ccf(is_auto, x, y)
    cc_re = cv(real.(x), real.(y)) + cv(imag.(x), imag.(y))
    cc_im = cv(real.(x), imag.(y)) - cv(imag.(x), real.(y))
    cc_re + im*cc_im
end

r_autocor(x::AbstractVector) = r_crosscor(x, x, true)

function r_crosscor(x::AbstractVector, y::AbstractVector, is_auto=false)
    cc = r_crosscov(x, y, is_auto)
    sc(z) = _r_ccf(real.(z), real.(z), 0) + _r_ccf(imag.(z), imag.(z), 0)
    cc/sqrt(sc(x)*sc(y))
end
