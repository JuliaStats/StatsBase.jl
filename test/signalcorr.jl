# Test corr.jl
#
# Many test cases are migrated from test/01.jl in the old version
# The reference results are generated from R.
#

using StatsBase
using Test

# random data for testing

x = [0.5316732188954801 + 0.015032026010878084im -0.8051436916816284 - 0.730238071058262im
     -0.26236085703683865 + 0.1561582684776427im -0.7733975759329091 + 0.7517417038318899im
     -0.12750317530769195 + 0.8359780109572358im -0.1426497126379726 - 0.1976368614625515im
     -0.3779012097765515 - 0.15588720329528907im  0.4625694946681199 - 0.6687893051682382im
     -0.14198553599197566 - 0.911613885474317im  -0.38003834312962825 - 0.16244828932805647im
     0.6809910683979451 - 1.1951424604960583im    0.36127363787490463 + 1.3578919531660747im
     -0.7024909626802256 - 0.2915841698461933im   0.8432078641246675 - 1.373657703019187im
     -0.5684618031428246 + 0.3729428437838556im  -0.18853209977543475 - 0.9698921354142658im
     1.2992421672413557 - 0.17943637636017656im  -1.7801079183259816 - 1.3063359194617465im
     -0.26481424468298875 - 1.1636864449113618im -0.2605299632484854 - 0.7356290016171836im]

x1 = view(x, :, 1)
x2 = view(x, :, 2)
cmplx_x = convert(AbstractMatrix{Complex}, x)
cmplx_x1 = convert(AbstractVector{Complex}, x1)
cmplx_x2 = convert(AbstractVector{Complex}, x2)
# autocov & autocorr

@test autocov([1:5;])          ≈ [2.0, 0.8, -0.2, -0.8, -0.8]
@test autocor([1, 2, 3, 4, 5]) ≈ [1.0, 0.4, -0.1, -0.4, -0.4]
@test_throws MethodError autocov([1 missing 2 3 4 5])
@test_throws MethodError autocor([1 missing 2 3 4 5])

acovx1 =  [0.755284179631112 + 0.0im,
           -0.005333112170365584 - 0.18633291805458002im,
           -0.28638133506011604 + 0.1397225175604478im,
           -0.05476759020476188 + 0.12991227531087618im,
           -0.00499902418309206 + 0.023463421186162473im,
           0.11854989361702409 - 0.0028772373657594066im,
           -0.005123812003979093 - 0.08041878599400383im,
           -0.14090692583679154 + 0.035230042290069444im,
           0.039899180711674635 + 0.004918236900383659im,
           -0.03857936468514856 - 0.04063999068714769im]

@test autocov(x1) ≈ acovx1
@test autocov(cmplx_x1) ≈ acovx1
@test autocov(x)  ≈ [autocov(x1) autocov(x2)]
@test autocov(cmplx_x)  ≈ [autocov(cmplx_x1) autocov(cmplx_x2)]

acorx1 = [1.0 + 0.0im,
          -0.007061066965510023 - 0.24670570770539213im,
          -0.3791703080554228 + 0.18499330626611243im,
          -0.07251256107537019 + 0.17200449686941222im,
          -0.006618732813301651 + 0.031065686027770673im,
          0.1569606471499575 - 0.0038094765432061303im,
          -0.00678395250709688 - 0.106474871528861im,
          -0.18656146869859214 + 0.04664475073114351im,
          0.05282671316002114 + 0.0065117700502952056im,
          -0.051079270194684986 - 0.0538075492419246im]

@test autocor(x1) ≈ acorx1
@test autocor(cmplx_x1) ≈ acorx1
@test autocor(x)  ≈ [autocor(x1) autocor(x2)]
@test autocor(cmplx_x)  ≈ [autocor(cmplx_x1) autocor(cmplx_x2)]


# crosscov & crosscor

rcov0 = [0.96 + 0.32im,
         -0.96 - 0.32im,
         0.24 + 0.08im,
         -1.44 - 0.48im,
         0.0 + 0.0im,
         1.44 + 0.48im,
         -0.24 - 0.08im,
         0.96 + 0.32im,
         -0.96 - 0.32im]

@test crosscov([1+9im, 2+7im, 3+5im, 4+3im, 5+1im], [1-im, -1+im, 1-im, -1+im, 1-im]) ≈ rcov0
@test crosscov([1:5 ;], [1:5;]) ≈ [-0.8, -0.8, -0.2, 0.8, 2.0, 0.8, -0.2, -0.8, -0.8]

c11 = crosscov(x1, x1)
c12 = crosscov(x1, x2)
c21 = crosscov(x2, x1)
c22 = crosscov(x2, x2)
@test crosscov(cmplx_x1, cmplx_x2) ≈ c12

@test crosscov(x,  x1) ≈ [c11 c21]
@test crosscov(cmplx_x, cmplx_x1) ≈ [c11 c21]
@test crosscov(x1, x)  ≈ [c11 c12]
@test crosscov(cmplx_x1, cmplx_x)  ≈ [c11 c12]
@test crosscov(x,  x)  ≈ cat([c11 c21], [c12 c22], dims=3)
@test crosscov(cmplx_x,  cmplx_x)  ≈ cat([c11 c21], [c12 c22], dims=3)

# issue #805: avoid converting one input to the other's eltype
@test crosscov([34566.5345, 3466.4566], Float16[1, 10]) ≈
    crosscov(Float16[1, 10], [34566.5345, 3466.4566]) ≈
    crosscov([34566.5345, 3466.4566], Float16[1, 10])

rcor0 = [0.230940107675850,
        -0.230940107675850,
         0.057735026918963,
        -0.346410161513775,
         0.000000000000000,
         0.346410161513775,
        -0.057735026918963,
         0.230940107675850,
        -0.230940107675850]

@test crosscor([1, 2, 3, 4, 5], [1, -1, 1, -1, 1]) ≈ rcor0
@test crosscor([1:5;], [1:5;]) ≈ [-0.4, -0.4, -0.1, 0.4, 1.0, 0.4, -0.1, -0.4, -0.4]

c11 = crosscor(x1, x1)
c12 = crosscor(x1, x2)
c21 = crosscor(x2, x1)
c22 = crosscor(x2, x2)
@test crosscor(cmplx_x1, cmplx_x2) ≈ c12

@test crosscor(x,  x1) ≈ [c11 c21]
@test crosscor(cmplx_x, cmplx_x1) ≈ [c11 c21]
@test crosscor(x1, x)  ≈ [c11 c12]
@test crosscor(cmplx_x1, cmplx_x)  ≈ [c11 c12]
@test crosscor(x,  x)  ≈ cat([c11 c21], [c12 c22], dims=3)
@test crosscor(cmplx_x, cmplx_x)  ≈ cat([c11 c21], [c12 c22], dims=3)

# issue #805: avoid converting one input to the other's eltype
@test crosscor([34566.5345, 3466.4566], Float16[1, 10]) ≈
    crosscor(Float16[1, 10], [34566.5345, 3466.4566]) ≈
    crosscor([34566.5345, 3466.4566], Float16[1, 10])


## pacf least squares
pacf_ls = [-1.598495044296996e-03 - 2.915104118351207e-01im
           -5.560162016912027e-01 + 2.950837739894279e-01im
           -2.547001916363494e-02 + 2.326084658014266e-01im
           -5.427443903358727e-01 + 3.146715147305132e-01im]

@test pacf(x1, 1:4) ≈ pacf_ls[1:4]

## pacf Yule-Walker

function yulewalker_qr(v::AbstractVector)
    A = toeplitz(v)
    b = v[2:end]
    x = -A\b
end
function toeplitz(v::AbstractVector{T}) where T
    N = length(v)
    A = zeros(T, N - 1, N - 1)
    for n in 1:N-1
        A[n, n+1:end] = conj(v[2:N-n])
        A[n, 1:n] = reverse(v[1:n])
    end
    return A
end
# durbin solver 
acf = autocor(x1)
p_qr = [yulewalker_qr(acf[1:n])[n-1] for n in 2:length(acf)]
y = similar(acf[2:end])
pd = similar(acf[2:end])
StatsBase.durbin!(acf[2:end], y, pd)
@test p_qr ≈ pd

@test pacf(x1, 1:4, method=:yulewalker) ≈ -p_qr[1:4]
