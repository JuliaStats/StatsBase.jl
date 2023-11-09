# Test corr.jl
#
# Many test cases are migrated from test/01.jl in the old version
# The reference results are generated from R.
#

using StatsBase
using Test

# random data for testing

x = [-2.133252557240862    -.7445937365828654;
       .1775816414485478   -.5834801838041446;
      -.6264517920318317   -.68444205333293;
      -.8809042583216906    .9071671734302398;
       .09251017186697393 -1.0404476733379926;
      -.9271887119115569   -.620728578941385;
      3.355819743178915    -.8325051361909978;
      -.2834039258495755   -.22394811874731657;
       .5354280026977677    .7481337671592626;
       .39182285417742585   .3085762550821047]

x1 = view(x, :, 1)
x2 = view(x, :, 2)
realx = convert(AbstractMatrix{Real}, x)
realx1 = convert(AbstractVector{Real}, x1)
realx2 = convert(AbstractVector{Real}, x2)

# autocov & autocorr

@test autocov([1:5;])          ≈ [2.0, 0.8, -0.2, -0.8, -0.8]
@test autocor([1, 2, 3, 4, 5]) ≈ [1.0, 0.4, -0.1, -0.4, -0.4]

racovx1 =  [1.839214242630635709475,
           -0.406784553146903871124,
            0.421772254824993531042,
            0.035874943792884653182,
           -0.255679775928512320604,
            0.231154400105831353551,
           -0.787016960267425180753,
            0.039909287349160660341,
           -0.110149697877911914579,
           -0.088687020167434751916]

@test autocov(x1) ≈ racovx1
@test autocov(realx1) ≈ racovx1
@test autocov(x)  ≈ [autocov(x1) autocov(x2)]
@test autocov(realx)  ≈ [autocov(realx1) autocov(realx2)]

racorx1 = [0.999999999999999888978,
          -0.221173011668873431557,
           0.229321981664153962122,
           0.019505581764945757045,
          -0.139015765538446717242,
           0.125681062460244019618,
          -0.427909344123907742219,
           0.021699096507690283225,
          -0.059889541590524189574,
          -0.048220059475281865091]

@test autocor(x1) ≈ racorx1
@test autocor(realx1) ≈ racorx1
@test autocor(x)  ≈ [autocor(x1) autocor(x2)]
@test autocor(realx)  ≈ [autocor(realx1) autocor(realx2)]


# crosscov & crosscor

rcov0 = [0.320000000000000006661,
        -0.319999999999999951150,
         0.080000000000000029421,
        -0.479999999999999982236,
         0.000000000000000000000,
         0.479999999999999982236,
        -0.080000000000000029421,
         0.319999999999999951150,
        -0.320000000000000006661]

@test crosscov([1, 2, 3, 4, 5], [1, -1, 1, -1, 1]) ≈ rcov0
@test crosscov([1:5;], [1:5;]) ≈ [-0.8, -0.8, -0.2, 0.8, 2.0, 0.8, -0.2, -0.8, -0.8]

c11 = crosscov(x1, x1)
c12 = crosscov(x1, x2)
c21 = crosscov(x2, x1)
c22 = crosscov(x2, x2)
@test crosscov(realx1, realx2) ≈ c12

@test crosscov(x,  x1) ≈ [c11 c21]
@test crosscov(realx, realx1) ≈ [c11 c21]
@test crosscov(x1, x)  ≈ [c11 c12]
@test crosscov(realx1, realx)  ≈ [c11 c12]
@test crosscov(x,  x)  ≈ cat([c11 c21], [c12 c22], dims=3)
@test crosscov(realx,  realx)  ≈ cat([c11 c21], [c12 c22], dims=3)

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
@test crosscor(realx1, realx2) ≈ c12

@test crosscor(x,  x1) ≈ [c11 c21]
@test crosscor(realx, realx1) ≈ [c11 c21]
@test crosscor(x1, x)  ≈ [c11 c12]
@test crosscor(realx1, realx)  ≈ [c11 c12]
@test crosscor(x,  x)  ≈ cat([c11 c21], [c12 c22], dims=3)
@test crosscor(realx, realx)  ≈ cat([c11 c21], [c12 c22], dims=3)

# issue #805: avoid converting one input to the other's eltype
@test crosscor([34566.5345, 3466.4566], Float16[1, 10]) ≈
    crosscor(Float16[1, 10], [34566.5345, 3466.4566]) ≈
    crosscor([34566.5345, 3466.4566], Float16[1, 10])


## pacf

rpacfr =  [-0.218158122381419,
            0.195015316828711,
            0.144315804606139,
           -0.199791229449779]

@test pacf(x[:,1], 1:4) ≈ rpacfr

rpacfy = [-0.221173011668873,
           0.189683314308021,
           0.111857020733719,
          -0.175020669835420]

@test pacf(x[:,1], 1:4, method=:yulewalker) ≈ rpacfy

rho, sigma = yule_walker([1.0, 2, 3]; order=1)
@test rho == [0.0]

rho, sigma = yule_walker([1.0, 2, 3]; order=2)
@test rho == [0.0, -1.5]

x = [0.9901178, -0.74795127, 0.44612542, 1.1362954, -0.04040932]
rho, sigma = yule_walker(x; order=3, method="mle")
@test rho ≈ [-0.9418963, -0.90335955, -0.33267884]
@test sigma ≈ 0.44006365345695164

rho, sigma = yule_walker(x; order=3)
@test isapprox(rho, [0.10959317, 0.05242324, 1.06587676], atol=tol)
@test isapprox(sigma, 0.15860522671108127, atol=tol)

rho, sigma = yule_walker(x; order=5, method="mle")
@test isapprox(
    rho, [-1.24209771, -1.56893346, -1.16951484, -0.79844781, -0.27598787], atol=tol
)
@test isapprox(sigma, 0.3679474002175471, atol=tol)

x = [
    0.9901178,
    -0.74795127,
    0.44612542,
    1.1362954,
    -0.04040932,
    0.28625813,
    0.88901716,
    -0.1079814,
    -0.33231995,
    0.4607741,
]

rho, sigma = yule_walker(x; order=3, method="mle")
@test isapprox(
    rho, [-0.4896151627237206, -0.5724647370433921, 0.09083516892540627], atol=tol
)
@test isapprox(sigma, 0.4249693094713215, atol=tol)

x = [
    0.9901178,
    -0.74795127,
    0.44612542,
    1.1362954,
    -0.04040932,
    0.28625813,
    0.88901716,
    -0.1079814,
    -0.33231995,
    0.4607741,
    0.7729643,
    -1.0998684,
    1.098167,
    1.0105597,
    -1.3370227,
    1.239718,
    -0.01393661,
    -0.4790918,
    1.5009186,
    -1.1647809,
]

rho, sigma = yule_walker(x; order=3, method="mle")
@test isapprox(rho, [-0.82245705, -0.57029742, 0.12166898], atol=tol)
@test isapprox(sigma, 0.5203501608988023, atol=tol)

rho, sigma = yule_walker(x; order=3)
@test isapprox(rho, [-0.93458149, -0.68653741, 0.10161722], atol=tol)
@test isapprox(sigma, 0.4269012058667671, atol=tol)

rho, sigma = yule_walker(x; order=5, method="mle")
@test isapprox(
    rho, [-0.83107755, -0.56407764, 0.20950143, 0.1232321, 0.10249279], atol=tol
)
@test isapprox(sigma, 0.5172269743102993, atol=tol)

rho, sigma = yule_walker(x; order=5)
@test isapprox(
    rho, [-0.96481241, -0.65359486, 0.31587079, 0.28403115, 0.1913565], atol=tol
)
@test isapprox(sigma, 0.41677565377507053, atol=tol)
