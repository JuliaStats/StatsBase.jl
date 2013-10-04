# Test taken from http://www.stanford.edu/~clint/bench/wilk.txt

using Test
using DataFrames
using Stats
using GLM

testeps = sqrt(eps())

nasty = DataFrame(quote
					label = ["One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"]
					x = [1:9]
					zero = repeat([0],9)
					miss = repeat([NA], 9)
					big = 99999990 + [1:9]
					little = (99999990 + [1:9])/10^8
					huge = [1.:9]*1e12
					tiny = [1.:9]*1e-12
					round = [0.5:8.5]
				end)

println(nasty)
println("\nII Real Numbers:\nII A")
print("Test rounding: ")
@test [@sprintf("%1.0f", x) for x in nasty["round"]] == ["1","2","3","4","5","6","7","8","9"]
println("OK")
print("Test math: ")
@test int(2.6*7 - 0.2) == 18
@test 2 - int(exp(log(sqrt(2)*sqrt(2)))) == 0
@test int(3 - exp(log(sqrt(2)*sqrt(2)))) == 1
println("OK")

print("Test means: ")
for vars in colnames(nasty)[2:]
	if vars == "miss"
		@test isna(mean(nasty[vars]))
	else
		@test_approx_eq mean(nasty[vars]) nasty[vars][5]
	end
end
println("OK")

print("Test standard deviation: ")
for vars in colnames(nasty)[[2,5:9]]
	@test (@sprintf("%.9e", std(nasty[vars])))[1:10] == "2.73861278"
end
println("OK")
# Note: Origin had one more digit but little fails then. R gives the same as Julia. Stata and SAS less precise returning 2.738612714381 and 2.738612780003404")

println("\nII D")
print("Test correlation: ")
cn = colnames(nasty)[[2,5:9]]
for i in 1:5
	for j = i+1:6
	@test_approx_eq cor(nasty[cn[i]], nasty[cn[j]]) 1
	end
end
println("OK")

print("Test autocorrelation: ")
# x was sampled by calling x = randn(10, 2)
x = [-2.133252557240862 -.7445937365828654;
.1775816414485478 -.5834801838041446;
-.6264517920318317 -.68444205333293;
-.8809042583216906 .9071671734302398;
.09251017186697393 -1.0404476733379926;
-.9271887119115569 -.620728578941385;
3.355819743178915 -.8325051361909978;
-.2834039258495755 -.22394811874731657;
.5354280026977677 .7481337671592626;
.39182285417742585 .3085762550821047]
# Set the number of printed digits in R to 20 by running the command options(digits=20)
# "racfx11" was computed by calling R's acf function on x: acf(x[, 1], plot=FALSE)$acf[, 1, 1]
racfx11 = [0.999999999999999888978, -0.221173011668873431557,  0.229321981664153962122,  0.019505581764945757045,
  -0.139015765538446717242, 0.125681062460244019618, -0.427909344123907742219,  0.021699096507690283225,
  -0.059889541590524189574, -0.048220059475281865091]
# "jacfx11" are computed by calling Julia's acf function on x
jacfx11 = acf(x[:, 1])
for i =1:length(racfx11)
@test_approx_eq racfx11[i] jacfx11[i]
end
println("OK")

print("Test cross-correlation: ")
# "racfx12" was computed by calling R's acf function on x: acf(x, plot=FALSE)$acf[, 2, 1]
racfx12 = [-0.0785530595168460604727, 0.1363402071384945957178, 0.5695846902378886023044, -0.1015722302688646383473,
  0.1601773904236522549915, -0.0251707633078918115166, 0.0079618741584954952351]
# "jacfx12" are computed by calling Julia's acf function on x
jacfx12 = ccf(x[:, 1], x[:, 2], 0:6)
for i =1:length(racfx12)
@test_approx_eq racfx12[i] jacfx12[i]
end
println("OK")

print("Test autocorrelation vs cross-correlation: ")
@test acf(x[:, 1]) == ccf(x[:, 1], x[:, 1])

print("Test autocovariance: ")
# "racvx11" was computed by calling R's acf function on x: acf(x[, 1], plot=FALSE, type="covariance")$acf[, 1, 1]
racvx11 =  [1.839214242630635709475, -0.406784553146903871124, 0.421772254824993531042, 0.035874943792884653182,
  -0.255679775928512320604, 0.231154400105831353551, -0.787016960267425180753, 0.039909287349160660341,
  -0.110149697877911914579, -0.088687020167434751916]
# "jacvx11" are computed by calling Julia's acf function on x
jacvx11 = acf(x[:, 1], correlation=false)
for i =1:length(racvx11)
@test_approx_eq racvx11[i] jacvx11[i]
end
println("OK")

print("Test cross-covariance: ")
# "racvx12" was computed by calling R's acf function on x: acf(x, plot=FALSE, type="covariance")$acf[, 2, 1]
racvx12 = [-0.0697521748336961816550, 0.1210649976420989648584, 0.5057698724968141545943, -0.0901923363334168753935,
  0.1422315236345411126884, -0.0223506951065748533936, 0.0070698460597599819752]
# "jacfx12" are computed by calling Julia's acf function on x
jacvx12 = ccf(x[:, 1], x[:, 2], 0:6, correlation=false)
for i =1:length(racvx12)
@test_approx_eq racvx12[i] jacvx12[i]
end
println("OK")

print("Test autocorrelation with demean=false: ")
# "racfx11nodemean" was computed by calling R's acf function on x: acf(x[, 1], plot=FALSE, demean=FALSE)$acf[, 1, 1]
racfx11nodemean = [1.000000000000000000000, -0.223355812053329189082, 0.228124838820096320635, 0.016984315796356772021,
  -0.137403649655494369819, 0.125861757190160739039, -0.428765059298546857836, 0.024683201099522964622,
  -0.058291459424761618568, -0.045424485824113770838]
# "jacfx11nodemean" are computed by calling Julia's acf function on x
jacfx11nodemean = acf(x[:, 1], demean=false)
for i =1:length(racfx11nodemean)
@test_approx_eq racfx11nodemean[i] jacfx11nodemean[i]
end
println("OK")

print("Test cross-correlation with demean=false: ")
# "racfx12nodemean" was computed by calling R's acf function on x: acf(x, plot=FALSE, demean=FALSE)$acf[, 2, 1]
racfx12nodemean = [-0.063791841616291797279, 0.143906622654901311664, 0.557311279399980707971, -0.070174737610974355362,
  0.270818343664623484290, 0.071161936583426677050, 0.103265547537476284901]
# "jacfx12" are computed by calling Julia's acf function on x
jacfx12nodemean = ccf(x[:, 1], x[:, 2], 0:6, demean=false)
for i =1:length(racfx12nodemean)
@test_approx_eq racfx12nodemean[i] jacfx12nodemean[i]
end
println("OK")

print("Test autocorrelation vs cross-correlation with demean=false: ")
@test acf(x[:, 1], demean=false) == ccf(x[:, 1], x[:, 1], demean=false)

print("Test autocovariance with demean=false: ")
# "racvx11" was computed by calling R's acf function on x:
# acf(x[, 1], plot=FALSE, type="covariance", demean=FALSE)$acf[, 1, 1]
racvx11nodemean =  [1.840102514084350771029, -0.410997591294682773633, 0.419773089437946556046, 0.031252882196878647991,
  -0.252836801175440550882, 0.231598535832688912084, -0.788971663566781833410, 0.045419620398881817291,
  -0.107262261037149780885, -0.083585710565940704586]
# "jacvx11" are computed by calling Julia's acf function on x
jacvx11nodemean = acf(x[:, 1], correlation=false, demean=false)
for i =1:length(racvx11nodemean)
@test_approx_eq racvx11nodemean[i] jacvx11nodemean[i]
end
println("OK")

print("Test cross-covariance with demean=false: ")
# "racvx12" was computed by calling R's acf function on x:
# acf(x, plot=FALSE, type="covariance", demean=FALSE)$acf[, 2, 1]
racvx12nodemean = [-0.061507621146693322589, 0.138753699571784711031, 0.537355407299582976677, -0.067661962183297230666,
  0.261121040867466125412, 0.068613812119833542114, 0.099567875993390383971]
# "jacfx12" are computed by calling Julia's acf function on x
jacvx12nodemean = ccf(x[:, 1], x[:, 2], 0:6, correlation=false, demean=false)
for i =1:length(racvx12nodemean)
@test_approx_eq racvx12nodemean[i] jacvx12nodemean[i]
end
println("OK")

print("Test spearman correlation: ")
cn = colnames(nasty)[[2,5:9]]
for i in 1:5
	for j = i+1:6
	@test_approx_eq cor_spearman(nasty[cn[i]], nasty[cn[j]]) 1
	end
end
println("OK")

println("\nII F")
print("Testing regression: ")
@test_approx_eq coef(lm(:(big~x), nasty)) [99999990, 1]
println("OK")

println("\nIV Regression:\nIV A")
nasty["x1"] = nasty["x"]
nasty["x2"] = nasty["x"].^2
nasty["x3"] = nasty["x"].^3
nasty["x4"] = nasty["x"].^4
nasty["x5"] = nasty["x"].^5
nasty["x6"] = nasty["x"].^6
nasty["x7"] = nasty["x"].^7
nasty["x8"] = nasty["x"].^8
nasty["x9"] = nasty["x"].^9
lm(:(x1~x2+x3+x4+x5+x6+x7+x8+x9), nasty)
@test_approx_eq coef(lm(:(x~x), nasty)) [0,1]
println("OK")