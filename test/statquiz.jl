# Test taken from http://www.stanford.edu/~clint/bench/wilk.txt

using Base.Test
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
# "racfx11" was computed by calling R's acf function on x: acf(x[, 1], plot=FALSE)$acf[, 1, 1]
racfx11 = [1.00000000, -0.22117301, 0.22932198, 0.01950558, -0.13901577,
  0.12568106, -0.42790934, 0.02169910, -0.05988954, -0.04822006]
# "jacfx11" are computed by calling Julia's acf function on x
jacfx11 = acf(x[:, 1])
for i =1:length(racfx11)
@test_approx_eq round(racfx11[i], 6) round(jacfx11[i], 6)
end
println("OK")

print("Test cross-correlation: ")
# "racfx12" was computed by calling R's acf function on x: acf(x, plot=FALSE)$acf[, 2, 1]
racfx12 = [-0.078553060, 0.136340207, 0.569584690, -0.101572230, 0.160177390, -0.025170763, 0.007961874]
# "jacfx12" are computed by calling Julia's acf function on x
jacfx12 = acf(x[:, 1], x[:, 2], 0:6)
for i =1:length(racfx12)
@test_approx_eq round(jacfx12[i], 6) round(racfx12[i], 6)
end
println("OK")

print("Test autocorrelation vs cross-correlation: ")
@test acf(x[:, 1]) == acf(x[:, 1], x[:, 1])

print("Test autocovariance: ")
# "racvx11" was computed by calling R's acf function on x: acf(x[, 1], plot=FALSE, type="covariance")$acf[, 1, 1]
racvx11 = [1.83921424, -0.40678455, 0.42177225, 0.03587494, -0.25567978,
0.23115440, -0.78701696, 0.03990929, -0.11014970, -0.08868702]
# "jacvx11" are computed by calling Julia's acf function on x
jacvx11 = acf(x[:, 1], correlation=false)
for i =1:length(racvx11)
@test_approx_eq round(racvx11[i], 6) round(jacvx11[i], 6)
end
println("OK")

print("Test cross-covariance: ")
# "racvx12" was computed by calling R's acf function on x: acf(x, plot=FALSE, type="covariance")$acf[, 2, 1]
racvx12 = [-0.069752175, 0.121064998, 0.505769872, -0.090192336, 0.142231524, -0.022350695, 0.007069846]
# "jacfx12" are computed by calling Julia's acf function on x
jacvx12 = acf(x[:, 1], x[:, 2], 0:6, correlation=false)
for i =1:length(racvx12)
@test_approx_eq round(jacvx12[i], 6) round(racvx12[i], 6)
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