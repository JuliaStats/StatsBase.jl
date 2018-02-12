var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "StatsBase.jl Documentation",
    "title": "StatsBase.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#StatsBase.jl-Documentation-1",
    "page": "StatsBase.jl Documentation",
    "title": "StatsBase.jl Documentation",
    "category": "section",
    "text": "CurrentModule = StatsBaseStatsBase.jl is a Julia package that provides basic support for statistics. Particularly, it implements a variety of statistics-related functions, such as scalar statistics, high-order moment computation, counting, ranking, covariances, sampling, and empirical density estimation.Pages = [\"weights.md\", \"means.md\", \"scalarstats.md\", \"robust.md\", \"deviation.md\", \"cov.md\", \"counts.md\", \"ranking.md\", \"sampling.md\", \"empirical.md\", \"signalcorr.md\", \"misc.md\", \"statmodels.md\"]\nDepth = 2"
},

{
    "location": "weights.html#",
    "page": "Weight Vectors",
    "title": "Weight Vectors",
    "category": "page",
    "text": ""
},

{
    "location": "weights.html#Weight-Vectors-1",
    "page": "Weight Vectors",
    "title": "Weight Vectors",
    "category": "section",
    "text": "In statistical applications, it is not uncommon to assign weights to samples. To facilitate the use of weight vectors, we introduce the abstract type AbstractWeights for the purpose of representing weight vectors, which has two advantages:A different type AbstractWeights distinguishes the role of the weight vector from other data vectors in the input arguments.\nStatistical functions that utilize weights often need the sum of weights for various purposes. The weight vector maintains the sum of weights, so that it needn't be computed repeatedly each time the sum of weights is needed.note: Note\nThe weight vector is a light-weight wrapper of the input vector. The input vector is NOT copied during construction.\nThe weight vector maintains the sum of weights, which is computed upon construction. If the value of the sum is pre-computed, one can supply it as the second argument to the constructor and save the time of computing the sum again."
},

{
    "location": "weights.html#Implementations-1",
    "page": "Weight Vectors",
    "title": "Implementations",
    "category": "section",
    "text": "Several statistical weight types are provided which subtype AbstractWeights. The choice of weights impacts how bias is corrected in several methods. See the var, std and cov docstrings for more details."
},

{
    "location": "weights.html#AnalyticWeights-1",
    "page": "Weight Vectors",
    "title": "AnalyticWeights",
    "category": "section",
    "text": "Analytic weights describe a non-random relative importance (usually between 0 and 1) for each observation. These weights may also be referred to as reliability weights, precision weights or inverse variance weights. These are typically used when the observations being weighted are aggregate values (e.g., averages) with differing variances.w = AnalyticWeights([0.2, 0.1, 0.3])\nw = aweights([0.2, 0.1, 0.3])"
},

{
    "location": "weights.html#FrequencyWeights-1",
    "page": "Weight Vectors",
    "title": "FrequencyWeights",
    "category": "section",
    "text": "Frequency weights describe the number of times (or frequency) each observation was observed. These weights may also be referred to as case weights or repeat weights.w = FrequencyWeights([2, 1, 3])\nw = fweights([2, 1, 3])"
},

{
    "location": "weights.html#ProbabilityWeights-1",
    "page": "Weight Vectors",
    "title": "ProbabilityWeights",
    "category": "section",
    "text": "Probability weights represent the inverse of the sampling probability for each observation, providing a correction mechanism for under- or over-sampling certain population groups. These weights may also be referred to as sampling weights.w = ProbabilityWeights([0.2, 0.1, 0.3])\nw = pweights([0.2, 0.1, 0.3])"
},

{
    "location": "weights.html#Weights-1",
    "page": "Weight Vectors",
    "title": "Weights",
    "category": "section",
    "text": "The Weights type describes a generic weights vector which does not support all operations possible for FrequencyWeights, AnalyticWeights and ProbabilityWeights.w = Weights([1., 2., 3.])\nw = weights([1., 2., 3.])"
},

{
    "location": "weights.html#Methods-1",
    "page": "Weight Vectors",
    "title": "Methods",
    "category": "section",
    "text": "AbstractWeights implements the following methods:eltype\nlength\nisempty\nvalues\nsum"
},

{
    "location": "means.html#",
    "page": "Mean Functions",
    "title": "Mean Functions",
    "category": "page",
    "text": ""
},

{
    "location": "means.html#StatsBase.geomean",
    "page": "Mean Functions",
    "title": "StatsBase.geomean",
    "category": "Function",
    "text": "geomean(a)\n\nReturn the geometric mean of a real-valued array.\n\n\n\n"
},

{
    "location": "means.html#StatsBase.harmmean",
    "page": "Mean Functions",
    "title": "StatsBase.harmmean",
    "category": "Function",
    "text": "harmmean(a)\n\nReturn the harmonic mean of a real-valued array.\n\n\n\n"
},

{
    "location": "means.html#StatsBase.genmean",
    "page": "Mean Functions",
    "title": "StatsBase.genmean",
    "category": "Function",
    "text": "genmean(a, p)\n\nReturn the generalized/power mean with exponent p of a real-valued array, i.e. left( frac1n sum_i=1^n a_i^p right)^frac1p, where n = length(a). It is taken to be the geometric mean when p == 0.\n\n\n\n"
},

{
    "location": "means.html#Base.mean-Tuple{AbstractArray,StatsBase.AbstractWeights}",
    "page": "Mean Functions",
    "title": "Base.mean",
    "category": "Method",
    "text": "mean(A::AbstractArray, w::AbstractWeights[, dim::Int])\n\nCompute the weighted mean of array A with weight vector w (of type AbstractWeights). If dim is provided, compute the weighted mean along dimension dim.\n\nExamples\n\nw = rand(n)\nmean(x, weights(w))\n\n\n\n"
},

{
    "location": "means.html#Base.mean!-Tuple{AbstractArray,AbstractArray,StatsBase.AbstractWeights,Int64}",
    "page": "Mean Functions",
    "title": "Base.mean!",
    "category": "Method",
    "text": "mean(R::AbstractArray, , A::AbstractArray, w::AbstractWeights[, dim::Int])\n\nCompute the weighted mean of array A with weight vector w (of type AbstractWeights) along dimension dim, and write results to R.\n\n\n\n"
},

{
    "location": "means.html#Mean-Functions-1",
    "page": "Mean Functions",
    "title": "Mean Functions",
    "category": "section",
    "text": "The package provides functions to compute means of different kinds.geomean\nharmmean\ngenmeanThe mean and mean! functions are also extended to accept a weight vector of type AbstractWeights to compute weighted mean.Base.mean(A::AbstractArray, w::AbstractWeights)\nBase.mean!(R::AbstractArray, A::AbstractArray, w::AbstractWeights, dim::Int)"
},

{
    "location": "scalarstats.html#",
    "page": "Scalar Statistics",
    "title": "Scalar Statistics",
    "category": "page",
    "text": ""
},

{
    "location": "scalarstats.html#Scalar-Statistics-1",
    "page": "Scalar Statistics",
    "title": "Scalar Statistics",
    "category": "section",
    "text": "The package implements functions for computing various statistics over an array of scalar real numbers."
},

{
    "location": "scalarstats.html#Base.var-Tuple{AbstractArray{T,N} where N where T<:Real,StatsBase.AbstractWeights}",
    "page": "Scalar Statistics",
    "title": "Base.var",
    "category": "Method",
    "text": "var(x, w::AbstractWeights, [dim]; mean=nothing, corrected=false)\n\nCompute the variance of a real-valued array x, optionally over a dimension dim. Observations in x are weighted using weight vector w. The uncorrected (when corrected=false) sample variance is defined as:\n\nfrac1sumw sum_i=1^n w_ileft(x_i - right)^2 \n\nwhere n is the length of the input and  is the mean. The unbiased estimate (when corrected=true) of the population variance is computed by replacing frac1sumw with a factor dependent on the type of weights used:\n\nAnalyticWeights: frac1sum w - sum w^2  sum w\nFrequencyWeights: frac1sumw - 1\nProbabilityWeights: fracn(n - 1) sum w where n equals count(!iszero, w)\nWeights: ArgumentError (bias correction not supported)\n\n\n\n"
},

{
    "location": "scalarstats.html#Base.std-Tuple{AbstractArray{T,N} where N where T<:Real,StatsBase.AbstractWeights}",
    "page": "Scalar Statistics",
    "title": "Base.std",
    "category": "Method",
    "text": "std(v, w::AbstractWeights, [dim]; mean=nothing, corrected=false)\n\nCompute the standard deviation of a real-valued array x, optionally over a dimension dim. Observations in x are weighted using weight vector w. The uncorrected (when corrected=false) sample standard deviation is defined as:\n\nsqrtfrac1sumw sum_i=1^n w_ileft(x_i - right)^2 \n\nwhere n is the length of the input and  is the mean. The unbiased estimate (when corrected=true) of the population standard deviation is computed by replacing frac1sumw with a factor dependent on the type of weights used:\n\nAnalyticWeights: frac1sum w - sum w^2  sum w\nFrequencyWeights: frac1sumw - 1\nProbabilityWeights: fracn(n - 1) sum w where n equals count(!iszero, w)\nWeights: ArgumentError (bias correction not supported)\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.mean_and_var",
    "page": "Scalar Statistics",
    "title": "StatsBase.mean_and_var",
    "category": "Function",
    "text": "mean_and_var(x, [w::AbstractWeights], [dim]; corrected=false) -> (mean, var)\n\nReturn the mean and variance of a real-valued array x, optionally over a dimension dim, as a tuple. Observations in x can be weighted using weight vector w. Finally, bias correction is be applied to the variance calculation if corrected=true. See var documentation for more details.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.mean_and_std",
    "page": "Scalar Statistics",
    "title": "StatsBase.mean_and_std",
    "category": "Function",
    "text": "mean_and_std(x, [w::AbstractWeights], [dim]; corrected=false) -> (mean, std)\n\nReturn the mean and standard deviation of a real-valued array x, optionally over a dimension dim, as a tuple. A weighting vector w can be specified to weight the estimates. Finally, bias correction is applied to the standard deviation calculation if corrected=true. See std documentation for more details.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.skewness",
    "page": "Scalar Statistics",
    "title": "StatsBase.skewness",
    "category": "Function",
    "text": "skewness(v, [wv::AbstractWeights], m=mean(v))\n\nCompute the standardized skewness of a real-valued array v, optionally specifying a weighting vector wv and a center m.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.kurtosis",
    "page": "Scalar Statistics",
    "title": "StatsBase.kurtosis",
    "category": "Function",
    "text": "kurtosis(v, [wv::AbstractWeights], m=mean(v))\n\nCompute the excess kurtosis of a real-valued array v, optionally specifying a weighting vector wv and a center m.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.moment",
    "page": "Scalar Statistics",
    "title": "StatsBase.moment",
    "category": "Function",
    "text": "moment(v, k, [wv::AbstractWeights], m=mean(v))\n\nReturn the kth order central moment of a real-valued array v, optionally specifying a weighting vector wv and a center m.\n\n\n\n"
},

{
    "location": "scalarstats.html#Moments-1",
    "page": "Scalar Statistics",
    "title": "Moments",
    "category": "section",
    "text": "Base.var(v::StatsBase.RealArray, w::AbstractWeights; mean=nothing, corrected::DepBool=nothing)\nBase.std(v::StatsBase.RealArray, w::AbstractWeights; mean=nothing, corrected::DepBool=nothing)\nmean_and_var\nmean_and_std\nskewness\nkurtosis\nmoment"
},

{
    "location": "scalarstats.html#StatsBase.span",
    "page": "Scalar Statistics",
    "title": "StatsBase.span",
    "category": "Function",
    "text": "span(x)\n\nReturn the span of an integer array, i.e. the range minimum(x):maximum(x). The minimum and maximum of x are computed in one-pass using extrema.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.variation",
    "page": "Scalar Statistics",
    "title": "StatsBase.variation",
    "category": "Function",
    "text": "variation(x, m=mean(x))\n\nReturn the coefficient of variation of an array x, optionally specifying a precomputed mean m. The coefficient of variation is the ratio of the standard deviation to the mean.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.sem",
    "page": "Scalar Statistics",
    "title": "StatsBase.sem",
    "category": "Function",
    "text": "sem(a)\n\nReturn the standard error of the mean of a, i.e. sqrt(var(a) / length(a)).\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.mad",
    "page": "Scalar Statistics",
    "title": "StatsBase.mad",
    "category": "Function",
    "text": "mad(v)\n\nCompute the median absolute deviation of v.\n\n\n\n"
},

{
    "location": "scalarstats.html#Measurements-of-Variation-1",
    "page": "Scalar Statistics",
    "title": "Measurements of Variation",
    "category": "section",
    "text": "span\nvariation\nsem\nmad"
},

{
    "location": "scalarstats.html#StatsBase.zscore",
    "page": "Scalar Statistics",
    "title": "StatsBase.zscore",
    "category": "Function",
    "text": "zscore(X, [μ, σ])\n\nCompute the z-scores of X, optionally specifying a precomputed mean μ and standard deviation σ. z-scores are the signed number of standard deviations above the mean that an observation lies, i.e. (x - )  .\n\nμ and σ should be both scalars or both arrays. The computation is broadcasting. In particular, when μ and σ are arrays, they should have the same size, and size(μ, i) == 1  || size(μ, i) == size(X, i) for each dimension.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.zscore!",
    "page": "Scalar Statistics",
    "title": "StatsBase.zscore!",
    "category": "Function",
    "text": "zscore!([Z], X, μ, σ)\n\nCompute the z-scores of an array X with mean μ and standard deviation σ. z-scores are the signed number of standard deviations above the mean that an observation lies, i.e. (x - )  .\n\nIf a destination array Z is provided, the scores are stored in Z and it must have the same shape as X. Otherwise X is overwritten.\n\n\n\n"
},

{
    "location": "scalarstats.html#Z-scores-1",
    "page": "Scalar Statistics",
    "title": "Z-scores",
    "category": "section",
    "text": "zscore\nzscore!"
},

{
    "location": "scalarstats.html#StatsBase.entropy",
    "page": "Scalar Statistics",
    "title": "StatsBase.entropy",
    "category": "Function",
    "text": "entropy(p, [b])\n\nCompute the entropy of an array p, optionally specifying a real number b such that the entropy is scaled by 1/log(b).\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.renyientropy",
    "page": "Scalar Statistics",
    "title": "StatsBase.renyientropy",
    "category": "Function",
    "text": "renyientropy(p, α)\n\nCompute the Rényi (generalized) entropy of order α of an array p.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.crossentropy",
    "page": "Scalar Statistics",
    "title": "StatsBase.crossentropy",
    "category": "Function",
    "text": "crossentropy(p, q, [b])\n\nCompute the cross entropy between p and q, optionally specifying a real number b such that the result is scaled by 1/log(b).\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.kldivergence",
    "page": "Scalar Statistics",
    "title": "StatsBase.kldivergence",
    "category": "Function",
    "text": "kldivergence(p, q, [b])\n\nCompute the Kullback-Leibler divergence of q from p, optionally specifying a real number b such that the divergence is scaled by 1/log(b).\n\n\n\n"
},

{
    "location": "scalarstats.html#Entropy-and-Related-Functions-1",
    "page": "Scalar Statistics",
    "title": "Entropy and Related Functions",
    "category": "section",
    "text": "entropy\nrenyientropy\ncrossentropy\nkldivergence"
},

{
    "location": "scalarstats.html#Quantile-and-Related-Functions-1",
    "page": "Scalar Statistics",
    "title": "Quantile and Related Functions",
    "category": "section",
    "text": "percentile\niqr\nnquantile\nquantile\nBase.median{W<:Real}(v::StatsBase.RealVector, w::AbstractWeights{W})"
},

{
    "location": "scalarstats.html#StatsBase.mode",
    "page": "Scalar Statistics",
    "title": "StatsBase.mode",
    "category": "Function",
    "text": "mode(a, [r])\n\nReturn the mode (most common number) of an array, optionally over a specified range r. If several modes exist, the first one (in order of appearance) is returned.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.modes",
    "page": "Scalar Statistics",
    "title": "StatsBase.modes",
    "category": "Function",
    "text": "modes(a, [r])::Vector\n\nReturn all modes (most common numbers) of an array, optionally over a specified range r.\n\n\n\n"
},

{
    "location": "scalarstats.html#Mode-and-Modes-1",
    "page": "Scalar Statistics",
    "title": "Mode and Modes",
    "category": "section",
    "text": "mode\nmodes"
},

{
    "location": "scalarstats.html#StatsBase.summarystats",
    "page": "Scalar Statistics",
    "title": "StatsBase.summarystats",
    "category": "Function",
    "text": "summarystats(a)\n\nCompute summary statistics for a real-valued array a. Returns a SummaryStats object containing the mean, minimum, 25th percentile, median, 75th percentile, and maxmimum.\n\n\n\n"
},

{
    "location": "scalarstats.html#StatsBase.describe",
    "page": "Scalar Statistics",
    "title": "StatsBase.describe",
    "category": "Function",
    "text": "describe(a)\n\nPretty-print the summary statistics provided by summarystats: the mean, minimum, 25th percentile, median, 75th percentile, and maximum.\n\n\n\n"
},

{
    "location": "scalarstats.html#Summary-Statistics-1",
    "page": "Scalar Statistics",
    "title": "Summary Statistics",
    "category": "section",
    "text": "summarystats\ndescribe"
},

{
    "location": "robust.html#",
    "page": "Robust Statistics",
    "title": "Robust Statistics",
    "category": "page",
    "text": ""
},

{
    "location": "robust.html#StatsBase.trim",
    "page": "Robust Statistics",
    "title": "StatsBase.trim",
    "category": "Function",
    "text": "trim(x; prop=0.0, count=0)\n\nReturn a copy of x with either count or proportion prop of the highest and lowest elements removed.  To compute the trimmed mean of x use mean(trim(x)); to compute the variance use trimvar(x) (see trimvar).\n\nExample\n\njulia> trim([1,2,3,4,5], prop=0.2)\n3-element Array{Int64,1}:\n 2\n 3\n 4\n\n\n\n"
},

{
    "location": "robust.html#StatsBase.winsor",
    "page": "Robust Statistics",
    "title": "StatsBase.winsor",
    "category": "Function",
    "text": "winsor(x; prop=0.0, count=0)\n\nReturn a copy of x with either count or proportion prop of the lowest elements of x replaced with the next-lowest, and an equal number of the highest elements replaced with the previous-highest.  To compute the Winsorized mean of x use mean(winsor(x)).\n\nExample\n\njulia> winsor([1,2,3,4,5], prop=0.2)\n5-element Array{Int64,1}:\n 2\n 2\n 3\n 4\n 4\n\n\n\n"
},

{
    "location": "robust.html#StatsBase.trimvar",
    "page": "Robust Statistics",
    "title": "StatsBase.trimvar",
    "category": "Function",
    "text": "trimvar(x; prop=0.0, count=0)\n\nCompute the variance of the trimmed mean of x. This function uses the Winsorized variance, as described in Wilcox (2010).\n\n\n\n"
},

{
    "location": "robust.html#Robust-Statistics-1",
    "page": "Robust Statistics",
    "title": "Robust Statistics",
    "category": "section",
    "text": "trim\nwinsor\ntrimvar"
},

{
    "location": "deviation.html#",
    "page": "Computing Deviations",
    "title": "Computing Deviations",
    "category": "page",
    "text": ""
},

{
    "location": "deviation.html#StatsBase.counteq",
    "page": "Computing Deviations",
    "title": "StatsBase.counteq",
    "category": "Function",
    "text": "counteq(a, b)\n\nCount the number of indices at which the elements of the arrays a and b are equal.\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.countne",
    "page": "Computing Deviations",
    "title": "StatsBase.countne",
    "category": "Function",
    "text": "countne(a, b)\n\nCount the number of indices at which the elements of the arrays a and b are not equal.\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.sqL2dist",
    "page": "Computing Deviations",
    "title": "StatsBase.sqL2dist",
    "category": "Function",
    "text": "sqL2dist(a, b)\n\nCompute the squared L2 distance between two arrays: sum_i=1^n a_i - b_i^2. Efficient equivalent of sumabs2(a - b).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.L2dist",
    "page": "Computing Deviations",
    "title": "StatsBase.L2dist",
    "category": "Function",
    "text": "L2dist(a, b)\n\nCompute the L2 distance between two arrays: sqrtsum_i=1^n a_i - b_i^2. Efficient equivalent of sqrt(sumabs2(a - b)).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.L1dist",
    "page": "Computing Deviations",
    "title": "StatsBase.L1dist",
    "category": "Function",
    "text": "L1dist(a, b)\n\nCompute the L1 distance between two arrays: sum_i=1^n a_i - b_i. Efficient equivalent of sum(abs, a - b).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.Linfdist",
    "page": "Computing Deviations",
    "title": "StatsBase.Linfdist",
    "category": "Function",
    "text": "Linfdist(a, b)\n\nCompute the L∞ distance, also called the Chebyshev distance, between two arrays: max_iin1n a_i - b_i. Efficient equivalent of maxabs(a - b).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.gkldiv",
    "page": "Computing Deviations",
    "title": "StatsBase.gkldiv",
    "category": "Function",
    "text": "gkldiv(a, b)\n\nCompute the generalized Kullback-Leibler divergence between two arrays: sum_i=1^n (a_i log(a_ib_i) - a_i + b_i). Efficient equivalent of sum(a*log(a/b)-a+b).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.meanad",
    "page": "Computing Deviations",
    "title": "StatsBase.meanad",
    "category": "Function",
    "text": "meanad(a, b)\n\nReturn the mean absolute deviation between two arrays: mean(abs(a - b)).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.maxad",
    "page": "Computing Deviations",
    "title": "StatsBase.maxad",
    "category": "Function",
    "text": "maxad(a, b)\n\nReturn the maximum absolute deviation between two arrays: maxabs(a - b).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.msd",
    "page": "Computing Deviations",
    "title": "StatsBase.msd",
    "category": "Function",
    "text": "msd(a, b)\n\nReturn the mean squared deviation between two arrays: mean(abs2(a - b)).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.rmsd",
    "page": "Computing Deviations",
    "title": "StatsBase.rmsd",
    "category": "Function",
    "text": "rmsd(a, b; normalize=false)\n\nReturn the root mean squared deviation between two optionally normalized arrays. The root mean squared deviation is computed as sqrt(msd(a, b)).\n\n\n\n"
},

{
    "location": "deviation.html#StatsBase.psnr",
    "page": "Computing Deviations",
    "title": "StatsBase.psnr",
    "category": "Function",
    "text": "psnr(a, b, maxv)\n\nCompute the peak signal-to-noise ratio between two arrays a and b. maxv is the maximum possible value either array can take. The PSNR is computed as 10 * log10(maxv^2 / msd(a, b)).\n\n\n\n"
},

{
    "location": "deviation.html#Computing-Deviations-1",
    "page": "Computing Deviations",
    "title": "Computing Deviations",
    "category": "section",
    "text": "This package provides functions to compute various deviations between arrays in a variety of ways:counteq\ncountne\nsqL2dist\nL2dist\nL1dist\nLinfdist\ngkldiv\nmeanad\nmaxad\nmsd\nrmsd\npsnrnote: Note\nAll these functions are implemented in a reasonably efficient way without creating any temporary arrays in the middle."
},

{
    "location": "cov.html#",
    "page": "Scatter Matrix and Covariance",
    "title": "Scatter Matrix and Covariance",
    "category": "page",
    "text": ""
},

{
    "location": "cov.html#StatsBase.scattermat",
    "page": "Scatter Matrix and Covariance",
    "title": "StatsBase.scattermat",
    "category": "Function",
    "text": "scattermat(X, [wv::AbstractWeights]; mean=nothing, vardim=1)\n\nCompute the scatter matrix, which is an unnormalized covariance matrix. A weighting vector wv can be specified to weight the estimate.\n\nArguments\n\nmean=nothing: a known mean value. nothing indicates that the mean is unknown, and the function will compute the mean. Specifying mean=0 indicates that the data are centered and hence there's no need to subtract the mean.\nvardim=1: the dimension along which the variables are organized. When vardim = 1, the variables are considered columns with observations in rows; when vardim = 2, variables are in rows with observations in columns.\n\n\n\n"
},

{
    "location": "cov.html#Base.cov",
    "page": "Scatter Matrix and Covariance",
    "title": "Base.cov",
    "category": "Function",
    "text": "cov(X, w::AbstractWeights; mean=nothing, vardim=1, corrected=false)\n\nCompute the weighted covariance matrix. Similar to var and std the biased covariance matrix (corrected=false) is computed by multiplying scattermat(X, w) by frac1sumw to normalize. However, the unbiased covariance matrix (corrected=true) is dependent on the type of weights used:\n\nAnalyticWeights: frac1sum w - sum w^2  sum w\nFrequencyWeights: frac1sumw - 1\nProbabilityWeights: fracn(n - 1) sum w where n equals count(!iszero, w)\nWeights: ArgumentError (bias correction not supported)\n\n\n\n"
},

{
    "location": "cov.html#Base.cor",
    "page": "Scatter Matrix and Covariance",
    "title": "Base.cor",
    "category": "Function",
    "text": "cor(X, w::AbstractWeights, vardim=1)\n\nCompute the Pearson correlation matrix of X along the dimension vardim with a weighting w .\n\n\n\n"
},

{
    "location": "cov.html#StatsBase.mean_and_cov",
    "page": "Scatter Matrix and Covariance",
    "title": "StatsBase.mean_and_cov",
    "category": "Function",
    "text": "mean_and_cov(x, [wv::AbstractWeights]; vardim=1, corrected=false) -> (mean, cov)\n\nReturn the mean and covariance matrix as a tuple. A weighting vector wv can be specified. vardim that designates whether the variables are columns in the matrix (1) or rows (2). Finally, bias correction is applied to the covariance calculation if corrected=true. See cov documentation for more details.\n\n\n\n"
},

{
    "location": "cov.html#StatsBase.cov2cor",
    "page": "Scatter Matrix and Covariance",
    "title": "StatsBase.cov2cor",
    "category": "Function",
    "text": "cov2cor(C, s)\n\nCompute the correlation matrix from the covariance matrix C and a vector of standard deviations s. Use Base.cov2cor! for an in-place version.\n\n\n\n"
},

{
    "location": "cov.html#StatsBase.cor2cov",
    "page": "Scatter Matrix and Covariance",
    "title": "StatsBase.cor2cov",
    "category": "Function",
    "text": "cor2cov(C, s)\n\nCompute the covariance matrix from the correlation matrix C and a vector of standard deviations s. Use StatsBase.cor2cov! for an in-place version.\n\n\n\n"
},

{
    "location": "cov.html#Scatter-Matrix-and-Covariance-1",
    "page": "Scatter Matrix and Covariance",
    "title": "Scatter Matrix and Covariance",
    "category": "section",
    "text": "This package implements functions for computing scatter matrix, as well as weighted covariance matrix.scattermat\ncov\ncor\nmean_and_cov\ncov2cor\ncor2cov"
},

{
    "location": "counts.html#",
    "page": "Counting Functions",
    "title": "Counting Functions",
    "category": "page",
    "text": ""
},

{
    "location": "counts.html#Counting-Functions-1",
    "page": "Counting Functions",
    "title": "Counting Functions",
    "category": "section",
    "text": "The package provides functions to count the occurences of distinct values."
},

{
    "location": "counts.html#StatsBase.counts",
    "page": "Counting Functions",
    "title": "StatsBase.counts",
    "category": "Function",
    "text": "counts(x, [wv::AbstractWeights])\ncounts(x, levels::UnitRange{<:Integer}, [wv::AbstractWeights])\ncounts(x, k::Integer, [wv::AbstractWeights])\n\nCount the number of times each value in x occurs. If levels is provided, only values falling in that range will be considered (the others will be ignored without raising an error or a warning). If an integer k is provided, only values in the range 1:k will be considered.\n\nIf a weighting vector wv is specified, the sum of the weights is used rather than the raw counts.\n\nThe output is a vector of length length(levels).\n\n\n\n"
},

{
    "location": "counts.html#StatsBase.proportions",
    "page": "Counting Functions",
    "title": "StatsBase.proportions",
    "category": "Function",
    "text": "proportions(x, levels=span(x), [wv::AbstractWeights])\n\nReturn the proportion of values in the range levels that occur in x. Equivalent to counts(x, levels) / length(x). If a weighting vector wv is specified, the sum of the weights is used rather than the raw counts.\n\n\n\nproportions(x, k::Integer, [wv::AbstractWeights])\n\nReturn the proportion of integers in 1 to k that occur in x.\n\n\n\n"
},

{
    "location": "counts.html#StatsBase.addcounts!-Tuple{AbstractArray,AbstractArray{T,N} where N where T<:Integer,UnitRange{T} where T<:Integer}",
    "page": "Counting Functions",
    "title": "StatsBase.addcounts!",
    "category": "Method",
    "text": "addcounts!(r, x, levels::UnitRange{<:Int}, [wv::AbstractWeights])\n\nAdd the number of occurrences in x of each value in levels to an existing array r. If a weighting vector wv is specified, the sum of weights is used rather than the raw counts.\n\n\n\n"
},

{
    "location": "counts.html#Counting-over-an-Integer-Range-1",
    "page": "Counting Functions",
    "title": "Counting over an Integer Range",
    "category": "section",
    "text": "counts\nproportions\naddcounts!(r::AbstractArray, x::StatsBase.IntegerArray, levels::StatsBase.IntUnitRange)"
},

{
    "location": "counts.html#StatsBase.countmap",
    "page": "Counting Functions",
    "title": "StatsBase.countmap",
    "category": "Function",
    "text": "countmap(x; alg = :auto)\n\nReturn a dictionary mapping each unique value in x to its number of occurrences.\n\n:auto (default): if StatsBase.radixsort_safe(eltype(x)) == true then use                    :radixsort, otherwise use :dict.\n:radixsort:      if radixsort_safe(eltype(x)) == true then use the                    radix sort                    algorithm to sort the input vector which will generally lead to                    shorter running time. However the radix sort algorithm creates a                    copy of the input vector and hence uses more RAM. Choose :dict                    if the amount of available RAM is a limitation.\n:dict:           use Dict-based method which is generally slower but uses less                    RAM and is safe for any data type.\n\n\n\n"
},

{
    "location": "counts.html#StatsBase.proportionmap",
    "page": "Counting Functions",
    "title": "StatsBase.proportionmap",
    "category": "Function",
    "text": "proportionmap(x)\n\nReturn a dictionary mapping each unique value in x to its proportion in x.\n\n\n\n"
},

{
    "location": "counts.html#StatsBase.addcounts!-Union{Tuple{Dict{T,V} where V,AbstractArray{T,N} where N}, Tuple{T}} where T",
    "page": "Counting Functions",
    "title": "StatsBase.addcounts!",
    "category": "Method",
    "text": "addcounts!(dict, x[, wv]; alg = :auto)\n\nAdd counts based on x to a count map. New entries will be added if new values come up. If a weighting vector wv is specified, the sum of the weights is used rather than the raw counts.\n\nalg can be one of:\n\n:auto (default): if StatsBase.radixsort_safe(eltype(x)) == true then use                    :radixsort, otherwise use :dict.\n:radixsort:      if radixsort_safe(eltype(x)) == true then use the                    radix sort                    algorithm to sort the input vector which will generally lead to                    shorter running time. However the radix sort algorithm creates a                    copy of the input vector and hence uses more RAM. Choose :dict                    if the amount of available RAM is a limitation.\n:dict:           use Dict-based method which is generally slower but uses less                    RAM and is safe for any data type.\n\n\n\n"
},

{
    "location": "counts.html#Counting-over-arbitrary-distinct-values-1",
    "page": "Counting Functions",
    "title": "Counting over arbitrary distinct values",
    "category": "section",
    "text": "countmap\nproportionmap\naddcounts!{T}(cm::Dict{T}, x::AbstractArray{T})"
},

{
    "location": "ranking.html#",
    "page": "Rankings and Rank Correlations",
    "title": "Rankings and Rank Correlations",
    "category": "page",
    "text": ""
},

{
    "location": "ranking.html#StatsBase.ordinalrank",
    "page": "Rankings and Rank Correlations",
    "title": "StatsBase.ordinalrank",
    "category": "Function",
    "text": "ordinalrank(x; lt = isless, rev::Bool = false)\n\nReturn the ordinal ranking (\"1234\" ranking) of an array. The lt keyword allows providing a custom \"less than\" function; use rev=true to reverse the sorting order. All items in x are given distinct, successive ranks based on their position in sort(x; lt = lt, rev = rev). Missing values are assigned rank missing.\n\n\n\n"
},

{
    "location": "ranking.html#StatsBase.competerank",
    "page": "Rankings and Rank Correlations",
    "title": "StatsBase.competerank",
    "category": "Function",
    "text": "competerank(x; lt = isless, rev::Bool = false)\n\nReturn the standard competition ranking (\"1224\" ranking) of an array. The lt keyword allows providing a custom \"less than\" function; use rev=true to reverse the sorting order. Items that compare equal are given the same rank, then a gap is left in the rankings the size of the number of tied items - 1. Missing values are assigned rank missing.\n\n\n\n"
},

{
    "location": "ranking.html#StatsBase.denserank",
    "page": "Rankings and Rank Correlations",
    "title": "StatsBase.denserank",
    "category": "Function",
    "text": "denserank(x)\n\nReturn the dense ranking (\"1223\" ranking) of an array. The lt keyword allows providing a custom \"less than\" function; use rev=true to reverse the sorting order. Items that compare equal receive the same ranking, and the next subsequent rank is assigned with no gap. Missing values are assigned rank missing.\n\n\n\n"
},

{
    "location": "ranking.html#StatsBase.tiedrank",
    "page": "Rankings and Rank Correlations",
    "title": "StatsBase.tiedrank",
    "category": "Function",
    "text": "tiedrank(x)\n\nReturn the tied ranking, also called fractional or \"1 2.5 2.5 4\" ranking, of an array. The lt keyword allows providing a custom \"less than\" function; use rev=true to reverse the sorting order. Items that compare equal receive the mean of the rankings they would have been assigned under ordinal ranking. Missing values are assigned rank missing.\n\n\n\n"
},

{
    "location": "ranking.html#StatsBase.corspearman",
    "page": "Rankings and Rank Correlations",
    "title": "StatsBase.corspearman",
    "category": "Function",
    "text": "corspearman(x, y=x)\n\nCompute Spearman's rank correlation coefficient. If x and y are vectors, the output is a float, otherwise it's a matrix corresponding to the pairwise correlations of the columns of x and y.\n\n\n\n"
},

{
    "location": "ranking.html#StatsBase.corkendall",
    "page": "Rankings and Rank Correlations",
    "title": "StatsBase.corkendall",
    "category": "Function",
    "text": "corkendall(x, y=x)\n\nCompute Kendall's rank correlation coefficient, τ. x and y must both be either matrices or vectors.\n\n\n\n"
},

{
    "location": "ranking.html#Rankings-and-Rank-Correlations-1",
    "page": "Rankings and Rank Correlations",
    "title": "Rankings and Rank Correlations",
    "category": "section",
    "text": "This package implements various strategies for computing ranks and rank correlations.ordinalrank\ncompeterank\ndenserank\ntiedrank\ncorspearman\ncorkendall"
},

{
    "location": "sampling.html#",
    "page": "Sampling from Population",
    "title": "Sampling from Population",
    "category": "page",
    "text": ""
},

{
    "location": "sampling.html#Sampling-from-Population-1",
    "page": "Sampling from Population",
    "title": "Sampling from Population",
    "category": "section",
    "text": ""
},

{
    "location": "sampling.html#StatsBase.sample",
    "page": "Sampling from Population",
    "title": "StatsBase.sample",
    "category": "Function",
    "text": "sample([rng], a, [wv::AbstractWeights])\n\nSelect a single random element of a. Sampling probabilities are proportional to the weights given in wv, if provided.\n\nOptionally specify a random number generator rng as the first argument (defaults to Base.GLOBAL_RNG).\n\n\n\nsample([rng], a, [wv::AbstractWeights], n::Integer; replace=true, ordered=false)\n\nSelect a random, optionally weighted sample of size n from an array a using a polyalgorithm. Sampling probabilities are proportional to the weights given in wv, if provided. replace dictates whether sampling is performed with replacement and order dictates whether an ordered sample, also called a sequential sample, should be taken.\n\nOptionally specify a random number generator rng as the first argument (defaults to Base.GLOBAL_RNG).\n\n\n\nsample([rng], a, [wv::AbstractWeights], dims::Dims; replace=true, ordered=false)\n\nSelect a random, optionally weighted sample from an array a specifying the dimensions dims of the output array. Sampling probabilities are proportional to the weights given in wv, if provided. replace dictates whether sampling is performed with replacement and order dictates whether an ordered sample, also called a sequential sample, should be taken.\n\nOptionally specify a random number generator rng as the first argument (defaults to Base.GLOBAL_RNG).\n\n\n\nsample([rng], wv::AbstractWeights)\n\nSelect a single random integer in 1:length(wv) with probabilities proportional to the weights given in wv.\n\nOptionally specify a random number generator rng as the first argument (defaults to Base.GLOBAL_RNG).\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.sample!",
    "page": "Sampling from Population",
    "title": "StatsBase.sample!",
    "category": "Function",
    "text": "sample!([rng], a, [wv::AbstractWeights], x; replace=true, ordered=false)\n\nDraw a random sample of length(x) elements from an array a and store the result in x. A polyalgorithm is used for sampling. Sampling probabilities are proportional to the weights given in wv, if provided. replace dictates whether sampling is performed with replacement and order dictates whether an ordered sample, also called a sequential sample, should be taken.\n\nOptionally specify a random number generator rng as the first argument (defaults to Base.GLOBAL_RNG).\n\n\n\n"
},

{
    "location": "sampling.html#Sampling-API-1",
    "page": "Sampling from Population",
    "title": "Sampling API",
    "category": "section",
    "text": "The package provides functions for sampling from a given population (with or without replacement).sample\nsample!"
},

{
    "location": "sampling.html#Algorithms-1",
    "page": "Sampling from Population",
    "title": "Algorithms",
    "category": "section",
    "text": "Internally, this package implements multiple algorithms, and the sample (and sample!) methods integrate them into a poly-algorithm, which chooses a specific algorithm based on inputs.Note that the choices made in sample are decided based on extensive benchmarking (see perf/sampling.jl and perf/wsampling.jl). It performs reasonably fast for most cases. That being said, if you know that a certain algorithm is particularly suitable for your context, directly calling an internal algorithm function might be slightly more efficient.Here are a list of algorithms implemented in the package. The functions below are not exported (one can still import them from StatsBase via using though)."
},

{
    "location": "sampling.html#Notations-1",
    "page": "Sampling from Population",
    "title": "Notations",
    "category": "section",
    "text": "a: source array representing the population\nx: the destination array\nwv: the weight vector (of type AbstractWeights), for weighted sampling\nn: the length of a\nk: the length of x. For sampling without replacement, k must not exceed n.\nrng: optional random number generator (defaults to Base.GLOBAL_RNG)All following functions write results to x (pre-allocated) and return x."
},

{
    "location": "sampling.html#StatsBase.direct_sample!-Tuple{AbstractRNG,AbstractArray,AbstractArray}",
    "page": "Sampling from Population",
    "title": "StatsBase.direct_sample!",
    "category": "Method",
    "text": "direct_sample!([rng], a::AbstractArray, x::AbstractArray)\n\nDirect sampling: for each j in 1:k, randomly pick i from 1:n, and set x[j] = a[i], with n=length(a) and k=length(x).\n\nThis algorithm consumes k random numbers.\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.samplepair",
    "page": "Sampling from Population",
    "title": "StatsBase.samplepair",
    "category": "Function",
    "text": "samplepair([rng], n)\n\nDraw a pair of distinct integers between 1 and n without replacement.\n\nOptionally specify a random number generator rng as the first argument (defaults to Base.GLOBAL_RNG).\n\n\n\nsamplepair([rng], a)\n\nDraw a pair of distinct elements from the array a without replacement.\n\nOptionally specify a random number generator rng as the first argument (defaults to Base.GLOBAL_RNG).\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.knuths_sample!",
    "page": "Sampling from Population",
    "title": "StatsBase.knuths_sample!",
    "category": "Function",
    "text": "knuths_sample!([rng], a, x)\n\nKnuth's Algorithm S for random sampling without replacement.\n\nReference: D. Knuth. The Art of Computer Programming. Vol 2, 3.4.2, p.142.\n\nThis algorithm consumes length(a) random numbers. It requires no additional memory space. Suitable for the case where memory is tight.\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.fisher_yates_sample!",
    "page": "Sampling from Population",
    "title": "StatsBase.fisher_yates_sample!",
    "category": "Function",
    "text": "fisher_yates_sample!([rng], a::AbstractArray, x::AbstractArray)\n\nFisher-Yates shuffling (with early termination).\n\nPseudo-code:\n\nn = length(a)\nk = length(x)\n\n# Create an array of the indices\ninds = collect(1:n)\n\nfor i = 1:k\n    # swap element `i` with another random element in inds[i:n]\n    # set element `i` in `x`\nend\n\nThis algorithm consumes k=length(x) random numbers. It uses an integer array of length n=length(a) internally to maintain the shuffled indices. It is considerably faster than Knuth's algorithm especially when n is greater than k. It is O(n) for initialization, plus O(k) for random shuffling\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.self_avoid_sample!",
    "page": "Sampling from Population",
    "title": "StatsBase.self_avoid_sample!",
    "category": "Function",
    "text": "self_avoid_sample!([rng], a::AbstractArray, x::AbstractArray)\n\nSelf-avoid sampling: use a set to maintain the index that has been sampled. Each time draw a new index, if the index has already been sampled, redraw until it draws an unsampled one.\n\nThis algorithm consumes about (or slightly more than) k=length(x) random numbers, and requires O(k) memory to store the set of sampled indices. Very fast when n  k, with n=length(a).\n\nHowever, if k is large and approaches n, the rejection rate would increase drastically, resulting in poorer performance.\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.seqsample_a!",
    "page": "Sampling from Population",
    "title": "StatsBase.seqsample_a!",
    "category": "Function",
    "text": "seqsample_a!([rng], a::AbstractArray, x::AbstractArray)\n\nRandom subsequence sampling using algorithm A described in the following paper (page 714): Jeffrey Scott Vitter. \"Faster Methods for Random Sampling\". Communications of the ACM, 27 (7), July 1984.\n\nThis algorithm consumes O(n) random numbers, with n=length(a). The outputs are ordered.\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.seqsample_c!",
    "page": "Sampling from Population",
    "title": "StatsBase.seqsample_c!",
    "category": "Function",
    "text": "seqsample_c!([rng], a::AbstractArray, x::AbstractArray)\n\nRandom subsequence sampling using algorithm C described in the following paper (page 715): Jeffrey Scott Vitter. \"Faster Methods for Random Sampling\". Communications of the ACM, 27 (7), July 1984.\n\nThis algorithm consumes O(k^2) random numbers, with k=length(x). The outputs are ordered.\n\n\n\n"
},

{
    "location": "sampling.html#Sampling-Algorithms-(Non-Weighted)-1",
    "page": "Sampling from Population",
    "title": "Sampling Algorithms (Non-Weighted)",
    "category": "section",
    "text": "StatsBase.direct_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)\nsamplepair\nStatsBase.knuths_sample!\nStatsBase.fisher_yates_sample!\nStatsBase.self_avoid_sample!\nStatsBase.seqsample_a!\nStatsBase.seqsample_c!"
},

{
    "location": "sampling.html#StatsBase.direct_sample!-Tuple{AbstractRNG,AbstractArray,StatsBase.AbstractWeights,AbstractArray}",
    "page": "Sampling from Population",
    "title": "StatsBase.direct_sample!",
    "category": "Method",
    "text": "direct_sample!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)\n\nDirect sampling.\n\nDraw each sample by scanning the weight vector.\n\nNoting k=length(x) and n=length(a), this algorithm:\n\nconsumes k random numbers\nhas time complexity O(n k), as scanning the weight vector each time takes O(n)\nrequires no additional memory space.\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.alias_sample!",
    "page": "Sampling from Population",
    "title": "StatsBase.alias_sample!",
    "category": "Function",
    "text": "alias_sample!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)\n\nAlias method.\n\nBuild an alias table, and sample therefrom.\n\nReference: Walker, A. J. \"An Efficient Method for Generating Discrete Random Variables with General Distributions.\" ACM Transactions on Mathematical Software 3 (3): 253, 1977.\n\nNoting k=length(x) and n=length(a), this algorithm takes O(n log n) time for building the alias table, and then O(1) to draw each sample. It consumes 2 k random numbers.\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.naive_wsample_norep!",
    "page": "Sampling from Population",
    "title": "StatsBase.naive_wsample_norep!",
    "category": "Function",
    "text": "naive_wsample_norep!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)\n\nNaive implementation of weighted sampling without replacement.\n\nIt makes a copy of the weight vector at initialization, and sets the weight to zero when the corresponding sample is picked.\n\nNoting k=length(x) and n=length(a), this algorithm consumes O(k) random numbers, and has overall time complexity O(n k).\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.efraimidis_a_wsample_norep!",
    "page": "Sampling from Population",
    "title": "StatsBase.efraimidis_a_wsample_norep!",
    "category": "Function",
    "text": "efraimidis_a_wsample_norep!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)\n\nWeighted sampling without replacement using Efraimidis-Spirakis A algorithm.\n\nReference: Efraimidis, P. S., Spirakis, P. G. \"Weighted random sampling with a reservoir.\" Information Processing Letters, 97 (5), 181-185, 2006. doi:10.1016/j.ipl.2005.11.003.\n\nNoting k=length(x) and n=length(a), this algorithm takes O(n + k log k) processing time to draw k elements. It consumes n random numbers.\n\n\n\n"
},

{
    "location": "sampling.html#StatsBase.efraimidis_ares_wsample_norep!",
    "page": "Sampling from Population",
    "title": "StatsBase.efraimidis_ares_wsample_norep!",
    "category": "Function",
    "text": "efraimidis_ares_wsample_norep!([rng], a::AbstractArray, wv::AbstractWeights, x::AbstractArray)\n\nImplementation of weighted sampling without replacement using Efraimidis-Spirakis A-Res algorithm.\n\nReference: Efraimidis, P. S., Spirakis, P. G. \"Weighted random sampling with a reservoir.\" Information Processing Letters, 97 (5), 181-185, 2006. doi:10.1016/j.ipl.2005.11.003.\n\nNoting k=length(x) and n=length(a), this algorithm takes O(k log(k) log(n  k)) processing time to draw k elements. It consumes n random numbers.\n\n\n\n"
},

{
    "location": "sampling.html#Weighted-Sampling-Algorithms-1",
    "page": "Sampling from Population",
    "title": "Weighted Sampling Algorithms",
    "category": "section",
    "text": "StatsBase.direct_sample!(rng::AbstractRNG, a::AbstractArray, wv::AbstractWeights, x::AbstractArray)\nStatsBase.alias_sample!\nStatsBase.naive_wsample_norep!\nStatsBase.efraimidis_a_wsample_norep!\nStatsBase.efraimidis_ares_wsample_norep!"
},

{
    "location": "empirical.html#",
    "page": "Empirical Estimation",
    "title": "Empirical Estimation",
    "category": "page",
    "text": ""
},

{
    "location": "empirical.html#Empirical-Estimation-1",
    "page": "Empirical Estimation",
    "title": "Empirical Estimation",
    "category": "section",
    "text": ""
},

{
    "location": "empirical.html#StatsBase.fit-Tuple{Type{StatsBase.Histogram},Vararg{Any,N} where N}",
    "page": "Empirical Estimation",
    "title": "StatsBase.fit",
    "category": "Method",
    "text": "fit(Histogram, data[, weight][, edges]; closed=:right, nbins)\n\nFit a histogram to data.\n\nArguments\n\ndata: either a vector (for a 1-dimensional histogram), or a tuple of vectors of equal length (for an n-dimensional histogram).\nweight: an optional AbstractWeights (of the same length as the data vectors), denoting the weight each observation contributes to the bin. If no weight vector is supplied, each observation has weight 1.\nedges: a vector (typically an AbstractRange object), or tuple of vectors, that gives the edges of the bins along each dimension. If no edges are provided, these are determined from the data.\n\nKeyword arguments\n\nclosed=:right: if :left, the bin intervals are left-closed [a,b); if :right (the default), intervals are right-closed (a,b].\nnbins: if no edges argument is supplied, the approximate number of bins to use along each dimension (can be either a single integer, or a tuple of integers).\n\nExamples\n\n# Univariate\nh = fit(Histogram, rand(100))\nh = fit(Histogram, rand(100), 0:0.1:1.0)\nh = fit(Histogram, rand(100), nbins=10)\nh = fit(Histogram, rand(100), weights(rand(100)), 0:0.1:1.0)\nh = fit(Histogram, [20], 0:20:100)\nh = fit(Histogram, [20], 0:20:100, closed=:left)\n\n# Multivariate\nh = fit(Histogram, (rand(100),rand(100)))\nh = fit(Histogram, (rand(100),rand(100)),nbins=10)\n\n\n\n"
},

{
    "location": "empirical.html#Histograms-1",
    "page": "Empirical Estimation",
    "title": "Histograms",
    "category": "section",
    "text": "The Histogram type represents data that has been tabulated into intervals (known as bins) along the real line, or in higher dimensions, over the real plane.Histograms can be fitted to data using the fit method.fit(::Type{Histogram}, args...; kwargs...)"
},

{
    "location": "empirical.html#StatsBase.ecdf",
    "page": "Empirical Estimation",
    "title": "StatsBase.ecdf",
    "category": "Function",
    "text": "ecdf(X)\n\nReturn an empirical cumulative distribution function (ECDF) based on a vector of samples given in X.\n\nNote: this is a higher-level function that returns a function, which can then be applied to evaluate CDF values on other samples.\n\n\n\n"
},

{
    "location": "empirical.html#Empirical-Cumulative-Distribution-Function-1",
    "page": "Empirical Estimation",
    "title": "Empirical Cumulative Distribution Function",
    "category": "section",
    "text": "ecdf"
},

{
    "location": "signalcorr.html#",
    "page": "Correlation Analysis of Signals",
    "title": "Correlation Analysis of Signals",
    "category": "page",
    "text": ""
},

{
    "location": "signalcorr.html#Correlation-Analysis-of-Signals-1",
    "page": "Correlation Analysis of Signals",
    "title": "Correlation Analysis of Signals",
    "category": "section",
    "text": "The package provides functions to perform correlation analysis of sequential signals. "
},

{
    "location": "signalcorr.html#StatsBase.autocov",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.autocov",
    "category": "Function",
    "text": "autocov(x, [lags]; demean=true)\n\nCompute the autocovariance of a vector or matrix x, optionally specifying the lags at which to compute the autocovariance. demean denotes whether the mean of x should be subtracted from x before computing the autocovariance.\n\nIf x is a vector, return a vector of the same length as x. If x is a matrix, return a matrix of size (length(lags), size(x,2)), where each column in the result corresponds to a column in x.\n\nWhen left unspecified, the lags used are the integers from 0 to min(size(x,1)-1, 10*log10(size(x,1))).\n\nThe output is not normalized. See autocor for a function with normalization.\n\n\n\n"
},

{
    "location": "signalcorr.html#StatsBase.autocov!",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.autocov!",
    "category": "Function",
    "text": "autocov!(r, x, lags; demean=true)\n\nCompute the autocovariance of a vector or matrix x at lags and store the result in r. demean denotes whether the mean of x should be subtracted from x before computing the autocovariance.\n\nIf x is a vector, r must be a vector of the same length as x. If x is a matrix, r must be a matrix of size (length(lags), size(x,2)), and where each column in the result will correspond to a column in x.\n\nThe output is not normalized. See autocor! for a method with normalization.\n\n\n\n"
},

{
    "location": "signalcorr.html#StatsBase.autocor",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.autocor",
    "category": "Function",
    "text": "autocor(x, [lags]; demean=true)\n\nCompute the autocorrelation function (ACF) of a vector or matrix x, optionally specifying the lags. demean denotes whether the mean of x should be subtracted from x before computing the ACF.\n\nIf x is a vector, return a vector of the same length as x. If x is a matrix, return a matrix of size (length(lags), size(x,2)), where each column in the result corresponds to a column in x.\n\nWhen left unspecified, the lags used are the integers from 0 to min(size(x,1)-1, 10*log10(size(x,1))).\n\nThe output is normalized by the variance of x, i.e. so that the lag 0 autocorrelation is 1. See autocov for the unnormalized form.\n\n\n\n"
},

{
    "location": "signalcorr.html#StatsBase.autocor!",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.autocor!",
    "category": "Function",
    "text": "autocor!(r, x, lags; demean=true)\n\nCompute the autocorrelation function (ACF) of a vector or matrix x at lags and store the result in r. demean denotes whether the mean of x should be subtracted from x before computing the ACF.\n\nIf x is a vector, r must be a vector of the same length as x. If x is a matrix, r must be a matrix of size (length(lags), size(x,2)), and where each column in the result will correspond to a column in x.\n\nThe output is normalized by the variance of x, i.e. so that the lag 0 autocorrelation is 1. See autocov! for the unnormalized form.\n\n\n\n"
},

{
    "location": "signalcorr.html#Autocovariance-and-Autocorrelation-1",
    "page": "Correlation Analysis of Signals",
    "title": "Autocovariance and Autocorrelation",
    "category": "section",
    "text": "autocov\nautocov!\nautocor\nautocor!"
},

{
    "location": "signalcorr.html#StatsBase.crosscov",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.crosscov",
    "category": "Function",
    "text": "crosscov(x, y, [lags]; demean=true)\n\nCompute the cross covariance function (CCF) between real-valued vectors or matrices x and y, optionally specifying the lags. demean specifies whether the respective means of x and y should be subtracted from them before computing their CCF.\n\nIf both x and y are vectors, return a vector of the same length as lags. Otherwise, compute cross covariances between each pairs of columns in x and y.\n\nWhen left unspecified, the lags used are the integers from -min(size(x,1)-1, 10*log10(size(x,1))) to min(size(x,1), 10*log10(size(x,1))).\n\nThe output is not normalized. See crosscor for a function with normalization.\n\n\n\n"
},

{
    "location": "signalcorr.html#StatsBase.crosscov!",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.crosscov!",
    "category": "Function",
    "text": "crosscov!(r, x, y, lags; demean=true)\n\nCompute the cross covariance function (CCF) between real-valued vectors or matrices x and y at lags and store the result in r. demean specifies whether the respective means of x and y should be subtracted from them before computing their CCF.\n\nIf both x and y are vectors, r must be a vector of the same length as lags. If either x is a matrix and y is a vector, r must be a matrix of size (length(lags), size(x, 2)); if x is a vector and y is a matrix, r must be a matrix of size (length(lags), size(y, 2)). If both x and y are matrices, r must be a three-dimensional array of size (length(lags), size(x, 2), size(y, 2)).\n\nThe output is not normalized. See crosscor! for a function with normalization.\n\n\n\n"
},

{
    "location": "signalcorr.html#StatsBase.crosscor",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.crosscor",
    "category": "Function",
    "text": "crosscor(x, y, [lags]; demean=true)\n\nCompute the cross correlation between real-valued vectors or matrices x and y, optionally specifying the lags. demean specifies whether the respective means of x and y should be subtracted from them before computing their cross correlation.\n\nIf both x and y are vectors, return a vector of the same length as lags. Otherwise, compute cross covariances between each pairs of columns in x and y.\n\nWhen left unspecified, the lags used are the integers from -min(size(x,1)-1, 10*log10(size(x,1))) to min(size(x,1), 10*log10(size(x,1))).\n\nThe output is normalized by sqrt(var(x)*var(y)). See crosscov for the unnormalized form.\n\n\n\n"
},

{
    "location": "signalcorr.html#StatsBase.crosscor!",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.crosscor!",
    "category": "Function",
    "text": "crosscor!(r, x, y, lags; demean=true)\n\nCompute the cross correlation between real-valued vectors or matrices x and y at lags and store the result in r. demean specifies whether the respective means of x and y should be subtracted from them before computing their cross correlation.\n\nIf both x and y are vectors, r must be a vector of the same length as lags. If either x is a matrix and y is a vector, r must be a matrix of size (length(lags), size(x, 2)); if x is a vector and y is a matrix, r must be a matrix of size (length(lags), size(y, 2)). If both x and y are matrices, r must be a three-dimensional array of size (length(lags), size(x, 2), size(y, 2)).\n\nThe output is normalized by sqrt(var(x)*var(y)). See crosscov! for the unnormalized form.\n\n\n\n"
},

{
    "location": "signalcorr.html#Cross-covariance-and-Cross-correlation-1",
    "page": "Correlation Analysis of Signals",
    "title": "Cross-covariance and Cross-correlation",
    "category": "section",
    "text": "crosscov\ncrosscov!\ncrosscor\ncrosscor!"
},

{
    "location": "signalcorr.html#StatsBase.pacf",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.pacf",
    "category": "Function",
    "text": "pacf(X, lags; method=:regression)\n\nCompute the partial autocorrelation function (PACF) of a real-valued vector or matrix X at lags. method designates the estimation method. Recognized values are :regression, which computes the partial autocorrelations via successive regression models, and :yulewalker, which computes the partial autocorrelations using the Yule-Walker equations.\n\nIf x is a vector, return a vector of the same length as lags. If x is a matrix, return a matrix of size (length(lags), size(x, 2)), where each column in the result corresponds to a column in x.\n\n\n\n"
},

{
    "location": "signalcorr.html#StatsBase.pacf!",
    "page": "Correlation Analysis of Signals",
    "title": "StatsBase.pacf!",
    "category": "Function",
    "text": "pacf!(r, X, lags; method=:regression)\n\nCompute the partial autocorrelation function (PACF) of a matrix X at lags and store the result in r. method designates the estimation method. Recognized values are :regression, which computes the partial autocorrelations via successive regression models, and :yulewalker, which computes the partial autocorrelations using the Yule-Walker equations.\n\nr must be a matrix of size (length(lags), size(x, 2)).\n\n\n\n"
},

{
    "location": "signalcorr.html#Partial-Autocorrelation-Function-1",
    "page": "Correlation Analysis of Signals",
    "title": "Partial Autocorrelation Function",
    "category": "section",
    "text": "pacf\npacf!"
},

{
    "location": "misc.html#",
    "page": "Miscellaneous Functions",
    "title": "Miscellaneous Functions",
    "category": "page",
    "text": ""
},

{
    "location": "misc.html#StatsBase.rle",
    "page": "Miscellaneous Functions",
    "title": "StatsBase.rle",
    "category": "Function",
    "text": "rle(v) -> (vals, lens)\n\nReturn the run-length encoding of a vector as a tuple. The first element of the tuple is a vector of values of the input and the second is the number of consecutive occurrences of each element.\n\nExamples\n\njulia> using StatsBase\n\njulia> rle([1,1,1,2,2,3,3,3,3,2,2,2])\n([1, 2, 3, 2], [3, 2, 4, 3])\n\n\n\n"
},

{
    "location": "misc.html#StatsBase.inverse_rle",
    "page": "Miscellaneous Functions",
    "title": "StatsBase.inverse_rle",
    "category": "Function",
    "text": "inverse_rle(vals, lens)\n\nReconstruct a vector from its run-length encoding (see rle). vals is a vector of the values and lens is a vector of the corresponding run lengths.\n\n\n\n"
},

{
    "location": "misc.html#StatsBase.levelsmap",
    "page": "Miscellaneous Functions",
    "title": "StatsBase.levelsmap",
    "category": "Function",
    "text": "levelsmap(a)\n\nConstruct a dictionary that maps each of the n unique values in a to a number between 1 and n.\n\n\n\n"
},

{
    "location": "misc.html#StatsBase.indexmap",
    "page": "Miscellaneous Functions",
    "title": "StatsBase.indexmap",
    "category": "Function",
    "text": "indexmap(a)\n\nConstruct a dictionary that maps each unique value in a to the index of its first occurrence in a.\n\n\n\n"
},

{
    "location": "misc.html#StatsBase.indicatormat",
    "page": "Miscellaneous Functions",
    "title": "StatsBase.indicatormat",
    "category": "Function",
    "text": "indicatormat(x, k::Integer; sparse=false)\n\nConstruct a boolean matrix I of size (k, length(x)) such that I[x[i], i] = true and all other elements are set to false. If sparse is true, the output will be a sparse matrix, otherwise it will be dense (default).\n\nExamples\n\njulia> using StatsBase\n\njulia> indicatormat([1 2 2], 2)\n2×3 Array{Bool,2}:\n  true  false  false\n false   true   true\n\n\n\nindicatormat(x, c=sort(unique(x)); sparse=false)\n\nConstruct a boolean matrix I of size (length(c), length(x)). Let ci be the index of x[i] in c. Then I[ci, i] = true and all other elements are false.\n\n\n\n"
},

{
    "location": "misc.html#Base.midpoints",
    "page": "Miscellaneous Functions",
    "title": "Base.midpoints",
    "category": "Function",
    "text": "midpoints(v)\n\nCalculate the midpoints (pairwise mean of consecutive elements).\n\n\n\n"
},

{
    "location": "misc.html#Miscellaneous-Functions-1",
    "page": "Miscellaneous Functions",
    "title": "Miscellaneous Functions",
    "category": "section",
    "text": "rle\ninverse_rle\nlevelsmap\nindexmap\nindicatormat\nmidpoints"
},

{
    "location": "statmodels.html#",
    "page": "Abstraction for Statistical Models",
    "title": "Abstraction for Statistical Models",
    "category": "page",
    "text": ""
},

{
    "location": "statmodels.html#StatsBase.adjr2",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.adjr2",
    "category": "Function",
    "text": "adjr2(obj::StatisticalModel, variant::Symbol)\nadjr²(obj::StatisticalModel, variant::Symbol)\n\nAdjusted coefficient of determination (adjusted R-squared).\n\nFor linear models, the adjusted R² is defined as 1 - (1 - (1-R^2)(n-1)(n-p)), with R^2 the coefficient of determination, n the number of observations, and p the number of coefficients (including the intercept). This definition is generally known as the Wherry Formula I.\n\nFor other models, one of the several pseudo R² definitions must be chosen via variant. The only currently supported variant is :MacFadden, defined as 1 - (log L - k)log L0. In this formula, L is the likelihood of the model, L0 that of the null model (the model including only the intercept). These two quantities are taken to be minus half deviance of the corresponding models. k is the number of consumed degrees of freedom of the model (as returned by dof).\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.aic",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.aic",
    "category": "Function",
    "text": "aic(obj::StatisticalModel)\n\nAkaike's Information Criterion, defined as -2 log L + 2k, with L the likelihood of the model, and k its number of consumed degrees of freedom (as returned by dof).\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.aicc",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.aicc",
    "category": "Function",
    "text": "aicc(obj::StatisticalModel)\n\nCorrected Akaike's Information Criterion for small sample sizes (Hurvich and Tsai 1989), defined as -2 log L + 2k + 2k(k-1)(n-k-1), with L the likelihood of the model, k its number of consumed degrees of freedom (as returned by dof), and n the number of observations (as returned by nobs).\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.bic",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.bic",
    "category": "Function",
    "text": "bic(obj::StatisticalModel)\n\nBayesian Information Criterion, defined as -2 log L + k log n, with L the likelihood of the model,  k its number of consumed degrees of freedom (as returned by dof), and n the number of observations (as returned by nobs).\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.coef",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.coef",
    "category": "Function",
    "text": "coef(obj::StatisticalModel)\n\nReturn the coefficients of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.coefnames",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.coefnames",
    "category": "Function",
    "text": "coefnames(obj::StatisticalModel)\n\nReturn the names of the coefficients.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.coeftable",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.coeftable",
    "category": "Function",
    "text": "coeftable(obj::StatisticalModel)\n\nReturn a table of class CoefTable with coefficients and related statistics.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.confint",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.confint",
    "category": "Function",
    "text": "confint(obj::StatisticalModel)\n\nCompute confidence intervals for coefficients.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.deviance",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.deviance",
    "category": "Function",
    "text": "deviance(obj::StatisticalModel)\n\nReturn the deviance of the model relative to a reference, which is usually when applicable the saturated model. It is equal, up to a constant, to -2 log L, with L the likelihood of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.dof",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.dof",
    "category": "Function",
    "text": "dof(obj::StatisticalModel)\n\nReturn the number of degrees of freedom consumed in the model, including when applicable the intercept and the distribution's dispersion parameter.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.fit",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.fit",
    "category": "Function",
    "text": "fit(Histogram, data[, weight][, edges]; closed=:right, nbins)\n\nFit a histogram to data.\n\nArguments\n\ndata: either a vector (for a 1-dimensional histogram), or a tuple of vectors of equal length (for an n-dimensional histogram).\nweight: an optional AbstractWeights (of the same length as the data vectors), denoting the weight each observation contributes to the bin. If no weight vector is supplied, each observation has weight 1.\nedges: a vector (typically an AbstractRange object), or tuple of vectors, that gives the edges of the bins along each dimension. If no edges are provided, these are determined from the data.\n\nKeyword arguments\n\nclosed=:right: if :left, the bin intervals are left-closed [a,b); if :right (the default), intervals are right-closed (a,b].\nnbins: if no edges argument is supplied, the approximate number of bins to use along each dimension (can be either a single integer, or a tuple of integers).\n\nExamples\n\n# Univariate\nh = fit(Histogram, rand(100))\nh = fit(Histogram, rand(100), 0:0.1:1.0)\nh = fit(Histogram, rand(100), nbins=10)\nh = fit(Histogram, rand(100), weights(rand(100)), 0:0.1:1.0)\nh = fit(Histogram, [20], 0:20:100)\nh = fit(Histogram, [20], 0:20:100, closed=:left)\n\n# Multivariate\nh = fit(Histogram, (rand(100),rand(100)))\nh = fit(Histogram, (rand(100),rand(100)),nbins=10)\n\n\n\nFit a statistical model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.fit!",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.fit!",
    "category": "Function",
    "text": "Fit a statistical model in-place.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.loglikelihood",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.loglikelihood",
    "category": "Function",
    "text": "loglikelihood(obj::StatisticalModel)\n\nReturn the log-likelihood of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.nobs",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.nobs",
    "category": "Function",
    "text": "nobs(obj::StatisticalModel)\n\nReturn the number of independent observations on which the model was fitted. Be careful when using this information, as the definition of an independent observation may vary depending on the model, on the format used to pass the data, on the sampling plan (if specified), etc.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.nulldeviance",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.nulldeviance",
    "category": "Function",
    "text": "nulldeviance(obj::StatisticalModel)\n\nReturn the deviance of the null model, that is the one including only the intercept.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.r2",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.r2",
    "category": "Function",
    "text": "r2(obj::StatisticalModel, variant::Symbol)\nr²(obj::StatisticalModel, variant::Symbol)\n\nCoefficient of determination (R-squared).\n\nFor a linear model, the R² is defined as ESSTSS, with ESS the explained sum of squares and TSS the total sum of squares, and variant can be omitted.\n\nFor other models, one of several pseudo R² definitions must be chosen via variant. Supported variants are:\n\n:MacFadden (a.k.a. likelihood ratio index), defined as 1 - log Llog L0.\n:CoxSnell, defined as 1 - (L0L)^2n\n:Nagelkerke, defined as (1 - (L0L)^2n)(1 - L0^2n), with n the number\n\nof observations (as returned by nobs).\n\nIn the above formulas, L is the likelihood of the model, L0 that of the null model (the model including only the intercept). These two quantities are taken to be minus half deviance of the corresponding models.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.stderr",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.stderr",
    "category": "Function",
    "text": "stderr(obj::StatisticalModel)\n\nReturn the standard errors for the coefficients of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.vcov",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.vcov",
    "category": "Function",
    "text": "vcov(obj::StatisticalModel)\n\nReturn the variance-covariance matrix for the coefficients of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.dof_residual",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.dof_residual",
    "category": "Function",
    "text": "dof_residual(obj::RegressionModel)\n\nReturn the residual degrees of freedom of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.fitted",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.fitted",
    "category": "Function",
    "text": "fitted(obj::RegressionModel)\n\nReturn the fitted values of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.modelmatrix",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.modelmatrix",
    "category": "Function",
    "text": "modelmatrix(obj::RegressionModel)\n\nReturn the model matrix (a.k.a. the design matrix).\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.model_response",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.model_response",
    "category": "Function",
    "text": "model_response(obj::RegressionModel)\n\nReturn the model response (a.k.a. the dependent variable).\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.predict",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.predict",
    "category": "Function",
    "text": "predict(obj::RegressionModel, [newX])\n\nForm the predicted response of model obj. An object with new covariate values newX can be supplied, which should have the same type and structure as that used to fit obj; e.g. for a GLM it would generally be a DataFrame with the same variable names as the original predictors.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.predict!",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.predict!",
    "category": "Function",
    "text": "predict!\n\nIn-place version of predict.\n\n\n\n"
},

{
    "location": "statmodels.html#StatsBase.residuals",
    "page": "Abstraction for Statistical Models",
    "title": "StatsBase.residuals",
    "category": "Function",
    "text": "residuals(obj::RegressionModel)\n\nReturn the residuals of the model.\n\n\n\n"
},

{
    "location": "statmodels.html#Abstraction-for-Statistical-Models-1",
    "page": "Abstraction for Statistical Models",
    "title": "Abstraction for Statistical Models",
    "category": "section",
    "text": "This package defines an abstract type StatisticalModel, and an abstract subtype RegressionModel.Particularly, instances of StatisticalModel implement the following methods.adjr2\naic\naicc\nbic\ncoef\ncoefnames\ncoeftable\nconfint\ndeviance\ndof\nfit\nfit!\nloglikelihood\nnobs\nnulldeviance\nr2\nstderr\nvcovRegressionModel extends StatisticalModel by implementing the following additional methods.dof_residual\nfitted\nmodelmatrix\nmodel_response\npredict\npredict!\nresiduals"
},

]}
