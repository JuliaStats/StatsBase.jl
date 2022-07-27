module StatsBase

import Base: length, size, isempty, values, sum, show, maximum, minimum, extrema
import Base.Cartesian: @nloops, @nref, @nextract
using Base: @irrational, @propagate_inbounds
using DataAPI
import DataAPI: describe
import DataStructures: heapify!, heappop!, percolate_down!
using SortingAlgorithms
using Missings
using LogExpFunctions: xlogx, xlogy

using Statistics
using LinearAlgebra
using Random
using Printf
using SparseArrays
import Random: rand, rand!
import LinearAlgebra: BlasReal, BlasFloat
import Statistics: mean, mean!, var, varm, varm!, std, stdm, cov, covm,
                   cor, corm, cov2cor!, unscaled_covzm, quantile, sqrt!,
                   median, middle
using StatsAPI: StatisticalModel, RegressionModel
import StatsAPI: pairwise, pairwise!, params, params!,
                 fitted, response, responsename, meanresponse, modelmatrix,
                 crossmodelmatrix, leverage, cooksdistance, residuals, predict,
                 predict!, dof_residual, coef, coefnames, coeftable, confint,
                 deviance, islinear, nulldeviance, loglikelihood, nullloglikelihood,
                 loglikelihood, loglikelihood, score, nobs, dof, mss, rss,
                 informationmatrix, stderror, vcov, weights, isfitted, fit, fit!,
                 aic, aicc, bic, r2, r², adjr2, adjr²

    ## tackle compatibility issues

export

    ## functions defined in Statistics
    cor,
    cov,
    mean,
    mean!,
    median,
    median!,
    quantile,
    quantile!,
    std,
    var,

    ## weights
    AbstractWeights,    # abstract type to represent any weight vector
    Weights,            # to represent a generic weight vector
    AnalyticWeights,    # to represent an analytic/precision/reliability weight vector
    FrequencyWeights,   # to representing a frequency/case/repeat weight vector
    ProbabilityWeights, # to representing a probability/sampling weight vector
    UnitWeights,        # to representing a uniform weight vector
    weights,            # construct a generic Weights vector
    aweights,           # construct an AnalyticWeights vector
    fweights,           # construct a FrequencyWeights vector
    pweights,           # construct a ProbabilityWeights vector
    eweights,           # construct an exponential Weights vector
    uweights,           # construct an UnitWeights vector
    wsum,               # weighted sum with vector as second argument
    wsum!,              # weighted sum across dimensions with provided storage

    ## moments
    skewness,       # (standardized) skewness
    kurtosis,       # (excessive) kurtosis
    moment,         # central moment of given order
    cumulant,       # cumulant of given order
    mean_and_var,   # (mean, var)
    mean_and_std,   # (mean, std)
    mean_and_cov,   # (mean, cov)

    ## scalarstats
    geomean,     # geometric mean
    harmmean,    # harmonic mean
    genmean,     # generalized/power mean
    middle,      # the mean of two real numbers
    mode,        # find a mode from data (the first one)
    modes,       # find all modes from data

    zscore,      # compute Z-scores
    zscore!,     # compute Z-scores inplace or to a pre-allocated array

    percentile,     # quantile using percentage (instead of fraction) as argument
    nquantile,      # quantiles at [0:n]/n
    quantilerank,   # quantile-position (0-1) of a value relative to a collection
    percentilerank, # percentile-position (0-100) of a value relative to a collection

    span,        # The range minimum(x):maximum(x)
    variation,   # ratio of standard deviation to mean
    sem,         # standard error of the mean, i.e. sqrt(var / n)
    mad,         # median absolute deviation
    iqr,         # interquatile range

    genvar,      # generalized variance
    totalvar,    # total variation

    entropy,        # the entropy of a probability vector
    renyientropy,   # the Rényi (generalised) entropy of a probability vector
    crossentropy,   # cross entropy between two probability vectors
    kldivergence,   # K-L divergence between two probability vectors

    summarystats,   # summary statistics
    describe,       # print the summary statistics

    # deviation
    counteq,        # count the number of equal pairs
    countne,        # count the number of non-equal pairs
    sqL2dist,       # squared L2 distance between two arrays
    L2dist,         # L2 distance between two arrays
    L1dist,         # L1 distance between two arrays
    Linfdist,       # L-inf distance between two arrays
    gkldiv,         # (Generalized) Kullback-Leibler divergence between two vectors
    meanad,         # mean absolute deviation
    maxad,          # maximum absolute deviation
    msd,            # mean squared deviation
    rmsd,           # root mean squared deviation
    psnr,           # peak signal-to-noise ratio (in dB)

    # cov
    scattermat,     # scatter matrix (i.e. unnormalized covariance)
    cov2cor,        # converts a covariance matrix to a correlation matrix
    cor2cov,        # converts a correlation matrix to a covariance matrix
    CovarianceEstimator, # abstract type for covariance estimators
    SimpleCovariance,    # simple covariance estimator

    ## counts
    addcounts!,     # add counts to an accumulating array or map
    counts,         # count integer values in given arrays
    proportions,    # proportions of integer values in given arrays
                    # (normalized version of counts)
    countmap,       # count distinct values and return a map
    proportionmap,  # proportions of distinct values returned as a map

    ## ranking
    ordinalrank,    # ordinal ranking ("1234" ranking)
    competerank,    # competition ranking ("1 2 2 4" ranking)
    denserank,      # dense ranking ("1 2 2 3" ranking)
    tiedrank,       # tied ranking ("1 2.5 2.5 4" ranking)

    ## rankcorr
    corspearman,       # spearman's rank correlation
    corkendall,        # kendall's rank correlation

    ## partialcor
    partialcor,        # partial correlation

    ## signalcorr
    autocov!, autocov,      # auto covariance
    autocor!, autocor,      # auto correlation
    crosscov!, crosscov,    # cross covariance
    crosscor!, crosscor,    # cross correlation
    pacf!, pacf,            # partial auto-correlation

    ## sampling
    samplepair,     # draw a pair of distinct elements   
    sample,         # sampling from a population
    sample!,        # sampling from a population, with pre-allocated output
    wsample,        # sampling from a population with weights
    wsample!,       # weighted sampling, with pre-allocated output

    ## empirical
    ecdf,           # empirical cumulative distribution function
    ECDF,           # type for empirical cumulative distribution function

    AbstractHistogram,
    Histogram,
    midpoints,
    # histrange,

    ## robust
    trim,           # trimmed set
    trim!,          # trimmed set
    winsor,         # Winsorized set
    winsor!,        # Winsorized set
    trimvar,        # variance of the mean of a trimmed set

    ## misc
    rle,            # run-length encoding
    inverse_rle,    # inverse run-length encoding
    indexmap,       # construct a map from element to index
    levelsmap,      # construct a map from n unique elements to [1, ..., n]
    indicatormat,   # construct indicator matrix
    pairwise,       # pairwise application of functions
    pairwise!,      # pairwise! application of functions

    # statistical models
    CoefTable,
    StatisticalModel,
    RegressionModel,

    adjr2,
    adjr²,
    aic,
    aicc,
    bic,
    coef,
    coefnames,
    coeftable,
    confint,
    cooksdistance,
    crossmodelmatrix,
    deviance,
    dof,
    dof_residual,
    fit,
    fit!,
    fitted,
    informationmatrix,
    isfitted,
    islinear,
    leverage,
    loglikelihood,
    meanresponse,
    modelmatrix,
    mss,
    response,
    responsename,
    nobs,
    nulldeviance,
    nullloglikelihood,
    rss,
    score,
    stderror,
    vcov,
    predict,
    predict!,
    residuals,
    r2,
    r²,

    ConvergenceException,

    # data standardization
    standardize,
    AbstractDataTransform, # the type to represent a abstract data transformation
    ZScoreTransform,       # the type to represent a z-score data transformation
    UnitRangeTransform,    # the type to represent a 0-1 data transformation

    # reliability
    CronbachAlpha,        # the type to represent Cronbach's alpha scores
    cronbachalpha         # function to compute Cronbach's alpha scores

# source files

include("common.jl")
include("weights.jl")
include("moments.jl")
include("scalarstats.jl")
include("robust.jl")
include("deviation.jl")
include("cov.jl")
include("counts.jl")
include("ranking.jl")
include("toeplitzsolvers.jl")
include("rankcorr.jl")
include("signalcorr.jl")
include("partialcor.jl")
include("empirical.jl")
include("hist.jl")
include("pairwise.jl")
include("reliability.jl")
include("misc.jl")

include("sampling.jl")
include("statmodels.jl")

include("transformations.jl")

include("deprecates.jl")

end # module
