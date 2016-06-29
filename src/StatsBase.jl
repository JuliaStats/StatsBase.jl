__precompile__(true)

module StatsBase
    using Compat
    import Compat.String
    import Compat.view
    using ArrayViews
    using StatsFuns

    import Base: length, isempty, eltype, values, sum, mean, mean!, show, quantile
    import Base: rand, rand!
    import Base.LinAlg: BlasReal, BlasFloat
    import Base.Cartesian: @nloops, @nref, @nextract

    ## tackle compatibility issues

    ## import mathfuns, which were migrated to StatsFuns
    import StatsFuns: xlogx, xlogy, logistic, logit,
                      softplus, invsoftplus,
                      logsumexp, softmax, softmax!

    if VERSION < v"0.5.0-dev"
        # reexport from ArrayViews
        export view
    end

    export

    ## mathfuns (TODO: removed after a certain period)
    xlogx,       # x * log(x)
    xlogy,       # x * log(y)
    logistic,    # 1 / (1 + exp(-x))
    logit,       # log(x / (1 - x))
    softplus,    # log(1 + exp(x))
    invsoftplus, # log(exp(x) - 1)
    logsumexp,   # log(exp(x) + exp(y)) or log(sum(exp(x)))
    softmax,
    softmax!,

    ## weights
    WeightVec,   # the type to represent a weight vector
    weights,     # construct a weight vector
    wsum,        # weighted sum with vector as second argument
    wsum!,       # weighted sum across dimensions with provided storage
    wmean,       # weighted mean
    wmean!,      # weighted mean across dimensions with provided storage
    wmedian,     # weighted median
    wquantile,   # weighted quantile

    ## moments
    skewness,       # (standardized) skewness
    kurtosis,       # (excessive) kurtosis
    moment,         # central moment of given order
    mean_and_var,   # (mean, var)
    mean_and_std,   # (mean, std)
    mean_and_cov,   # (mean, cov)

    ## scalarstats
    geomean,     # geometric mean
    harmmean,    # harmonic mean
    trimmean,    # trimmed mean
    middle,      # the mean of two real numbers
    mode,        # find a mode from data (the first one)
    modes,       # find all modes from data

    zscore,      # compute Z-scores
    zscore!,     # compute Z-scores inplace or to a pre-allocated array

    percentile,  # quantile using percentage (instead of fraction) as argument
    nquantile,   # quantiles at [0:n]/n

    span,        # The range minimum(x):maximum(x)
    variation,   # ratio of standard deviation to mean
    sem,         # standard error of the mean, i.e. sqrt(var / n)
    mad,         # median absolute deviation
    iqr,         # interquatile range

    entropy,        # the entropy of a probability vector
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

    AbstractHistogram,
    Histogram,
    hist,
    # histrange,
    midpoints,

    ## misc
    rle,            # run-length encoding
    inverse_rle,    # inverse run-length encoding
    indexmap,       # construct a map from element to index
    levelsmap,      # construct a map from n unique elements to [1, ..., n]
    findat,         # find the position within a for elements in b
    indicatormat,   # construct indicator matrix

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
    coeftable,
    confint,
    deviance,
    df,
    df_residual,
    fit,
    fit!,
    fitted,
    loglikelihood,
    nobs,
    nulldeviance,
    nullloglikelihood,
    stderr,
    vcov,
    predict,
    predict!,
    residuals,
    r2,
    r²,
    model_response

    # source files

    include("common.jl")
    include("weights.jl")
    include("moments.jl")
    include("scalarstats.jl")
    include("deviation.jl")
    include("cov.jl")
    include("counts.jl")
    include("ranking.jl")
    include("toeplitzsolvers.jl")
    include("rankcorr.jl")
    include("signalcorr.jl")
    include("rand.jl")
    include("empirical.jl")
    include("hist.jl")
    include("misc.jl")

    include("sampling.jl")
    include("statmodels.jl")

    include("deprecates.jl")

end # module
