module StatsBase
    using ArrayViews

    import Base: length, isempty, eltype, values, sum, mean, mean!, show, quantile
    import Base: rand, rand!
    import Base.LinAlg: BlasReal, BlasFloat
    import Base.Cartesian: @ngenerate, @nloops, @nref, @nextract

    export

    # common
    WeightVec, weights,

    # mathfuns
    xlogx,       # x * log(x)
    xlogy,       # x * log(y)
    logistic,    # 1 / (1 + exp(-x))
    logit,       # log(x / (1 - x))
    softplus,    # log(1 + exp(x))
    invsoftplus, # log(exp(x) - 1)
    logsumexp,   # log(exp(x) + exp(y)) or log(sum(exp(x)))
    softmax,
    softmax!,

    # means
    geomean,     # geometric mean
    harmmean,    # harmonic mean
    trimmean,    # trimmed mean
    wmean,       # weighted mean
    wsum,        # weighted sum with vector as second argument
    wsum!,       # in-place weighted sum across dimensions

    # scalar_stats 
    skewness,   # (standardized) skewness
    kurtosis,   # (excessive) kurtosis
    variation,  # ratio of standard deviation to mean
    sem,        # standard error of the mean, i.e. sqrt(var / n)
    mad,        # median absolute deviation
    middle,     # the mean of two real numbers
    entropy,    # the entropy of a probability vector
    crossentropy,   # cross entropy between two probability vectors
    kldivergence,   # K-L divergence between two probability vectors
    percentile, # quantile using percentage (instead of fraction) as argument
    iqr,        # interquatile range 
    nquantile,  # quantiles at [0:n]/n
    mode,       # find a mode from data 
    modes,      # find all modes from data
    summarystats,   # summary statistics
    describe,       # print the summary statistics

    # counts
    addcounts!,     # add counts to an accumulating array or map
    counts,         # count integer values in given arrays
    proportions,    # proportions of integer values in given arrays 
                    # (normalized version of counts)                    
    countmap,       # count distinct values and return a map
    proportionmap,  # proportions of distinct values returned as a map                     

    # ranking
    ordinalrank,    # ordinal ranking ("1234" ranking)
    competerank,    # competition ranking ("1 2 2 4" ranking)
    denserank,      # dense ranking ("1 2 2 3" ranking)
    tiedrank,       # tied ranking ("1 2.5 2.5 4" ranking)

    # rankcorr
    corspearman,       # spearman's rank correlation
    corkendall,        # kendall's rank correlation

    # signalcorr
    autocov!, autocov,      # auto covariance
    autocor!, autocor,      # auto correlation
    crosscov!, crosscov,    # cross covariance
    crosscor!, crosscor,    # cross correlation
    pacf!, pacf,  # partial auto-correlation

    # sampling
    samplepair,     # draw a pair of distinct elements   
    sample,         # sampling from a population 
    sample!,        # sampling from a population, with pre-allocated output
    wsample,        # sampling from a population with weights
    wsample!,       # weighted sampling, with pre-allocated output

    # empirical
    ecdf,           # empirical cumulative distribution function

    Histogram,
    hist,
    histrange,
    midpoints,

    # misc
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

    coef,
    coeftable,
    confint,
    deviance,
    fit,
    fitted,
    loglikelihood,
    nobs,
    stderr,
    vcov,
    predict,
    predict!,
    residuals,
    model_response

    # source files

    include("common.jl")
    include("mathfuns.jl")
    include("weights.jl")
    
    include("misc.jl")
    include("means.jl")
    include("scalarstats.jl")
    include("counts.jl")
    include("ranking.jl")
    include("toeplitzsolvers.jl")
    include("rankcorr.jl")
    include("signalcorr.jl")
    include("rand.jl")
    include("sampling.jl")
    include("empirical.jl")
    include("hist.jl")

    include("statmodels.jl")

    include("deprecates.jl")
    

end # module
