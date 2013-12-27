module Stats
    import Base: length, isempty, eltype, values, sum, mean, show, quantile
    import Base.LinAlg: BlasReal

    export

    # common
    WeightVec, weights,

    # means
    geomean, gmean,     # geometric mean
    harmmean, hmean,    # harmonic mean
    trimmean,           # trimmed mean
    wmean,              # weighted mean

    # scalar_stats 
    skewness,   # (standardized) skewness
    kurtosis,   # (excessive) kurtosis
    variation,  # ratio of standard deviation to mean
    sem,        # standard error of the mean, i.e. sqrt(var / n)
    mad,        # median absolute deviation
    minmax,     # obtain min & max in a single pass
    middle,     # the mean of two real numbers
    midrange,   # the mean of minimum and maximum
    range,      # the difference between maximum and minimum
    prctile,    # quantile using percentage (instead of fraction) as argument
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

    # corr
    acf,
    ccf,
    cor_spearman,
    cor_kendall,
    pacf,    

    # others
    ## Types
    StatisticalModel,
    RegressionModel,
    
    ## Functions
    coef,
    coeftable,
    confint,
    deviance,
    durbin,
    ecdf,
    findat,
    indicators,
    inverse_rle,
    levinson,
    loglikelihood,
    nobs,
    predict,
    residuals,
    model_response,
    rle,
    stderr,
    vcov

    include("common.jl")
    include("means.jl")
    include("scalarstats.jl")
    include("counts.jl")
    include("ranking.jl")
    include("toeplitzsolvers.jl")
    include("corr.jl")
    include("statmodels.jl")
    include("others.jl")

end # module
