module Stats
    import Base.quantile
    import Base.LinAlg: BlasReal

    export

    # means
    gmean,
    hmean,
    wmean,

    # scalar_stats 
    variation, 
    sem,
    mad,
    skewness,
    kurtosis,
    minmax,
    midrange,
    range,
    tiedrank,
    percentile,
    quantile,
    quartile,
    quintile,
    decile,
    iqr,
    table,
    mode,
    modes,
    summarystats,
    describe,

    # intstats
    addcounts!,
    addwcounts!,
    counts,
    wcounts,

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

    include("means.jl")
    include("scalarstats.jl")
    include("intstats.jl")
    include("ranking.jl")
    include("toeplitzsolvers.jl")
    include("corr.jl")
    include("statmodels.jl")
    include("others.jl")

end # module
