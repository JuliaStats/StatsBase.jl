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

    # others
    ## Types
    StatisticalModel,
    RegressionModel,
    
    ## Fcuntions
    coef,
    coeftable,
    confint,
    ecdf,
    findat,
    indicators,
    inverse_rle,
    loglikelihood,
    nobs,
    predict,
    residuals,
    model_response,
    rle,
    stderr,
    vcov

    include("means.jl")
    include("scalar_stats.jl")
    include("intstats.jl")
    include("corr.jl")
    include("others.jl")

end # module
