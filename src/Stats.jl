module Stats
    import Base.quantile

    export

    # scalar_stats
    gmean, 
    hmean, 
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
    describe,

    # weighted stats
    wmean,

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

    include("scalar_stats.jl")
    include("weighted_stats.jl")
    include("corr.jl")
    include("others.jl")

end # module
