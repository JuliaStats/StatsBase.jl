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
    autocor,
    cor_spearman,
    cor_kendall,    

    # others
    rle,
    inverse_rle,
    ecdf,
    findat
     

    include("scalar_stats.jl")
    include("weighted_stats.jl")
    include("corr.jl")
    include("others.jl")

end # module
