Correlation Analysis 
======================

The package provides functions to perform correlation analysis of sequential signals. 

Autocovariance and Autocorrelation
------------------------------------

- **autocov** (x, lags[; demean=true])

    Compute auto-covariance of ``x`` at specified lags. If ``x`` is a vector, it returns a vector of th same length of ``lags``. If ``x`` is a matrix, it returns a matrix of size ``(length(lags), size(x,2))``, where each column in the result corresponding to a column in ``x``. 

    Here, ``demean`` is a keyword argument (default value is ``true``), which means that the function will subtract each ``x`` from its mean before computing the results. Otherwise, ``x`` is considered as having been centered. 

- **autocov** (x[; demean=true])

    Compute auto-covariance at default lags.  

- **autocov!** (r, x, lags[; demean=true])

    Write the computed auto-covariance to ``r``.   
 
- **autocor** (x, lags[; demean=true])

    Compute auto-correlation of ``x`` at specified lags. 

- **autocor** (x[; demean=true])

    Compute auto-correlation at default lags.  

- **autocor!** (r, x, lags[; demean=true])

    Write the computed auto-correlation to ``r``.   


Cross-covariance and Cross-correlation
---------------------------------------

- **crosscov** (x, y, lags[; demean=true])

    Compute cross-covariance of ``x`` at specified lags. If both ``x`` and ``y`` are vectors, it returns a vector of th same length of ``lags``. Otherwise, it computes cross covariances between each pairs of columns in ``x`` and ``y``.

- **crosscov** (x, y[; demean=true])

    Compute cross-covariance at default lags.  

- **crosscov!** (r, x, y, lags[; demean=true])

    Write the computed cross-covariance to ``r``.   
 
- **crosscor** (x, y, lags[; demean=true])

    Compute cross-correlation of ``x`` at specified lags. 

- **crosscor** (x, y[; demean=true])

    Compute cross-correlation at default lags.  

- **crosscor!** (r, x, y, lags[; demean=true])

    Write the computed cross-correlation to ``r``.   


Partial Autocorrelation Function
----------------------------------

- **pacf** (x, lags[; method=:regression])

    Compute partial auto-correlation of ``x`` at specified lags. If ``x`` is a vector, it returns a vector of th same length of ``lags``. If ``x`` is a matrix, it returns a matrix of size ``(length(lags), size(x,2))``, where each column in the result corresponding to a column in ``x``. 

    Here, ``method`` is a keyword argument to specify the choice of algorithm, which can be either ``:regresion`` or ``:yulewalker``. The default value is ``:regression``. 


- **pacf!** (r, x, lags[; method=:regression])

    Write the computed partial auto-correlation to ``r``.


- **corspearman** (x, y)

    Spearman's rank correlation. Here, ``x`` and ``y`` can be either real vectors or matrices. When ``xx`` and ``y`` are matrices, it computes the correlations between their columns (pairwisely).

- **corkendall** (x, y)   

    Kendall's rank correlation. Here, ``x`` and ``y`` can be either real vectors or matrices. When ``xx`` and ``y`` are matrices, it computes the correlations between their columns (pairwisely).

