## StatsBase.jl Release Notes

### Changes from v0.5 to v0.6

* Support of weighted sampling without replacements

* Implement several new sampling algorithms and a new polyalgorithm for sampling (for both with and without replacements). The new polyalgorithms choose the proper routine depending on the size of the pool and the number of samples to be drawn. (See [documentation](http://statsbasejl.readthedocs.org/en/latest/sampling.html)). The choices are determined based on benchmarks.

* New function ``scattermat``: scatter matrix and weighted scatter matrix.

* Weighted covariance matrix.

* Weighted variance and standard deviation (both over an entire array or along a given dimension).

* The function ``middle`` is moved to Julia Base.

* Functions to compute certain combinations of statistics (more efficiently): ``mean_and_var``, ``mean_and_std``, and ``mean_and_cov``.

* New function ``zscore``: computation of z-scores.


### Changes from v0.6 to v0.7

* Remove functions deprecated before v0.5: ``gmean``, ``hmean``, ``kde``, ``acf``, ``ccf`` and ``indicators``.

* Add functions to compute deviations (moved from [MLBase](https://github.com/JuliaStats/MLBase.jl)): 
  ``countne``, ``counteq``, ``sqL2dist``, ``L2dist``, ``L1dist``, ``Linfdist``, ``gkldiv``, ``meanad``, ``maxad``, ``msd``, ``rmsd``, and ``psnr`` (See [documentation here](http://statsbasejl.readthedocs.org/en/latest/deviation.html)).

