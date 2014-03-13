Empirical Estimation
=====================

Empirical Cumulative Distribution Function
-------------------------------------------

.. py:function:: ecdf(x)

  Return an empirical cumulative distribution function based on a vector of samples given in ``x``. 

  **Note:** this is a higher-level function that returns a function, which can then be applied to evaluate CDF values on other samples.


Kernel Density Estimation
---------------------------

.. py:function:: kde(data[; width=NaN, npoints=2048]) 

  Kernel density estimation. This function returns an instance of ``UnivariateKDE`` defined as below:
  
  .. code-block:: julia

    immutable UnivariateKDE
        x::Vector{Float64}
        density::Vector{Float64}
    end

  **Keyword arguments***

  - ``width``: effective kernel width (default = ``NaN``, indicating that the width is automatically chosen)
  - ``npoints``: the number of sample points


.. py:function:: kde(x, y[; width=NaN, resolution=25])

  Two-dimensional kernel density estimation. This function returns an instance of ``BivariateKDE`` defined as below:
  
  .. code-block:: julia

    immutable BivariateKDE
        x::Vector{Float64}
        y::Vector{Float64}
        density::Matrix{Float64}
    end

