Abstraction for Statistical Models
===================================

This package defines an abstract type ``StatisticalModel``, and an abstract subtype ``RegressionModel``. 

Particularly, instances of ``StatisticalModel`` implement the following methods.  Note that the naming of the methods follow the convention of the `R language <http://www.r-project.org>`_.

- **coef** (model)

    Return the estimates of the coefficients, or, more generally, the parameters in the fitted model.
  
- **coeftable** (model)

    Return a `CoefTable` modelect for the fitted model with rows corresponding to parameters and columns usually including the estimates, their standard errors, a test statistic and a p-value, if appropriate.
  
- **confint** (model[, prob])

    A two-column matrix of confidence intervals on the parameters of the fitted model with coverage probability ``prob``, which defaults to `0.95`.
  
- **deviance** (model)

    The deviance of the fitted model at the parameter estimates.

- **fit** (Type{Model}, params...)

    Fit a model to data. If the model can be fit using a design matrix and response vector, the model should implement ``fit(Type{Model}, X, y, params...)``.
  
- **loglikelihood** (model)

    The log-likelihood of the fitted model at the parameter estimates.
  
- **stderr** (model)

    Standard errors of the parameters in the fitted model. For models based on a linear predictor these are for the coefficient parameters only.
  
- **vcov** (model)

    Estimated covariance matrix of the parameters in the fitted model.


``RegressionModel`` extends ``StatisticalModel``, which also implements the following methods

 
- **nobs** (model)

    The number of data values (observations) in the response to which the model was fit.
  
- **model_response** (model)

    The response vector to which the model was fit.
  
- **predict** (model)

    The fitted values from the model.
  
- **residuals** (model)

    The residuals at the parameter estimates for the fitted model.
