Abstraction for Statistical Models
===================================

This package defines an abstract type ``StatisticalModel``, and an abstract subtype ``RegressionModel``. 

Particularly, instances of ``StatisticalModel`` implement the following methods.  Note that the naming of the methods follow the convention of the `R language <http://www.r-project.org>`_.

.. function:: coef(model)

    Return the estimates of the coefficients, or, more generally, the parameters in the fitted model.
  
.. function:: coeftable(model)

    Return a `CoefTable` modelect for the fitted model with rows corresponding to parameters and columns usually including the estimates, their standard errors, a test statistic and a p-value, if appropriate.
  
.. function:: confint(model[, prob])

    A two-column matrix of confidence intervals on the parameters of the fitted model with coverage probability ``prob``, which defaults to `0.95`.
  
.. function:: deviance(model)

    The deviance of the fitted model at the parameter estimates.

.. function:: fit(Type{Model}, params...)

    Fit a model to data. If the model can be fit using a design matrix and response vector, the model should implement ``fit(Type{Model}, X, y, params...)``.
  
.. function:: loglikelihood(model)

    The log-likelihood of the fitted model at the parameter estimates.
  
.. function:: stderr(model)

    Standard errors of the parameters in the fitted model. For models based on a linear predictor these are for the coefficient parameters only.
  
.. function:: vcov(model)

    Estimated covariance matrix of the parameters in the fitted model.


``RegressionModel`` extends ``StatisticalModel``, which also implements the following methods

 
.. function:: nobs(model)

    The number of data values (observations) in the response to which the model was fit.
  
.. function:: model_response(model)

    The response vector to which the model was fit.
  
.. function:: predict(model)

    The fitted values from the model.
  
.. function:: residuals(model)

    The residuals at the parameter estimates for the fitted model.
