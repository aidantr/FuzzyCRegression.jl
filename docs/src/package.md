# FuzzyCRegression.jl Tutorial

## Installation

```julia
Pkg.add("FuzzyCRegression")
```

will install the package and its dependencies.

## Fitting the FCR model

An FCR model is fit using `fit(y,X,G,m,...)`. The arguments are
  - `y` a vector holding values of the dependent variable
  - `X` a matrix holding values of the independent variable(s) that will have heterogeneous coefficients (the default is a constant term, which will estimate a fixed effect for each group)
  - `Z` a matrix holding values of the independent variable(s) that will have homogeneous coefficient
  - `G` number of groups
  - `m` regularization parameter (greater than 1), where group assignment becomes binary as $m \rightarrow 1$



  
  
  
 ## Methods applied to fitted models
 
 The package provides several methods that can be applied to fitted models. The names are similar to those in the [GLM.jl](https://juliastats.org/GLM.jl/stable/) package.
 
- `adjr2`: adjusted R² for a linear model (an alias for `adjr²`)
- `aic`: Akaike's Information Criterion
- `aicc`: corrected Akaike's Information Criterion for small sample sizes
- `bic`: Bayesian Information Criterion
- `coef`: estimates of the coefficients in the model
- `confint`: confidence intervals for coefficients
- `cooksdistance`: [Cook's distance](https://en.wikipedia.org/wiki/Cook%27s_distance) for each observation
- `deviance`: measure of the model fit, weighted residual sum of squares for lm's
- `dispersion`: dispersion (or scale) parameter for a model's distribution
- `dof`: number of degrees of freedom consumed in the model
- `dof_residual`: degrees of freedom for residuals, when meaningful
- `fitted`: fitted values of the model
- `glm`: fit a generalized linear model (an alias for `fit(GeneralizedLinearModel, ...)`)
- `lm`: fit a linear model (an alias for `fit(LinearModel, ...)`)
- `loglikelihood`: log-likelihood of the model
- `modelmatrix`: design matrix
- `nobs`: number of rows, or sum of the weights when prior weights are specified
- `nulldeviance`: deviance of the model with all predictors removed
- `nullloglikelihood`: log-likelihood of the model with all predictors removed
- `predict`: predicted values of the dependent variable from the fitted model
- `r2`: R² of a linear model (an alias for `r²`)
- `residuals`: vector of residuals from the fitted model
- `response`: model response (a.k.a the dependent variable)
- `stderror`: standard errors of the coefficients
- `vcov`: variance-covariance matrix of the coefficient estimates
 
