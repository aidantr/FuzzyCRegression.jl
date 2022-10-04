# FuzzyCRegression.jl Tutorial

## Installation

```julia
Pkg.add("FuzzyCRegression")
```

will install the package and its dependencies, which include [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) for minimization and [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) for automatic differentiation.

## Fitting the FCR model
There are two ways to fit an FCR model, using DataFrames and using simple arrays. The model is fit on a DataFrame (from the DataFrames.jl package) using `fit(df,y_name,x_name,G,m,...)`

An FCR model is fit using `fit(y,X,G,m,...)`. The arguments are
  - `y`: a vector holding values of the dependent variable
  - `X`: a matrix holding values of the independent variable(s) that will have heterogeneous coefficients (the default is a constant term, which will estimate a fixed effect for each group)
  - `Z`: a matrix holding values of the independent variable(s) that will have homogeneous coefficients
  - `G`: number of groups
  - `m`: regularization parameter (greater than 1), where group assignment becomes binary as $m \rightarrow 1$
  - `startvals`: number of starting values for the minimization routine (default = 100)

 
 ## Methods applied to fitted models
 
 The package provides several methods that can be applied to fitted models. The names are similar to those in the [GLM.jl](https://juliastats.org/GLM.jl/stable/) package.
 
- `aic`: Akaike's Information Criterion
- `bic`: Bayesian Information Criterion
- `coef`: estimates of the coefficients in the model
- `confint`: confidence intervals for coefficients
- `predict`: obtain predicted values of the dependent variable from the fitted model, using modal group membership
- `residuals`: vector of residuals from the fitted model, using modal group membership
- `stderror`: standard errors of the coefficients
- `vcov`: variance-covariance matrix of the coefficient estimates

## Simple example 

To illustrate the package's functionality, we start with a simple example estimating grouped fixed effects with no controls. 

## More complicated example

We now turn to a more complicated example where the data has a panel structure, we estimate heterogeneous coefficients on the independent variable for each time period, and we include a common control variable. 

## Choosing $G$

## Choosing $m$ 
 
