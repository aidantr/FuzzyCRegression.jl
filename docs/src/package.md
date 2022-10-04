# FuzzyCRegression.jl Tutorial

## Installation

```julia
Pkg.add("FuzzyCRegression")
```
will install the package and its dependencies, which include [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) for minimization and [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) for automatic differentiation.

## Fitting the FCR model
There are two ways to fit an FCR model, using DataFrames or using simple arrays. 

If the dataset is stored in a [DataFrame](https://dataframes.juliadata.org/stable/), the model can be fit using `fit(df,y_name,x_names,G,m,...)`, where the independent and dependent variables are referenced using their column names. For example, using the iris dataset from the [RDatasets](https://github.com/JuliaStats/RDatasets.jl) package:

```julia
using FuzzyCRegression, RDatasets

iris = dataset("datasets", "iris")
fcr_model = fit(iris,y_names = "SepalWidth", x_names = c("Constant","SepalLength), G=3, m=1.5)
```

where "Constant" specifies a vector of 1's. An advantage of this approach is that the coefficients are labeled by variable name.

```julia
coef(fcr_model)
```

Alternatively, the data can be passed as arrays

```julia
using FuzzyCRegression, RDatasets

iris = dataset("datasets", "iris")
y = iris.SepalWidth
X = [ones(length(y)) iris.SepalLength]
fcr_model = fit(y, X, G=3, m=1.5)
```

The optional arguments for fitting a model are:
  - `Z`: a list of column names or matrix holding values of the independent variable(s) that will have homogeneous coefficients
  - `G`: number of groups
  - `m`: regularization parameter (greater than 1), where group assignment becomes binary as $m \rightarrow 1$
  - `startvals`: number of starting values for the minimization routine (default = 100)
  - `cores`: number of parallel workers (default = 1)

 
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
 
