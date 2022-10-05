# FuzzyCRegression.jl tutorial

## Installation

```julia
using Pkg; Pkg.add("FuzzyCRegression")
```
will install the package and its dependencies, which include [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) for minimization and [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) for automatic differentiation.

## Fitting the FCR model
There are two ways to fit an FCR model, depending on whether the data is stored as a DataFrame or as a set of arrays.

If the dataset is stored as a [DataFrame](https://dataframes.juliadata.org/stable/), the model can be fit using `fit(df,y,X,G,m,...)`, where the variables are referenced by their column names. `G` specifies the number of groups and `m` represents the regularization parameter, where group membership becomes less fuzzy as `m` approaches 1. For example, using the iris dataset from [RDatasets](https://github.com/JuliaStats/RDatasets.jl):

```julia
using FuzzyCRegression, RDatasets

iris = dataset("datasets", "iris")

fcr_model = fit(iris,y = "SepalWidth", x = ["1","SepalLength"], G = 3, m = 1.5)
summarize(fcr_model)
```
where "1" specifies a constant term. An advantage of this approach is that the estimated coefficients are labeled by variable name.

Alternatively, the data can be passed directly as arrays:

```julia
y = iris.SepalWidth
X = [ones(length(y)) iris.SepalLength]

fcr_model = fit(y, X, G = 3, m = 1.5)
summarize(fcr_model)
```

The arguments for fitting the model are:
  - `df`: name of dataframe (if missing, data must be passed as arrays)
  - `y`: column name or array holding values of the dependent variable
  - `X`: a list of column names or a matrix holding values of the independent variable(s) with heterogeneous coefficients (defaults to a constant)
  - `Z`: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients
  - `G`: number of groups
  - `m`: regularization parameter (greater than 1), where group assignment becomes binary as $m \rightarrow 1$
  - `unit`: column name or array with unit identifier (if panel structure)
  - `time`: column name or array with time indicators (if panel structure)
  - `startvals`: number of starting values for the minimization routine (default = 10)
  - `cores`: number of parallel workers (default = 1)

 ## Methods applied to fitted models
 
 The package provides several methods that can be applied to fitted models. The names are similar to those in the [GLM.jl](https://juliastats.org/GLM.jl/stable/) package. Full documentation for these methods can be found [here](https://aidantr.github.io/FuzzyCRegression.jl/dev/API/).
 
- `aic`: Akaike's Information Criterion
- `bic`: Bayesian Information Criterion
- `coef`: estimates of the coefficients in the model
- `confint`: confidence intervals for coefficients
- `distribution`: distribution of coefficients using group weights
- `predict`: obtain predicted values of the dependent variable from the fitted model, using modal group membership
- `residuals`: vector of residuals from the fitted model, using modal group membership
- `stderror`: standard errors of the coefficients
- `summarize`: summarize model results
- `vcov`: variance-covariance matrix of the coefficient estimates

For example, to plot the distribution of coefficients on SepalLength from the previous fitted model:

```julia
SepalWidth_coefs = distribution(fcr_model,"SepalWidth")

using Gadfly
plot(SepalWidth_coefs, Geom.hist)
```
![](assets/iris_plot.svg)

## Simple example 

To illustrate the package's functionality, we start with a simple example estimating grouped fixed effects with no controls. 

## More complicated example

We now turn to a more complicated example where the data has a panel structure, we estimate heterogeneous coefficients on the independent variable for each time period, and we include a common control variable. 

## Choosing $G$

## Choosing $m$ 
 
