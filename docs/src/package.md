# FuzzyCRegression Tutorial

## Installation

```julia
Pkg.add("FuzzyCRegression")
```
will install the package and its dependencies, which include [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) for minimization and [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) for automatic differentiation.

## Fitting the FCR model
There are two ways to fit an FCR model, using DataFrames or arrays. 

If the data is stored as a [DataFrame](https://dataframes.juliadata.org/stable/), the model can be fit using `fit(df,y,X,G,m,...)`, where the variables are referenced by their column names. For example, using the iris dataset from [RDatasets](https://github.com/JuliaStats/RDatasets.jl):

```julia
using FuzzyCRegression, RDatasets

iris = dataset("datasets", "iris")

fcr_model = fit(iris,y = "SepalWidth", x = ["1","SepalLength"], G=3, m=1.5)
summarize(fcr_model)
```
where "1" specifies a constant term. An advantage of this approach is that the estimated coefficients are labeled by variable name.

Alternatively, the data can be passed directly as arrays:

```julia
using FuzzyCRegression, RDatasets

iris = dataset("datasets", "iris")
y = iris.SepalWidth
X = [ones(length(y)) iris.SepalLength]

fcr_model = fit(y, X, G=3, m=1.5)
summarize(fcr_model)
```

The optional arguments for fitting the model are:
  - `Z`: a list of column names or a matrix holding values of the independent variable(s) with homogeneous coefficients
  - `G`: number of groups
  - `m`: regularization parameter (greater than 1), where group assignment becomes binary as $m \rightarrow 1$
  - `unit`: unit identifier (if panel structure)
  - `time`: vector of time indicators (if panel structure)
  - `startvals`: number of starting values for the minimization routine (default = 100)
  - `cores`: number of parallel workers (default = 1)

 ## Methods applied to fitted models
 
 The package provides several methods that can be applied to fitted models. The names are similar to those in the [GLM.jl](https://juliastats.org/GLM.jl/stable/) package.
 
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

For example, to plot the distribution of coefficients on SepalLength from the previous example:

```julia
using Gadfly

SepalLength_coefs = distribution(fcr_model,"SepalLength")
plot(SepalLength_coefs, geom.hist)

```
## Simple example 

To illustrate the package's functionality, we start with a simple example estimating grouped fixed effects with no controls. 

## More complicated example

We now turn to a more complicated example where the data has a panel structure, we estimate heterogeneous coefficients on the independent variable for each time period, and we include a common control variable. 

## Choosing $G$

## Choosing $m$ 
 
