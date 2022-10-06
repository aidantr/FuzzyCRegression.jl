```@meta
CurrentModule = FuzzyCRegression
```

# FuzzyCRegression.jl manual

This package implements the heterogeneous effects estimator from [Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)](https://drive.google.com/file/d/1U_MJHtJcB7H1Edv3xceilU_HJoxhLssP/view) in Julia. 

## Fuzzy C-Regression estimator

Fuzzy C-Regression (FCR) is a method for estimating heterogeneous effects in settings with grouped patterns of unobserved heterogeneity. It extends the "Fuzzy C-Means" clustering algorithm to regression settings.

FCR can be used to estimate ["grouped fixed effects"](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA11319) (a constant term for each group) as well as heterogeneous coefficients.

In panel settings, each unit belongs to the same group over time, meaning that separate coefficients are estimated for each group-time pair.

#### FCR objective function
Consider a linear model with grouped heterogeneity:

$$y = \sum_{g=1}^G\mu_{g} \theta_{g}X+\varepsilon$$
where $X$ are covariates (which could simply be a constant), $\theta_{g}$ represent group-specific coefficients for groups $g=1,\ldots,G$, and $\mu_{g}$ represent group weights which sum to 1. FCR is concerned with jointly estimating $\theta$ and $\mu$ for each group.

The FCR objective function takes the form:

$$L^{FCR}_m\left(\theta,\mu\right)=\mathbb{E}\left[\sum_{g=1}^{G}\mu_{g}^{m}\left\Vert y-\theta_g X\right\Vert ^{2}\right]$$

where $m > 1$ is a regularization parameter governing the "fuzziness" of the FCR clusters (group membership becomes binary as $m \rightarrow 1$. The weights are defined as

$$\mu_{g}\left(y,X;\theta,m\right)=\left(\sum_{h=1}^{G}\frac{\left\Vert y-\theta_g X\right\Vert ^{2/\left(m-1\right)}}{\left\Vert y-\theta_h X\right\Vert ^{2/\left(m-1\right)}}\right)^{-1} \text{ for } g=1,\ldots,G$$


Combing these two equations, we can write the FCR objective as a continuous function of only the group-specific errors $||y-\theta_gX||$:

$$J^{FCR}\left(\theta\right)=\mathbb{E}\left[\left(\sum_{g=1}^{G}\left\Vert y-\theta_g X\right\Vert ^{-2/\left(m-1\right)}\right)^{1-m}\right]$$

Thus, for fixed $m$, the FCR function is differentiable and can be written as a standard GMM problem. 

#### Useful properties

  - __Fast:__ FCR can be solved in a single step through standard non-linear minimization. This makes it substantially faster than previous approaches, which require iteration over all possible groupings of units. While FCR remains quite efficient with large datasets, computation times for iterative algorithms become prohibitive.

  - __Customizable:__ The “fuzziness” of the FCR groups is governed by the regularization parameter $m$, where group assignment becomes binary as $m \rightarrow 1$. Choosing $m$ allows the estimator to better accommodate the uncertainty of group membership in realistic datasets, where noise means that cluster membership cannot be ascertained with certainty. Moreover, it means that FCR can recover the full distribution of effects.

  - __Inference:__ Since FCR is a GMM problem, it's asymptotic properties follow from standard theory and we are able to derive analytic standard errors.

## Package installation

FuzzyCRegressions.jl can be installed using:

```julia
Pkg.add("FuzzyCRegression")
```
This adds the latest version of the package and its dependencies, which include [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) for minimization and [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) for automatic differentiation. 

## Fitting the FCR model
There are two ways to fit an FCR model, depending on whether the data is stored as a DataFrame or as a set of arrays.

If the dataset is stored as a [DataFrame](https://dataframes.juliadata.org/stable/), the model can be fit using `fit(df,y,X,G,m,...)`, where the variables are referenced by their column names. `G` specifies the number of groups and `m` sets the regularization parameter, where group membership becomes binary as $m \rightarrow 1^+$ (tips for selecting these options are discussed below). For example, using the iris dataset from [RDatasets](https://github.com/JuliaStats/RDatasets.jl):

```julia
using FuzzyCRegression, RDatasets

iris = dataset("datasets", "iris")

fcr_model = fit(df = iris,y = "SepalLength", x = ["SepalWidth", "PetalWidth"], G = 3, m = 1.5)
summarize(fcr_model)
```
An advantage of this approach is that the estimated coefficients in the regression output are labeled by variable name.

Alternatively, the data can be passed directly as arrays:

```julia
fcr_model = fit(y = iris.SepalLength, X = [iris.SepalWidth iris.PetalWidth], G = 3, m = 1.5)
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
 
 The package provides several methods that can be applied to fitted models. The names are similar to those in [GLM.jl](https://juliastats.org/GLM.jl/stable/). Full documentation for these methods can be found [here](https://aidantr.github.io/FuzzyCRegression.jl/dev/API/).
 
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

For example, to plot the distribution of coefficients on SepalWidth from the fitted model above:

```julia
SepalWidth_coefs = distribution(fcr_model,"SepalWidth")

using Gadfly
plot(SepalWidth_coefs, Geom.hist)
```
![](assets/iris_plot.svg)

## Choosing $G$

## Choosing $m$ 




