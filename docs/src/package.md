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
  - 
  
  
  
 ## Methods applied to fitted models
 
 The package provides a number of methods that can be applied to fitted models:
 
