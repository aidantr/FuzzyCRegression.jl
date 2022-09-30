# Julia Implementation

To install the FuzzyCRegression.jl Julia package, type:

```julia
julia> using Pkg
julia> Pkg.add("FuzzyCRegression")
```

The package is organized around the function [FuzzyCRegression.fit](https://aidantr.github.io/FuzzyCRegression.jl/dev/#FuzzyCRegression.fit-Tuple{}), which estimates the model. It returns the type `FCRModel` which stores the estimation results. The functions `distribution` and `information` can then calculate the distribution of coefficients or the information (e.g. BIC) used to select the optimal number of groups. 

