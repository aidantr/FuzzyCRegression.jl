# FuzzyCRegression.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aidantr.github.io/FuzzyCRegression.jl/dev/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aidantr.github.io/FuzzyCRegression.jl/dev/)
[![Build Status](https://github.com/aidantr/FuzzyCRegression.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aidantr/FuzzyCRegression.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aidantr/FuzzyCRegression.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aidantr/FuzzyCRegression.jl)

This package implements the heterogeneous effects estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers in Julia. 

The [documentation](https://aidantr.github.io/FuzzyCRegression.jl/dev/) provides an introduction to the estimator, describes the FuzzyCRegression.jl package, and shows several examples.

## Installation 

```julia
julia> Pkg.add("FuzzyCRegression")
```

## Citation

If you use `FuzzyCRegression.jl` in your work, please cite the following:

```tex
@article{lewis2022,
  author  = {Lewis, Daniel and Davide Melcangi and Laura Pilossph and Aidan Toner-Rodgers},
  title   = {Approximating Grouped Fixed Effects Estimation via Fuzzy Clustering Regression},
  year    = {2022},
}
```
