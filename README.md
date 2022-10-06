# FuzzyCRegression.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aidantr.github.io/FuzzyCRegression.jl/dev/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aidantr.github.io/FuzzyCRegression.jl/dev/)
[![Build Status](https://github.com/aidantr/FuzzyCRegression.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aidantr/FuzzyCRegression.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aidantr/FuzzyCRegression.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aidantr/FuzzyCRegression.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package implements the heterogeneous effects estimator from [Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)](https://drive.google.com/file/d/1U_MJHtJcB7H1Edv3xceilU_HJoxhLssP/view) in Julia. 

The [documentation](https://aidantr.github.io/FuzzyCRegression.jl/dev/) briefly introduces the estimator and describes the `FuzzyCRegression.jl` package.

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
