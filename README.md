# FuzzyCRegression.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aidantr.github.io/FuzzyCRegression.jl/dev/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aidantr.github.io/FuzzyCRegression.jl/dev/)
[![Build Status](https://github.com/aidantr/FuzzyCRegression.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aidantr/FuzzyCRegression.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aidantr/FuzzyCRegression.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aidantr/FuzzyCRegression.jl)




| Documentation | CI Status | Coverage 
|:-----------------:|:------------------:|:-----------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][ci-img]][ci-url] | [![][codecov-img]][codecov-url] |
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://aidantr.github.io/FuzzyCRegression.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://aidantr.github.io/FuzzyCRegression.jl/dev/

[ci-img]: https://github.com/JuliaStats/FuzzyCRegression.jl/workflows/CI-stable/badge.svg
[ci-url]: https://github.com/aidantr/FuzzyCRegression.jl/actions/workflows/CI.yml/badge.svg?branch=main

[codecov-img]: https://codecov.io/gh/aidantr/FuzzyCRegression.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/aidantr/FuzzyCRegression.jl








This package implements the heterogeneous effects estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers in Julia. 

The [documentation](https://aidantr.github.io/FuzzyCRegression.jl/dev/) briefly describes the estimator and introduces the `FuzzyCRegression.jl` package through several examples.

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
