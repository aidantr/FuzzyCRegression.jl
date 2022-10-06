module FuzzyCRegression

using DataFrames
using Distributions
using ForwardDiff
using LinearAlgebra
using Optim
using PrettyTables
using Statistics
using StatsModels

include("fit.jl")
include("methods.jl")

export fit
export aic, bic
export coef, weights, coefnames
export stderror, vcov
export predict, residuals
export distribution
export summarize

end

