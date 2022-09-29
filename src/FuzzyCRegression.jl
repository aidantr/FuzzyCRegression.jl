module FuzzyCRegression

using DataFrames, Statistics, StatsModels, LinearAlgebra, Optim, ForwardDiff

include("fit.jl")
include("methods.jl")

end
